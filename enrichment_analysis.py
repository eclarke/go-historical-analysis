#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __init__ import fetch
from numpy import array
from numpy import intersect1d
from collections import defaultdict
import json
import scipy.stats as stats
import urllib
import urllib2

import os

PVAL_CUTOFF = 0.05
DIFF_AVG_CUTOFF = 1
UNIPROT = 'http://www.uniprot.org/mapping/'
mapfile = 'data/uniprot2entrez.json'


def main(file_or_accn, annotation_files):
    uniprot2entrez_map = json.load(open(mapfile))
    assert len(uniprot2entrez_map) > 27000
    # import the dataset
    dataset = fetch(file_or_accn, destdir='/tmp/eanalysis/')
    dataset = dataset.to_numeric()
    dataset.filter().log2xform()

    # acquire the platform used from the dataset metadata
    platform = fetch(dataset.meta['platform'], destdir='/tmp/eanalysis/')

    # import the annotation files (in JSON format)
    annotations = [json.load(open(f)) for f in annotation_files]

    # autovivified tree
    def tree():
        return defaultdict(tree)

    enriched_terms = tree()

    for factor in dataset.factors:
        for subset in dataset.factors[factor]:
            diffexp = dataset.diffexpressed(subset, factor,
                PVAL_CUTOFF, DIFF_AVG_CUTOFF)

            # convert the probe names to entrez
            diffexp = map2entrez(platform, probes=diffexp)  # [x.decode('ascii') for x in map2entrez(platform, probes=diffexp)]
            background = map2entrez(platform)   # [x.decode('ascii') for x in map2entrez(platform)]
            not_diffexp = [gene for gene in background if gene not in diffexp]
            if len(diffexp) == 0:
                print("Warning: no differentially expressed genes found for " +
                    "{}:{}".format(factor, subset))

            # run an enrichment analysis on the enriched genes for each
            # annotation set
            for annotation in annotations:
                year = annotation['meta']['year']
                print "({}) {}: {}".format(year, factor,
                    subset)
                try:
#                    for term in annotation['anno']:
                    for term in progress.bar(annotation['anno']):
                        pval = _fexact(diffexp, not_diffexp, background,
                            annotation['anno'][term], uniprot2entrez_map)
                        enriched_terms[factor][subset][year][term] = pval
                except KeyboardInterrupt:
                    with open('.intermediate_results.json', 'wb') as out:
                        json.dump(enriched_terms, out)
                        print "\n--break: enriched terms dumped to ", out.name

                results = enriched_terms[factor][subset][year]
                top = sorted([(results[term], term) for term in results],
                    key=lambda x: x[0])[0:10]
                if 1.0 in [r[0] for r in top]:
                    print "No enriched terms found.\n"
                else:
                    for r in top:
                        print "{}\t{}".format(r[0],
                            annotation['anno'][r[1]]['name'])
                    print
    return enriched_terms


def _fexact(diffexp, not_diffexp, background, term, uniprot2entrez_map):
    """
    Conducts a modified Fisher's exact test (aka EASE score) on the differentially 
    expressed genes and the provided gene list for the term. The test is modified
    to be more conservative in that it subtracts one from the 'intersection' cell
    of the contingency table. For rationale, see PMCID: PMC328459. Near-identical 
    to DAVID's Functional Annotation Tool.

    Returns the p-value (the likelihood of that term being enriched against a null
    distribution) for the term's gene set.

    Arguments:
    diffexp: a list of differentially expressed genes (in Entrez Gene id format)
    not_diffexp: a list of genes *not* differentially expressed in the background
    background: all genes (often all genes tested by the probe set)
    term: a dict of the form {'name':'term name', 'genes':['1234','1235',...]}
    """

    if not diffexp:
        return 1.0

    # convert Uniprot ids to Entrez Gene ids
    term_genes = map_uniprot(term['genes'], uniprot2entrez_map)
    
    

    # contingency table for fisher's exact test:
    # |  g_e - 1  | g_ne  |  (e = diff. expressed, ne = not diff. expressed,
    # +–––––-----–+–-–––––+   g = geneset,  ng = not geneset)
    # |    ng_e   | ng_ne |

    g_e = len(intersect1d(term_genes, diffexp)) - 1
    if g_e < 1:
        return 1.0
    g_ne = len(intersect1d(term_genes, not_diffexp))
    not_term_genes = [gene for gene in background if gene not in term_genes]
    ng_e = len(intersect1d(not_term_genes, diffexp))
    ng_ne = len(intersect1d(not_term_genes, not_diffexp))
    table = array([[g_e, g_ne], [ng_e, ng_ne]])
    
    odds, pval = stats.fisher_exact(table)
    
    return pval


def map_uniprot(uniprots, uniprot2entrez_map):
    if uniprot2entrez_map:
        return [uniprot2entrez_map[x] for x in uniprots if x in uniprot2entrez_map]
    else:
        print "fetching from web... :<"
        params = {
            'from': 'ACC+ID',
            'to': 'P_ENTREZGENEID',
            'format': 'list',
            'query': ' '.join(uniprots)
        }
        data = urllib.urlencode(params)
        request = urllib2.Request(UNIPROT, data)
        response = urllib2.urlopen(request).read().split('\n')
        return response


def map2entrez(platform, probes=None):
    """Returns the Entrez Gene ID for the probes. Some probes may not map (i.e,
        controls, etc)- these are simply empty strings in the returned array.
    """
    try:
        entrez_column = platform.table[0].index('ENTREZ_GENE_ID')
        if probes == None:
            return [x[entrez_column] for x in platform.table]
        else:
            return [x[entrez_column] for x in platform.table if x[0] in probes]
    except ValueError:
        print "Error: could not find Entrez mappings for this platform."
        print probes
        raise


def test(file_or_accn, annotation_file):
    try:
        uniprot2entrez_map = json.load(open(mapfile))
        assert len(uniprot2entrez_map) > 27000
        # import the dataset
        dataset = fetch(file_or_accn)
        dataset = dataset.to_numeric()
        dataset.filter().log2xform()

        # acquire the platform used from the dataset metadata
        if str(dataset.meta['platform']) == 'GPL570':
            # ugly hack, remove
            platform = fetch('/var/folders/xc/jvc059n56y91z1n9_48flnmm0000gp/T/tmpVjr6zS/GPL570.soft')
        else:
            platform = fetch(dataset.meta['platform'])

        # import the annotation files (in JSON format)
        annotation = json.load(open(annotation_file))

        # autovivified tree
        def tree():
            return defaultdict(tree)

        enriched_terms = tree()

        factor = 'disease state'
        subset = 'glioblastomas'
        diffexp = dataset.diffexpressed(subset, factor,
            PVAL_CUTOFF, DIFF_AVG_CUTOFF)

        # convert the probe names to entrez
        diffexp = map2entrez(platform, probes=diffexp)
        background = map2entrez(platform)
        not_diffexp = [gene for gene in background if gene not in diffexp]
        if len(diffexp) == 0:
            print("Warning: no differentially expressed genes found for " +
                "{}:{}".format(factor, subset))

        # run an enrichment analysis on the enriched genes for each
        # annotation set

        year = annotation['meta']['year']
        print "({}) {}: {}".format(year, factor,
            subset)
        try:
            for term in progress.bar(annotation['anno']):
                pval = _fexact(diffexp, not_diffexp, background,
                    annotation['anno'][term], uniprot2entrez_map)
                enriched_terms[factor][subset][year][term] = pval
        except KeyboardInterrupt:
            with open('.intermediate_results.json', 'wb') as out:
                json.dump(enriched_terms, out)
                print "\n--break: enriched terms dumped to ", out.name

        results = enriched_terms[factor][subset][year]
        top = sorted([(results[term], term) for term in results],
            key=lambda x: x[0])[0:10]
        if 1.0 in [r[0] for r in top]:
            print "No enriched terms found.\n"
        else:
            for r in top:
                print "{}\t{}".format(r[0],
                    annotation['anno'][r[1]]['name'])
            print
        return enriched_terms
    except KeyboardInterrupt:
        return

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print """
Usage: enrichment.py <dataset filename or accession> <annotation JSON 1>
        <annotation JSON 2> <annotation JSON ...> ..."""
    else:
        results = main(sys.argv[1], sys.argv[2:])
        import pprint
        pprint.pprint(results)
