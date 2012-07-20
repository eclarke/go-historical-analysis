#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __init__ import fetch
from numpy import array
from collections import defaultdict
import json
import scipy.stats as stats
import urllib
import urllib2

import os

UNIPROT = 'http://www.uniprot.org/mapping/'
mapfile = 'data/uniprot2entrez.json'

def _fexact(diffexp, not_diffexp, background, term, uniprot2entrez_map, EASE=True):
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
    EASE: if false, do a traditional Fisher's exact test, not the EASE modification
    """

    if not diffexp:
        return 1.0

    # convert Uniprot ids to Entrez Gene ids
    term_genes = map_uniprot(term['genes'], uniprot2entrez_map)

    # convert all to sets
    term_genes = set(term_genes)
    diffexp = set(diffexp)
    not_diffexp = set(not_diffexp)
    background = set(background)
    
    # contingency table for fisher's exact test:
    # |  g_e - 1  | g_ne  |  (e = diff. expressed, ne = not diff. expressed,
    # +–––––-----–+–-–––––+   g = geneset,  ng = not geneset)
    # |    ng_e   | ng_ne |

    g_e = len(term_genes.intersection(diffexp))
    g_e = g_e - 1 if EASE else g_e

    if g_e < 1:
        return 1.0

    g_ne = len(term_genes.intersection(not_diffexp)) 
    not_term_genes = set([gene for gene in background if gene not in term_genes])
    ng_e = len(not_term_genes.intersection(diffexp))
    ng_ne = len(not_term_genes.intersection(not_diffexp))
    table = array([[g_e, g_ne], [ng_e, ng_ne]])

    odds, pval = stats.fisher_exact(table, alternative='greater')
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
        controls, etc)- these are removed.
    """
    if 'ENTREZ_GENE_ID' in platform.table[0]:
        entrez_column = platform.table[0].index('ENTREZ_GENE_ID')
    elif 'GENE' in platform.table[0]:
        entrez_column = platform.table[0].index('GENE')
    else:
        raise ValueError('Cannot find Entrez mappings for this platform!')
    if probes == None:
        result = [x[entrez_column] for x in platform.table]
    else:
        result = [x[entrez_column] for x in platform.table if x[0] in probes]

    return [x for x in result if x and '/' not in x]
