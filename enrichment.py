#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import json
import multiprocessing
import hashlib
import time

from multiprocessing import Process
from collections import defaultdict
from ConfigParser import ConfigParser

from statsmodels.stats import multitest
import MySQLdb as mysql

from __init__ import fetch
import enrichment_analysis as ea
from Annotations import parse_flat


config = ConfigParser()
config.read('settings.cfg')

ANNO_MAX_SIZE = config.getint('Annotations', 'max size')
ANNO_MIN_SIZE = config.getint('Annotations', 'min size')

QVAL_CUTOFF = config.getfloat('FDR', 'cutoff')
NCORES = multiprocessing.cpu_count()


# MySQL settings
MYUSER = config.get('MySQL', 'user')
MYHOST = config.get('MySQL', 'host')
MYPASS = config.get('MySQL', 'pass')
MYDB = config.get('MySQL', 'db')

# MySQL commands (reference results_db_schema.sql)
store_results_sql = """
insert or replace into results (_id, goid, term, pval, dataset, factor, subset, year, num_genes)
 values (:id, :goid, :term, :pval, :dataset, :factor, :subset, :year, :numgenes)
"""

select_results_sql = """
select _id, pval from results where dataset=? and subset = ? and year=? order by pval
"""

insert_postinfo_sql = """
update results set qval = :qval, num_annos= :numannos, anno_max= :annomax, anno_min = :annomin where _id = :id
"""

mapfile = 'data/uniprot2entrez.json'

def split(annotation_dict, blocks=8):
    returnlist = [dict() for b in xrange(blocks)]
    b = 0
    for k, v in annotation_dict.iteritems():
        returnlist[b][k] = v
        b = b + 1 if b < blocks - 1 else 0
    return returnlist


def filter_annos(annotation_dict, _max=ANNO_MAX_SIZE, _min=ANNO_MIN_SIZE):
    """Removes annotations that are composed of more than 'max' genes
    or fewer than 'min'.
    This returns a filtered copy."""
    
    returndict = dict(annotation_dict)
    for k, v in annotation_dict.iteritems():
        if len(v['genes']) > _max or len(v['genes']) < _min:
            del returndict[k]
    print "Removed %d terms from annotation set." % (len(annotation_dict) - len(returndict))
    return returndict


def restrict_subontology(annotation_dict, ontology, year):
    go_id, name = {'MF': ('GO:0003674', 'Molecular Function'),
                   'CC': ('GO:0005575', 'Cellular Component'),
                   'BP': ('GO:0008150', 'Biological Process')}[ontology]
    try:
        ontology = parse_flat("data/go-%s.flat" % year)
    except IOError:
        print("Warning: Flattened ontology file not found for year: %s. Create flattened ontology using go_flattener.jar." % year)
        print("No subontology restriction done; using all terms from all ontologies.")
        return annotation_dict
    returndict = dict(annotation_dict)
    for term in annotation_dict:
        if not go_id in ontology[term]:
            del returndict[term]
    print "Restricting to %s removed %d terms." % (name, (len(annotation_dict) - len(returndict)))
    return returndict
            

def store_in_db(fn):
    def store(dataset, platform, factor, subset, annotations, year, u2emap):
        results, diffexp = fn(dataset, platform, factor, subset, annotations, year, u2emap)
        p = multiprocessing.current_process()
        with mysql.connect(MYHOST, MYUSER, MYPASS, MYDB) as conn:
            conn.executemany(store_results_sql, ({
                        # the id field is there to prevent redundant entries
                        'id': hashlib.md5('|'.join([t,dataset.id,factor,subset,year])).hexdigest(),
                        'goid': t,
                        'term': annotations[t]['name'],
                        'pval': pval,
                        'dataset': dataset.id,
                        'factor': factor,
                        'subset': subset,
                        'year': year,
                        'numgenes': len(diffexp)} for t, pval in results.iteritems()))
            conn.commit()
        print("<%s> DONE: Stored %d terms in db" % (p.name, len(results)))
    return store


@store_in_db
def enriched(dataset, platform, factor, subset, annotations, 
                year, uniprot2entrez_map):
    diffexp = dataset.diffexpressed(subset, factor,
        QVAL_CUTOFF)
    diffexp = ea.map2entrez(platform, probes=diffexp)
    background = ea.map2entrez(platform)
    not_diffexp = [x for x in background if x not in diffexp]
    p = multiprocessing.current_process()
    total = len(annotations)
    if len(diffexp) == 0:
        print("Warning: no differentially expressed genes found for " +
            "%s:%s" % (factor, subset))

    results = {}
    for i, term in enumerate(annotations):
        pval = ea._fexact(diffexp, not_diffexp, background, annotations[term],
                             uniprot2entrez_map)
        if pval < PVAL_CUTOFF:
            print "<{name}>: ({i}/{total}) {pval}\t{term}".format(name=p.name,
                i=i, pval=pval, term=term, total=total)
        results[term] = pval
    return results, diffexp


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python enrichment-hpc.py <GDS file or accn> <{MF, CC, BP}> <annotation file 1> [more anno files...]")
        print("\n See README.md for more information.")
        sys.exit(1)
    file_or_accn = sys.argv[1]
    ontology = sys.argv[2]
    annotation_files = sys.argv[3:]
    
    # this file can be downloaded from Uniprot's mapping service
    uniprot2entrez_map = json.load(open(mapfile))
    assert len(uniprot2entrez_map) > 27000

    # import the dataset
    dataset = fetch(file_or_accn, destdir='data')
    dataset = dataset.to_numeric()
    dataset.filter().log2xform()

    # import the annotation files (in JSON format)
    annotation_years = (json.load(open(f)) for f in annotation_files)

    # acquire the platform used from the dataset metadata
    platform = fetch(dataset.meta['platform'], destdir='data')

    print("Detected %d cores, splitting into %d subprocesses..." % (NCORES, NCORES-1))

    for annotations in annotation_years:
        year = annotations['meta']['year']
        print "Filtering out annotation gene sets greater than %d and less than %d" % (ANNO_MAX_SIZE, ANNO_MIN_SIZE)
        filtered_annotations = filter_annos(annotations['anno'])
        filtered_annotations = restrict_subontology(filtered_annotations, ontology, year)
        blocks = split(filtered_annotations, blocks=NCORES-1)
        print("Split %d annotations into %d blocks of ~%d terms each..." % (len(filtered_annotations), len(blocks), len(blocks[0])))
        # We're actually going to only look at "disease state" factors for this analysis.
        factor = 'disease state'
        for subset in dataset.factors[factor]:
            print("-- [year: %s] [dataset: %s] [%s: %s] --" % (year, dataset.id, 
                                                               factor, subset))
            jobs = []
            for block in blocks:
                p = Process(target=enriched,
                            args=(dataset, platform, factor, subset, 
                                  block, year, uniprot2entrez_map))
                jobs.append(p)
                p.start()
            [p.join() for p in jobs] # wait for them all to finish

        # Calculate the q-values now that we have all the p-values, and also add the number of terms
        for subset in dataset.factors[factor]:
            with mysql.connect(MYHOST, MYUSER, MYPASS, MYDB) as conn:
                cursor = conn.execute(select_results_sql, (dataset.id,subset,year))
                results = list(cursor.fetchall())    # list of tuples [(id, pval), ...]
                pvals = [x[1] for x in results] 
                ids = [x[0] for x in results]
                rejected, qvals = multitest.fdrcorrection(pvals)
                results = zip(ids,qvals)
                conn.executemany(insert_postinfo_sql, ({'qval': qval,
                                                        'id': id,
                                                        'numannos': len(filtered_annotations),
                                                        'annomax': ANNO_MAX_SIZE,
                                                        'annomin': ANNO_MIN_SIZE
                                                        } for id, qval in results))
                conn.commit()
                print("Main: Inserted q-values, num. terms for annotation year %s" % year)

                                     
                    
                    
                                 
                        

                    
