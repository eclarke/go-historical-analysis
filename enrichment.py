#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sqlite3
import sys
import os
import json
import multiprocessing
import hashlib
import time

from multiprocessing import Process
from collections import defaultdict

from statsmodels.stats import multitest

from __init__ import fetch
import enrichment_analysis as ea
from Annotations import parse_flat

ANNO_MAX_SIZE=1000
ANNO_MIN_SIZE=3

PVAL_CUTOFF = 0.05
DIFF_AVG_CUTOFF = 1
NCORES = multiprocessing.cpu_count()

db_filename = 'results.db'
mapfile = 'data/uniprot2entrez.json'

# SQLite commands (reference results_db_schema.sql)
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
        while True:
            try:
                with sqlite3.connect(db_filename, timeout=30) as conn:
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
                break
            except sqlite3.OperationalError as e:
                print("<%s> OperationalError: %s. Sleeping for 2 seconds and retrying." % (p.name, e))
                time.sleep(2)
        print("<%s> DONE: Stored %d terms in db" % (p.name, len(results)))
    return store


@store_in_db
def enriched(dataset, platform, factor, subset, annotations, 
                year, uniprot2entrez_map):
    diffexp = dataset.diffexpressed(subset, factor,
        PVAL_CUTOFF, DIFF_AVG_CUTOFF)
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
    
    # check that database exists and has results table 
    try:
        assert os.path.isfile(db_filename)
        with sqlite3.connect(db_filename) as conn:
            conn.execute('select * from results limit 1')
    except (AssertionError, sqlite3.OperationalError):
        print "Fatal: SQLite database '%s' not configured correctly. Run `python init_results_db.py force` to create database with correct schema." % db_filename

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
            try:
                with sqlite3.connect(db_filename, timeout=30) as conn:
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
            except sqlite3.OperationalError as e:
                print("Main: OperationalError: %s. Sleeping for 2 seconds and retrying." % e)
                time.sleep(2)

                                     
                    
                    
                                 
                        

                    
