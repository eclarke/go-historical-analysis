#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sqlite3
import sys
import os
import json
import multiprocessing
from multiprocessing import Process
from __init__ import fetch
import enrichment_analysis as ea
import hashlib
import time


ANNO_MAX_SIZE=2000
ANNO_MIN_SIZE=10

PVAL_CUTOFF = 0.05
DIFF_AVG_CUTOFF = 1
NCORES = 8

db_filename = 'results.db'
mapfile = 'data/uniprot2entrez.json'
store_results_sql = """
insert into results (_id, goid, term, pval, dataset, factor, subset, year)
 values (:id, :goid, :term, :pval, :dataset, :factor, :subset, :year)
"""


def split(annotation_dict, blocks=8):
    returnlist = [dict() for b in xrange(blocks)]
    b = 0
    for k, v in annotation_dict.iteritems():
        returnlist[b][k] = v
        b = b + 1 if b < blocks - 1 else 0
    return returnlist


def filter_size(annotation_dict, _max=ANNO_MAX_SIZE, _min=ANNO_MIN_SIZE):
    """Removes annotations that are composed of more than 'max' genes
    or fewer than 'min'. This returns a filtered copy."""
    returndict = dict(annotation_dict)
    for k, v in annotation_dict.iteritems():
        if len(v['genes']) > _max or len(v['genes']) < _min:
            del returndict[k]
    print "Removed %d terms from annotation set." % (len(annotation_dict) - len(returndict))
    return returndict


def store_in_db(fn):
    def store(dataset, platform, factor, subset, annotations, year, u2emap):
        results = fn(dataset, platform, factor, subset, annotations, year, u2emap)
        p = multiprocessing.current_process()
        while True:
            try:
                with sqlite3.connect(db_filename) as conn:
                    conn.executemany(store_results_sql, ({
                                # the id field is there to prevent redundant entries
                                'id': hashlib.md5('|'.join([t,dataset.id,factor,subset,year])).hexdigest(),
                                'goid': t,
                                'term': annotations[t]['name'],
                                'pval': pval,
                                'dataset': dataset.id,
                                'factor': factor,
                                'subset': subset,
                                'year': year} for t, pval in results.iteritems()))
                    conn.commit()
                break
            except sqlite3.DatabaseError as e:
                print("<%s> DatabaseError: %s. Sleeping for 2 seconds and retrying." % e.message)
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
    return results


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python enrichment-hpc.py <GDS file or accn> <annotation file 1> [more anno files...]")
        print("\n See README.md for more information.")
        sys.exit(1)
    file_or_accn = sys.argv[1]
    annotation_files = sys.argv[2:]
    
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

    # acquire the platform used from the dataset metadata
    platform = fetch(dataset.meta['platform'], destdir='data')

    # import the annotation files (in JSON format)
    annotation_years = [json.load(open(f)) for f in annotation_files]

    for factor in dataset.factors:
        for subset in dataset.factors[factor]:
            for annotations in annotation_years:
                year = annotations['meta']['year']
                filtered_annotations = filter_size(annotations['anno'])
                blocks = split(filtered_annotations, blocks=NCORES)
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
    

