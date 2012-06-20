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
from contextlib import closing
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
set_isolation_level = """set session transaction isolation level read committed"""

store_results_sql = """
replace into results (_id, _subid, ontology, goid, term, pval, dataset, factor, subset, year, num_annos, num_genes, anno_min, anno_max)
 values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
"""

select_pvals_sql = """
select _subid, pval from results where _id=%s
"""

insert_qval_sql = """
update results set qval=%s where _subid=%s
"""

mapfile = 'data/uniprot2entrez.json'


def split(iterable, blocks=8):
    if isinstance(iterable, dict):
        returnlist = [dict() for b in xrange(blocks)]
        b = 0
        for k, v in iterable.iteritems():
            returnlist[b][k] = v
            b = b + 1 if b < blocks - 1 else 0
        return returnlist
    else:
        returnlist = [list() for b in xrange(blocks)]
        b = 0
        for i in iterable:
            returnlist[b].append(i)
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


def get_connection(max_retries=30):
    i = 0
    while True and i <= max:
        try:
            db = mysql.connect(MYHOST, MYUSER, MYPASS, MYDB)
            with closing(db.cursor()) as c:
                c.execute(set_isolation_level)
                db.commit()
            return db
        except mysql.OperationalError:
            print("Operational error, sleeping for 5 seconds...")
            time.sleep(5)
            i += 1
    print("Max retries reached, aborting...")
    raise mysql.OperationalError("Failed to connect after %d retries" % max_retries)


def md5hash(*args):
    return hashlib.md5(''.join(args)).hexdigest()


def store_in_db(fn):
    def store(dataset, platform, factor, subset, annotations, year, num_annos, ontology, u2emap):
        results, diffexp = fn(dataset, platform, factor, subset, annotations, year, u2emap)
        p = multiprocessing.current_process()
        db = get_connection(100)
        with closing(db.cursor()) as c:
            col_results = []
            for goid, pval in results.iteritems():
                _id = md5hash(dataset.id, factor, subset, year, ontology)
                _subid = md5hash(goid, _id)
                col_results.append((_id, _subid, ontology, goid, annotations[goid]['name'], pval, dataset.id,
                                    factor, subset, year, num_annos, len(diffexp), ANNO_MIN_SIZE, ANNO_MAX_SIZE))
            c.executemany(store_results_sql, col_results)
            db.commit()
        db.close()
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
        if len(diffexp) == 0:
            pval = 1
        else:
            pval = ea._fexact(diffexp, not_diffexp, background, annotations[term],
                              uniprot2entrez_map)
        if pval < QVAL_CUTOFF:
            print "<{name}>: ({i}/{total}) {pval}\t{term}".format(name=p.name,
                                                                  i=i, pval=pval, 
                                                                  term=term, 
                                                                  total=total)
        results[term] = pval
    return results, diffexp


def multitest_correction(dataset, ontology, annotation_files):
    annotation_years = (json.load(open(f)) for f in annotation_files)
    factor = 'disease state'
    db = get_connection(100)
    for annotations in annotation_years:
        year = annotations['meta']['year']
        for subset in dataset.factors[factor]:
            print "[%s]-[%s]-[%s]-[%s]:" % (dataset.id, year, ontology, subset),
            with closing(db.cursor()) as c:
                _id = md5hash(dataset.id, factor, subset, year, ontology)
                print "selecting pvals... ",
                c.execute(select_pvals_sql, _id)
                results = list(c.fetchall())    # list of tuples [(_subid, pval), ...]
            pvals = [x[1] for x in results] 
            subids = [x[0] for x in results]
            print "calculating FDR... ",
            rejected, qvals = multitest.fdrcorrection(pvals)
            results = zip(qvals, subids)
            with closing(db.cursor()) as c:
                print "inserting %d qvals... " % len(results),
                c.executemany(insert_qval_sql, results)
                db.commit()
            print "done."
    db.close()


def _usage():
    print("Inserts into mysql database specified in settings.cfg results from the enrichment analysis of the given dataset "
          "against the specified annotation files, using one of the sub-ontologies specified.\n")
    print("Usage: python fdr_correction.py <GDS file or accn> <[MF, CC, BP]> <anno file 1> [more anno files...]")
    print("\n See README.md for more information.")
    sys.exit(1)


if __name__ == '__main__':
    if len(sys.argv) < 3 or len(sys.argv[2]) != 2:
        _usage()

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

    print("Detected %d cores, splitting into %d subprocesses..." % (NCORES, NCORES))

    for annotations in annotation_years:
        year = annotations['meta']['year']
        print "Filtering out annotation gene sets greater than %d and less than %d" % (ANNO_MAX_SIZE, ANNO_MIN_SIZE)
        filtered_annotations = filter_annos(annotations['anno'])
        filtered_annotations = restrict_subontology(filtered_annotations, ontology, year)
        blocks = split(filtered_annotations, blocks=NCORES)
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
                                  block, year, len(filtered_annotations), ontology, uniprot2entrez_map))
                jobs.append(p)
                p.start()
            [p.join() for p in jobs] # wait for them all to finish



        # Calculate the q-values now that we have all the p-values, and also add the number of terms
        """ This should be in a separate process that runs *after* we have all the pvals 
        for subset in dataset.factors[factor]:
            db = get_connection(100)
            with closing(db.cursor()) as c:
                c.execute(select_results_sql, (dataset.id, subset, int(year)))
                results = list(c.fetchall())    # list of tuples [(id, pval), ...]
                
            pvals = [x[1] for x in results] 
            ids = [x[0] for x in results]
            rejected, qvals = multitest.fdrcorrection(pvals)
            results = zip(ids,qvals)
            postinfo = [(qval, ontology, len(filtered_annotations), ANNO_MAX_SIZE, ANNO_MIN_SIZE, id)
                        for id, qval in results]
            print("Main: Inserting (or waiting to insert) additional information into results for subset...")
            for postinfo_chunk in split(postinfo, 4):
                with closing(db.cursor()) as c:
                    c.executemany(insert_postinfo_sql, postinfo)
                    db.commit()
            db.close()

        print("Main: Updated results for annotation year %s" % year)
        """
                                     
                    
                    
                                 
                        

                    
