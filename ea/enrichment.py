#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import multiprocessing
import time

from multiprocessing import Process
from ConfigParser import ConfigParser
from optparse import OptionParser
from contextlib import closing
from statsmodels.stats import multitest
import MySQLdb as mysql

from __init__ import fetch
import enrichment_analysis as ea
from Annotations import parse_flat


# MySQL commands (reference results_db_schema.sql)
set_isolation_level = """set session transaction isolation level 
read committed"""


store_results_sql = """
replace into {table} (ontology, goid, term, pval, dataset, factor, subset, 
    year, num_annos, num_genes, anno_min, anno_max, min_depth, max_depth, 
    min_var, filter_similar, filter_size, filter_depth, shuffled)
 values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
"""

select_pvals_sql = """
select goid, pval from {table} where dataset=%s and subset=%s and year=%s 
and ontology=%s and shuffled=%s
"""

insert_qval_sql = """
update {table} set qval=%s where dataset=%s and subset=%s and year=%s
 and ontology=%s and shuffled=%s and goid=%s
"""


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


def filter_annos(annotations, _max, _min):
    """Removes annotations that are composed of more than 'max' genes
    or fewer than 'min'.
    This returns a filtered copy."""
    
    returned = dict(annotations)
    for k, v in annotations.iteritems():
        if len(v['genes']) > _max or len(v['genes']) < _min:
            del returned[k]
    print "Removed %d terms from annotation set." % (len(annotations) - 
        len(returned))
    return returned


def filter_annos_by_depth(annotations, min_depth, max_depth):
    returned = dict(annotations)
    for k, v in annotations.iteritems():
        if len(v['parents']) < min_depth or len(v['parents']) > max_depth:
            del returned[k]
    print "Removed %d terms from annotation set." % (len(annotations) -
        len(returned))
    return returned


def filter_similar_terms(annotations, min_variance):
    returned = dict(annotations)
    for k, v in annotations.iteritems():
        parents = v['parents']
        genes = v['genes']
        for p in parents:
            if p in annotations:
                pgenes = annotations[p]['genes']
                if (len(pgenes) - len(genes)) < min_variance:
                    if k in returned:
                        del returned[k]
                    continue
    print "Removed %d terms from annotation set." % (len(annotations) -
        len(returned))
    return returned


def restrict_subontology(annotation_dict, ontology, year):
    go_id, name = {'MF': ('GO:0003674', 'Molecular Function'),
                   'CC': ('GO:0005575', 'Cellular Component'),
                   'BP': ('GO:0008150', 'Biological Process')}[ontology]
    try:
        ontology = parse_flat("data/go-%s.flat" % year)
    except IOError:
        print("Warning: Flattened ontology file not found for year: %s."
            " Create flattened ontology using go_flattener.jar." % year)
        print("No subontology restriction done; using all terms from all "
            "ontologies.")
        return annotation_dict
    returndict = dict(annotation_dict)
    for term in annotation_dict:
        if not go_id in ontology[term]:
            del returndict[term]
    print "Restricting to %s removed %d terms." % (name, 
        (len(annotation_dict) - len(returndict)))
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
    raise mysql.OperationalError("Failed to connect after %d retries" 
        % max_retries)


def store_in_db(fn):
    def store(dataset, platform, factor, subset, annotations,
     year, shuffled, num_annos, ontology, u2emap):
        results, diffexp = fn(dataset, platform, factor, subset, annotations,
         year, u2emap)
        p = multiprocessing.current_process()
        db = get_connection(100)
        with closing(db.cursor()) as c:
            col_results = []
            for goid, pval in results.iteritems():
                if pval == 1:
                    continue    # we don't need to store pvals of 1
                col_results.append((ontology, goid, annotations[goid]['name'], 
                    pval, dataset.id, factor, subset, year, num_annos,
                    len(diffexp), ANNO_MIN_SIZE, ANNO_MAX_SIZE, MIN_DEPTH,
                    MAX_DEPTH, MIN_VARIANCE, FILTER_SIMILAR, FILTER_BY_SIZE,
                    FILTER_BY_DEPTH, shuffled))
            assert '{table}' not in store_results_sql
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
            pval = ea._fexact(diffexp, not_diffexp, background,
                annotations[term], uniprot2entrez_map)
        if pval < QVAL_CUTOFF:
            print "<{name}>: ({i}/{total}) {pval}\t{term}".format(name=p.name,
                i=i, pval=pval, term=annotations[term]['name'], total=total)
        results[term] = pval
    return results, diffexp


def multitest_correction(dataset, ontology, annotation_files):
    annotation_years = (json.load(open(f)) for f in annotation_files)
    factor = 'disease state'
    db = get_connection(100)
    for annotations in annotation_years:
        year = annotations['meta']['year']
        # If we didn't shuffle the annotations, the shuffle level is 0
        shuffled = annotations['meta'].get('shuffled', 0.0)
        for subset in dataset.factors[factor]:
            print"[%s]-[%s]-[%s]-[%s]-[%f]:" % (dataset.id, year, 
                ontology, subset, shuffled),
            with closing(db.cursor()) as c:
                _id = (dataset.id, subset, year, ontology, shuffled)
                print "selecting pvals... ",
                c.execute(select_pvals_sql, _id)
                # list of tuples [(_subid, pval), ...]
                results = list(c.fetchall())
            pvals = [x[1] for x in results] 
            subids = [_id + (x[0],) for x in results]
            print "calculating FDR... ",
            rejected, qvals = multitest.fdrcorrection(pvals)
            results = [(qvals[i],) + subids[i] for i, v in enumerate(qvals)]
            with closing(db.cursor()) as c:
                print "inserting %d qvals... " % len(results),
                c.executemany(insert_qval_sql, results)
                db.commit()
            print "done."
    db.close()


def main(file_or_accn, annotation_files, ontology):
    
    # this file can be downloaded from Uniprot's mapping service
    uniprot2entrez_map = json.load(open(MAPFILE))
    assert len(uniprot2entrez_map) > 27000

    # import the dataset
    dataset = fetch(file_or_accn, destdir='data')
    dataset = dataset.to_numeric()
    dataset.filter().log2xform()

    # import the annotation files (in JSON format)
    annotation_years = (json.load(open(f)) for f in annotation_files)

    # acquire the platform used from the dataset metadata
    platform = fetch(dataset.meta['platform'], destdir='data')

    print("Detected %d cores, splitting into %d subprocesses..." 
        % (NCORES, NCORES))

    if FDR_CORRECTION:
        jobs = []
        for annofile in annotation_files:
            p = Process(target=multitest_correction,
                       args=(dataset, ontology, [annofile]))
            jobs.append(p)
            p.start()
        [p.join() for p in jobs]
        return

    for annotations in annotation_years:
        year = annotations['meta']['year']
        annos = annotations['anno']
        shuffled = annotations['meta'].get('shuffled', 0.0)
        if FILTER_SIMILAR:
            print("Filtering out terms with less than a %d-gene "
                "difference from their parents" % MIN_VARIANCE)
            annos = filter_similar_terms(annos, MIN_VARIANCE)
        if FILTER_BY_DEPTH:
            print("Filtering out terms with fewer than %d "
                "or greater than %d parents" % (MIN_DEPTH, MAX_DEPTH))
            annos = filter_annos_by_depth(annos, MIN_DEPTH, MAX_DEPTH)
        if FILTER_BY_SIZE:
            print("Filtering out annotation gene sets greater than %d "
                "and less than %d" % (ANNO_MAX_SIZE, ANNO_MIN_SIZE))
            annos = filter_annos(annos, ANNO_MAX_SIZE, ANNO_MIN_SIZE)
        filtered_annotations = restrict_subontology(annos, ontology, year)
        blocks = split(filtered_annotations, blocks=NCORES)
        print("Split %d annotations into %d blocks of ~%d terms each..." 
            % (len(filtered_annotations), len(blocks), len(blocks[0])))
        # We're only looking at one factor for this analysis
        # Iterate over factors if this is no longer true
        factor = 'disease state'
        for subset in dataset.factors[factor]:
            print("-- [year: %s] [dataset: %s] [%s: %s] --" 
                % (year, dataset.id, factor, subset))
            jobs = []
            for block in blocks:
                p = Process(target=enriched, 
                    args=(dataset, platform, factor, subset, block, year, 
                        shuffled, len(filtered_annotations), ontology, 
                        uniprot2entrez_map))
                jobs.append(p)
                p.start()
            [p.join() for p in jobs]  # wait for them all to finish


def print_usage():
    import textwrap
    usage = """Conducts an enrichment analysis on the given dataset
using the specified GO sub-ontology and annotation files. The enriched
terms are then stored in a MySQL database along with p-values and other data.
Specific options are defined in config files, which can be passed to the
program with the --config <CFG_FILE> option. By default, the script looks for
configs/settings.cfg."""

    print(textwrap.fill(usage, 79))
    print("\nSee README.md for more detailed information.\n")

if __name__ == '__main__':
    parser = OptionParser(
        usage='%prog [options] <GEO dataset> <anno file 1 [anno file 2...]>')
    config = ConfigParser()

    if '--config' in sys.argv:
        cfg_file = sys.argv[sys.argv.index('--config') + 1]
        config.read(cfg_file)
    else:
        config.read("configs/settings.cfg")

    parser.add_option('--config', action='store', dest='config', 
        default='configs/settings.cfg', help="Alternate configuration file")
    parser.add_option('-o', action='store', type='choice', 
        choices=['MF', 'CC', 'BP'], dest='ontology', 
        help="REQUIRED: GO sub-ontology to use (MF, CC, BP)")
    parser.add_option('--fdr_correction', action='store_true',
        default=False, dest='fdrcorr', 
        help="Calculate p-values instead of doing EA")
    parser.add_option('--use_shuffled', action='store_true', 
        default=False, dest='shuffled', help="Work with shuffled annotations")
    parser.add_option('--filter_by_size', action='store_true', 
        default=False, dest='filter_size', 
        help=("Filter out term that have more or less than the specified max "
            "and min gene set sizes"))
    parser.add_option('--filter_by_depth', action='store_true', 
        default=False, dest='filter_depth', 
        help=("Filter terms by ontology depth (terms above/below a certain "
            "depth are omitted)"))
    parser.add_option('--filter_similar', action='store_true', 
        default=False, dest='filter_similar',
        help=("Filter out terms that have less than a certain number of genes "
            "not in common with their parents"))
    parser.add_option('--max_depth', action='store', type=int, 
        dest='max_depth', default=config.getint('Annotations', 'max depth'), 
        help="Maximum ontology depth of terms")
    parser.add_option('--min_depth', action='store', type=int, 
        dest='min_depth', default=config.getint('Annotations', 'min depth'), 
        help="Minimum ontology depth of terms")
    parser.add_option('--max_size', action='store', type=int, dest='max_size', 
        default=config.getint('Annotations', 'max size'), 
        help="Maximum number of genes annotated to term")
    parser.add_option('--min_size', action='store', type=int, dest='min_size', 
        default=config.getint('Annotations', 'min size'), 
        help="Minimum number of genes annotated to term")
    parser.add_option('--min_var', action='store', type=int, 
        dest='min_variance', default=config.getint('Annotations', 
            'min variance'), 
        help=("Minimum number of genes a child term must have different from a"
            " parent"))
    parser.add_option('--max_fdr', action='store', type=float, dest='max_fdr', 
        default=config.getfloat('FDR', 'cutoff'), 
        help="FDR q-value cutoff for defining differentially expressed genes")
    parser.add_option('--sql_table', action='store', dest='sql_table', 
        default=config.get('MySQL', 'table'), 
        help=("Table to store results (other MySQL options specified in "
            "config file)"))

    opts, args = parser.parse_args()

    if not opts.ontology:
        print_usage()
        parser.print_help()
        sys.exit(1)

    file_or_accn = args[0]
    annotation_files = args[1:]

    print annotation_files

    ontology = opts.ontology

    SHUFFLED = opts.shuffled
    FDR_CORRECTION = opts.fdrcorr

    FILTER_BY_SIZE = opts.filter_size
    FILTER_BY_DEPTH = opts.filter_depth
    FILTER_SIMILAR = opts.filter_similar

    ANNO_MAX_SIZE = opts.max_size
    ANNO_MIN_SIZE = opts.min_size

    MAX_DEPTH = opts.max_depth
    MIN_DEPTH = opts.min_depth

    MIN_VARIANCE = opts.min_variance
    
    QVAL_CUTOFF = opts.max_fdr
    NCORES = multiprocessing.cpu_count()

    MAPFILE = 'data/uniprot2entrez.json'

    # MySQL settings
    MYUSER = config.get('MySQL', 'user')
    MYHOST = config.get('MySQL', 'host')
    MYPASS = config.get('MySQL', 'pass')
    MYDB = config.get('MySQL', 'db')
    table = opts.sql_table

    # set SQL table to insert results into
    store_results_sql = store_results_sql.format(table=table)
    select_pvals_sql = select_pvals_sql.format(table=table)
    insert_qval_sql = insert_qval_sql.format(table=table)

    main(file_or_accn, annotation_files, ontology)
