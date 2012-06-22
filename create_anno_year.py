#!/usr/bin/env python
# -*- coding: utf-8 -*-
import Annotations
from Annotations import import_annotations
import json
import os
import glob

Annotations.verbose = True

# Downloaded from Uniprot, lists each obsolete Uniprot ID and its fate. (deleted, merged, demerged, etc)
obsoletes = 'data/obsolete.list'
# The name of the flat ontology structure from a given year, produced by go_flattener.jar (eats .owl files)
flatfile = 'data/go-%d.flat'
# A current GO structure (in obo format) to reference GOID->names
go_ext_obo = 'data/gene_ontology_ext.obo'

assert all([os.path.isfile(x) for x in (obsoletes, flatfile, go_ext_obo)])

def main(years, namefmt, keep_iea):

    for year in years:
        goafile = glob.glob('goa_human-*-%d' % year)
        if goafile:
            goafile = goafile[0]
            print("Found goafile for %d at %s" % (year, goafile))
        else:
            print("Could not find human goa file for %d, skipping" % year)
            continue
        annotations = import_annotations(goafile, obsfile, flatfile % year, go_ext_obo, str(year), keep_iea)
        with open(namefmt % year, 'wb') as out:
            json.dump(annotations, out)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 4:
        print """Usage: python create_anno_year.py <beginning year> <end year> <outputname (sub %d for year)> --keep_iea"""
        sys.exit(1)
    years = xrange(int(sys.argv[1]), int(sys.argv[2])+1)
    keep_iea = '--keep_iea' in sys.argv
    main(years, sys.argv[3], keep_iea)
