#!/usr/bin/env python
import sys
import enrichment as ea
from multiprocessing import Process
from __init__ import fetch
def _usage():
    print("Corrects results from enrichment analysis for multiple comparison errors. Inserts q-values into results table.")
    print("Usage: python fdr_correction.py <GDS file or accn> <[MF, CC, BP]> <anno file 1> [more anno files...]")
    print("\n See README.md for more information.")
    sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) < 3 or len(sys.argv[2]) != 2:
        _usage()
    jobs = []
    dataset = fetch(sys.argv[1], destdir='data')

    for year in sys.argv[3:]:
        p = Process(target=ea.multitest_correction,
                    args=(dataset, sys.argv[2], [year]))
        jobs.append(p)
        p.start()
    [p.join() for p in jobs] # wait for them all to finish
