#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Jobs Spawner: Spawns a new PBS job for each dataset in a given list of GDS accession numbers.
The job conducts an enrichment analysis of the dataset using historical Gene Ontology annotations,
storing the results in a SQLite database (generally 'results.db').

The number of concurrent jobs spawned can be defined on the command line.

Usage: python jobs_spawner.py <file listing GDS accns, one per line> [max concurrent jobs]

Author: eclarke@scripps.edu
"""

import sys
import subprocess

# Job-specific options
NODES = 1            # we don't have MPI capabilities so no point in increasing this
PROC_PER_NODE = 8    # max cores for garibaldi machines
WALLTIME = 3         # in hours
JOBNAME = "{GDS}-ea.job"
SCRIPT_FILE = "jobs/{jobname}.sh"
COMMAND = "qsub {seriesopts} {scriptfile}"

SCRIPT_CONTENTS = """
#PBS -l nodes={nodes}:ppn={proc_per_node}
#PBS -l mem=20gb
#PBS -l walltime={hours}:00:00
#PBS -N {name}
#PBS -j oe

cd eanalysis
python enrichment-hpc.py {gds} data/goa-*.json
"""

def spawn(gds, after=None, dryrun=False):
    args = {'nodes': NODES,
            'proc_per_node': PROC_PER_NODE,
            'hours': WALLTIME,
            'name': JOBNAME.format(GDS=gds),
            'script': SCRIPT_FILE.format(jobname=gds),
            'gds':gds}

    script_file = create_job_script(args)
    series_options = '-W depend=afterany:'+after if after else ''
    command = COMMAND.format(seriesopts=series_options,
                             scriptfile=script_file)
    jobid = subprocess.check_output(command.split(' ')) if not dryrun else script_file
    return jobid, command


def create_job_script(args):
    script = SCRIPT_CONTENTS.format(**args)
    with open(args['script'], 'wb') as scriptfile:
        scriptfile.write(script)
    return args['script']


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python jobs_spawner.py <file listing GDS accns, one per line> [max concurrent jobs] [dryrun]")
        sys.exit(1)

    with open(sys.argv[1]) as infile:
        accns = [x.strip('\n') for x in infile.readlines()]

    max_concurrent = int(sys.argv[2]) if len(sys.argv) > 2 else 10

    dryrun = 'dryrun' in sys.argv

    jobs = []
    after = None
    for i, accn in enumerate(accns):
        after = jobs[i-1] if i > 0 and (i % max_concurrent) == 0 else after
        job, command = spawn(accn, after=after, dryrun=dryrun)
        print "launched with command: ", command
        jobs.append(job)
#    for job in jobs:
 #       print job



    
        
