#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Jobs Spawner: Spawns a new PBS job for each dataset in a given list of GDS accession numbers.
The specifics of the job are specified in a config file (default: config/job.settings.cfg) which can
be specified using the --config option.

To be sure that the config file contains the correct fields, run with --dryrun first. This runs 
through the entire script except for submitting the job to the cluster.

Usage: python jobs_spawner <list of datasets> [--config CONFIG_FILE] [--dryrun]

Author: eclarke@scripps.edu
"""

import sys
import subprocess
from ConfigParser import SafeConfigParser
from argparse import ArgumentParser


def spawn(gds, config, dryrun):
    pyscript = config.get('Job', 'pyscript')
    jobname = config.get('Job', 'name').format(gds=gds)
    template = open(config.get('Job', 'template')).read()
    args = dict(config.items('Template')+config.items('Job'))
    args['gds'] = gds
    outfile = config.get('Job', 'jobscript')
    
    script_file = create_job_script(template, args, outfile)

    command = config.get('Job', 'command').format(gds=gds)
    jobid = subprocess.check_output([x for x in command.split(' ') if x]) if not dryrun else script_file
    return jobid, command


def create_job_script(script, args, outfile):
    script = script.format(**args)
    with open(outfile, 'wb') as out:
        out.write(script)
    return outfile


def main():
    parser = ArgumentParser(description="Spawn jobs on Garibaldi using the Torque queue system. Job options specified in config file.")
    parser.add_argument('--config', action='store', dest="cfg_file", help="File with configuration options", default='configs/job.settings.cfg')
    parser.add_argument('--dryrun', action='store_true', default=False, dest="dryrun", help="Create job script in jobs/ but do not launch on cluster.")
    parser.add_argument('-d', nargs='+', action='store', dest="datasets", help="One or more GEO dataset accessions")
    parser.add_argument('-f', action='store', dest="datasets_file", help="File listing GEO datasets, one per line", type=file)

    args = parser.parse_args()

    config = SafeConfigParser()
    config.read(args.cfg_file)

    if args.datasets:
        accessions = args.datasets
    elif args.datasets_file:
        accessions = [x.strip("\n") for x in args.datasets_file]
        args.datasets_file.close()
    else:
        print "Must specify GEO datasets, either in file or on command line."
        parser.print_usage()
        return

    for i, accn in enumerate(accessions):
        job, command = spawn(accn, config, args.dryrun)
        print "%s launched with command: '%s'" % (job, command)
    

if __name__ == "__main__":
    main()
