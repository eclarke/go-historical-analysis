#!/usr/bin/env python
import sys

NODES = 1            # we don't have MPI capabilities so no point in increasing this
PROC_PER_NODE = 8    # max cores for garibaldi machines
WALLTIME = 15        # in hours

JOBNAME = "{GDS}-{pyscript}-ea.job"
SCRIPT_FILE = "jobs/{jobname}.sh"
COMMAND = "qsub {seriesopts} {scriptfile}"

SCRIPT_CONTENTS = """
#PBS -l nodes={nodes}:ppn={proc_per_node}
#PBS -l mem=20gb
#PBS -l walltime={hours}:00:00
#PBS -N {name}
#PBS -j oe

cd go
python {pyscriptname} {gds} {ontology} {years}
"""

def spawn(gds, scriptname, ontology, years, after, dryrun):
    jobname = JOBNAME.format(GDS=gds, pyscript=scriptname)
    years = ' '.join(['data/goa-%s.json' % x for x in years])
    args = {'nodes': NODES,
            'proc_per_node': PROC_PER_NODE,
            'hours': WALLTIME,
            'name': jobname,
            'script': SCRIPT_FILE.format(jobname=jobname),
            'gds':gds,
            'pyscriptname': scriptname,
            'ontology': ontology,
            'years': years}

    script_file = create_job_script(args)
    series_options = '-W depend=afterany:'+after if after else ''
    command = COMMAND.format(seriesopts=series_options,
                             scriptfile=script_file)
    jobid = check_output([x for x in command.split(' ') if x]) if not dryrun else script_file
    return jobid, command

def create_job_script(args):
    script = SCRIPT_CONTENTS.format(**args)
    with open(args['script'], 'wb') as scriptfile:
        scriptfile.write(script)
    return args['script']


def _usage():
    print("Usage: python spawn_custom.py [-df] <file listing GDS accns, one per line> <ontology> <years>")
    print("""
    -d: dryrun
    -f: FDR correction (default enrichment analysis)
    (combine flags after one hyphen)""")

    sys.exit(1)

if __name__ == "__main__":
    try:
        from subprocess import check_output
    except ImportError:
        print "Could not import check_output from subprocess. Did you load Python 2.7?"
        _usage()

    if len(sys.argv) < 4:
        _usage()

    if sys.argv[1] in ["-df", "-d", "-f", "-fd"]:
        dryrun = "d" in sys.argv[1]
        fdr = "f" in sys.argv[1]
        infile = sys.argv[2]
    else:
        dryrun = fdr = False
        infile = sys.argv[1]
    

    with open(infile) as inf:
        accns = [x.strip('\n') for x in inf.readlines()]
        
    ontology = sys.argv[sys.argv.index(infile)+1]
    years = sys.argv[sys.argv.index(ontology)+1:]
    scriptname = 'enrichment.py' if not fdr else 'fdr_correction.py'
    for accn in accns:
        job, cmd = spawn(accn, scriptname, ontology, years, after=None, dryrun=dryrun)
        print "%s launched with command: '%s'" % (job, cmd)
        
        
        
