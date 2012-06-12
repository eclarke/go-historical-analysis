#PBS -l nodes=1:ppn=8
#PBS -l mem=20gb
#PBS -l walltime=3:00:00
#PBS -N gds1962-ea.job
#PBS -j oe

cd eanalysis
python enrichment-hpc.py GDS1962 data/goa-*.json

