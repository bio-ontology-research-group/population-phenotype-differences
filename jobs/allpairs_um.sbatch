#!/bin/bash
#SBATCH --job-name=pavs-allpairs
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --cpus-per-task=36
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --output=/data/hohndor/pavs-popmodel/allpairs-%j.out
set -u
W=/data/hohndor/pavs-popmodel
T=${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}
cd $W
F=results/fine
C="--cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_pavs.tsv --cases $F/cases_abs.tsv"
[ -f $W/simcache_v4.npy ] && cp $W/simcache_v4.npy $T/
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6
echo "start $(date)"
venv/bin/python population_model.py allpairs $C --hpo-dir hpo --cache $T/simcache_v4.npy --workers 36 \
  --groups results/groups_lastauthor.tsv --pops MENA,EUR,EASIA,SASIA,LATAM,AFR --perms 1000 --boot 200 \
  --pair-workers 6 --out results/fine/popmodel_fine
[ -f $W/simcache_v4.npy ] || cp $T/simcache_v4.npy $W/
echo "end $(date)"
echo ALLDONE
