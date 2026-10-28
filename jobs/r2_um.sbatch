#!/bin/bash
#SBATCH --job-name=pavs-r2
#SBATCH --partition=debug
#SBATCH --nodelist=node008
#SBATCH --cpus-per-task=36
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH --output=/data/hohndor/pavs-popmodel/r2-%j.out
# Variance decomposition of within-gene pair similarity for all 15 population pairs.
set -u
W=/data/hohndor/pavs-popmodel
T=${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}
cd $W
F=results/fine
C="--cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_pavs.tsv --cases $F/cases_abs.tsv"
cp $W/simcache_v4.npy $T/
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6
echo "start $(date)"
venv/bin/python population_model.py allpairs $C --hpo-dir hpo --cache $T/simcache_v4.npy --workers 8 \
  --groups results/groups_lastauthor.tsv --pops MENA,EUR,EASIA,SASIA,LATAM,AFR --r2-only --r2-boot 50 \
  --pair-workers 6 --out results/fine/popmodel_r2
echo "end $(date)"
echo ALLDONE
