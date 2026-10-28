#!/bin/bash
#SBATCH --job-name=pavs-allpairs
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=320G
#SBATCH --output=/ibex/scratch/projects/c2014/pavs-popmodel/allpairs-%j.out
set -u
module load python/3.11.0
export OMP_NUM_THREADS=4
cd /ibex/scratch/projects/c2014/pavs-popmodel
F=results/fine
C="--cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_pavs.tsv --cases $F/cases_abs.tsv"
python3 population_model.py allpairs $C --hpo-dir hpo --cache simcache_v4.npy --groups results/groups_lastauthor.tsv \
  --pops MENA,EUR,EASIA,SASIA,LATAM,AFR --perms 1000 --boot 200 --out results/fine/popmodel_fine 2>&1 | grep -v "sim rows"
echo ALLDONE
