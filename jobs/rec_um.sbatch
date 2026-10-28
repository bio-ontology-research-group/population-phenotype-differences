#!/bin/bash
#SBATCH --job-name=pavs-recovery
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --cpus-per-task=8
#SBATCH --mem=110G
#SBATCH --time=1-00:00:00
#SBATCH --output=/data/hohndor/pavs-popmodel/recovery-%j.out
set -u
W=/data/hohndor/pavs-popmodel
T=${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}
cd $W
cp $W/simcache_v4.npy $T/
F=results/fine
C="--cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_pavs.tsv --cases $F/cases_abs.tsv"
R="$C --hpo-dir hpo --cache $T/simcache_v4.npy --authors results/all_authors.tsv --groups results/groups_lastauthor.tsv --lib-corpus store,pmcp,pavs_lit,pubmed_abs --perms 2000"
for r in MENA EASIA SASIA LATAM AFR; do
  echo "########## gene recovery $r vs EUR  $(date)"
  venv/bin/python recovery_main.py $R --r1 $r --r2 EUR --out results/fine/recovery_$r 2>&1 | grep -v "sim rows\|queries [0-9]"
done
echo ALLDONE
