#!/bin/bash
#SBATCH --job-name=pavs-recv4
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=320G
#SBATCH --output=/ibex/scratch/projects/c2014/pavs-popmodel/recv4-%j.out
set -u
module load python/3.11.0
cd /ibex/scratch/projects/c2014/pavs-popmodel
N=$(( $(cat results/cases_store.tsv results/cases_pmcp_map2.tsv results/cases_pavs.tsv results/cases_abs.tsv | wc -l) - 4 ))
# wait until the first job has written the full similarity cache (header + N*N float32)
until [ -f simcache_v4.npy ] && [ $(stat -c %s simcache_v4.npy) -ge $((N*N*4)) ] && ! lsof simcache_v4.npy >/dev/null 2>&1; do sleep 120; done
sleep 120
C="--cases results/cases_store.tsv --cases results/cases_pmcp_map2.tsv --cases results/cases_pavs.tsv --cases results/cases_abs.tsv"
R="$C --hpo-dir hpo --cache simcache_v4.npy --authors results/all_authors.tsv --groups results/groups_lastauthor.tsv --lib-corpus store,pmcp,pavs_lit,pubmed_abs --perms 2000"
for r in MENA EASIA; do
  echo "########## [v4] gene recovery $r vs EUR"
  python3 recovery_main.py $R --r1 $r --r2 EUR --out results/recovery_v4_$r 2>&1 | grep -v "sim rows\|queries [0-9]"
done
echo ALLDONE
