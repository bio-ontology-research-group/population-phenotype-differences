#!/bin/bash
#SBATCH --job-name=pavs-popv4
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=320G
#SBATCH --output=/ibex/scratch/projects/c2014/pavs-popmodel/popv4-%j.out
set -u
module load python/3.11.0
export OMP_NUM_THREADS=4
cd /ibex/scratch/projects/c2014/pavs-popmodel
C="--cases results/cases_store.tsv --cases results/cases_pmcp_map2.tsv --cases results/cases_pavs.tsv --cases results/cases_abs.tsv"
CACHE=simcache_v4.npy
for r in MENA EASIA; do
  echo "########## [v4: + abstracts] $r vs EUR: decomposition"
  python3 population_model.py pairs $C --hpo-dir hpo --cache $CACHE --groups results/groups_lastauthor.tsv --r1 $r --r2 EUR --perms 1000 --boot 200 --out results/popmodel_v4_$r 2>&1 | grep -v "sim rows\|^cases"
done
LIT=store,pmcp,pavs_lit,pubmed_abs
RC="$C --hpo-dir hpo --cache $CACHE --authors results/all_authors.tsv --groups results/groups_lastauthor.tsv --perms 1000"
for spec in "MENA:lit:$LIT" "MENA:clinSA:pavs_clin" "EASIA:lit:$LIT" "MENA:abs:pubmed_abs" "EASIA:abs:pubmed_abs"; do
  IFS=: read r n q <<< "$spec"
  echo "########## [v4] retrieval $r vs EUR, queries=$q"
  python3 retrieval_perm.py $RC --r1 $r --r2 EUR --query-corpus $q --lib-corpus $LIT --out results/popretr_v4_${r}_$n 2>&1 | grep -v "sim rows\|^cases"
done
echo ALLDONE
