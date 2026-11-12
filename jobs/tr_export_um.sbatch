#!/bin/bash
#SBATCH --job-name=pavs-trexport
#SBATCH --partition=debug
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --exclude=node008
#SBATCH --output=/data/hohndor/pavs-popmodel/tool_restriction/log/export-%j.out
set -u
W=/data/hohndor/pavs-popmodel; cd $W
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
F=results/fine2
echo "TMPDIR=${TMPDIR:-unset}"
venv/bin/python export_profiles.py --cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_abs.tsv \
  --hpo-dir hpo --termclass results/termclass.jsonl --rand 'gda_rand/gdarand_chunk*.csv' --out tool_restriction/profiles.jsonl
echo "EXIT $?"
