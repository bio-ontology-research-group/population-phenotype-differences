#!/bin/bash
#SBATCH --job-name=pavs-terms
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --array=1-50%14
#SBATCH --output=/data/hohndor/pavs-popmodel/terms/log/%a.out
set -u
W=/data/hohndor/pavs-popmodel; cd $W
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
L=$(sed -n "${SLURM_ARRAY_TASK_ID}p" terms_tasks.txt); NAME=${L%% *}; ARGS=${L#* }
F=results/fine2
venv/bin/python population_model.py terms --cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_abs.tsv \
  --r2 EUR $ARGS --obo hpo/hp.obo --termdb ontofunc_hp --groups results/groups_lastauthor.tsv --termclass results/termclass.jsonl \
  --perms 2000 --show 10 --out terms/$NAME
echo "TASKDONE $NAME"
