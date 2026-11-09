#!/bin/bash
#SBATCH --job-name=pavs-zeroshare
#SBATCH --partition=debug
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=06:00:00
#SBATCH --array=1-5
#SBATCH --output=/data/hohndor/pavs-popmodel/terms_noclass/log/%a.out
set -u
W=/data/hohndor/pavs-popmodel; cd $W
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
L=$(sed -n "${SLURM_ARRAY_TASK_ID}p" zeroshare_tasks.txt); NAME=${L%% *}; ARGS=${L#* }
F=results/fine2
venv/bin/python population_model.py terms --cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_abs.tsv \
  --r2 EUR $ARGS --obo hpo/hp.obo --termdb ontofunc_hp --groups results/groups_lastauthor.tsv --termclass results/termclass.jsonl \
  --require-classified --perms 2000 --show 10 --out terms_noclass/$NAME
echo "TASKDONE $NAME"
