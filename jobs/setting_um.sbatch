#!/bin/bash
#SBATCH --job-name=pavs-setting
#SBATCH --partition=debug
#SBATCH --nodelist=node008
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --array=2-15%14
#SBATCH --output=/data/hohndor/pavs-popmodel/terms/log/setting_%a.out
set -u
W=/data/hohndor/pavs-popmodel; cd $W
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
L=$(sed -n "${SLURM_ARRAY_TASK_ID}p" setting_tasks.txt); NAME=${L%% *}; ARGS=${L#* }
F=results/setting2
venv/bin/python population_model.py terms --cases $F/cases_abs.tsv --cases $F/cases_pmcp_map2.tsv $ARGS --strata total --min-cases 5 \
  --obo hpo/hp.obo --termdb ontofunc_hp --groups results/groups_lastauthor.tsv --termclass results/termclass.jsonl \
  --perms 2000 --show 10 --out terms/$NAME
echo "TASKDONE $NAME"
