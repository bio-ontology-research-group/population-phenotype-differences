#!/bin/bash
#SBATCH --job-name=pavs-gdarand
#SBATCH --partition=debug
#SBATCH --array=0-39
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=/data/hohndor/pavs-popmodel/gda_rand/log/gdarand-%A_%a.out
set -u
W=/data/hohndor/pavs-popmodel; cd $W
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
F=results/fine2
venv/bin/python gda_impact.py rank --random-control --cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_abs.tsv \
  --hpo-dir hpo --termclass results/termclass.jsonl --chunks 40 --chunk $SLURM_ARRAY_TASK_ID --out $TMPDIR/gdarand
cp $TMPDIR/gdarand_chunk${SLURM_ARRAY_TASK_ID}.csv gda_rand/
echo "TASKDONE $SLURM_ARRAY_TASK_ID"
