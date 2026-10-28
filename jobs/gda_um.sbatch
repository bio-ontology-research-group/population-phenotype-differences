#!/bin/bash
#SBATCH --job-name=pavs-gda
#SBATCH --partition=debug
#SBATCH --array=0-39
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=gda/log/gda-%A_%a.out
set -u
W=/data/hohndor/pavs-popmodel; cd $W
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
F=results/fine2
venv/bin/python gda_impact.py rank --cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_abs.tsv \
  --hpo-dir hpo --termclass results/termclass.jsonl --chunks 40 --chunk $SLURM_ARRAY_TASK_ID --out $TMPDIR/gda
cp $TMPDIR/gda_chunk${SLURM_ARRAY_TASK_ID}.csv gda/
echo "TASKDONE $SLURM_ARRAY_TASK_ID"
