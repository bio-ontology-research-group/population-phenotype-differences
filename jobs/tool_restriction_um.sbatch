#!/bin/bash
# Rank paired cases with Exomiser / LIRICAL (phenotype-only) under all / bedside / random profiles.
# sbatch --array=0-(N-1)%150 --mem=.. jobs/tool_restriction_um.sbatch TOOL IDS TAG XMX
#SBATCH --job-name=pavs-tr
#SBATCH --partition=debug
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=12:00:00
#SBATCH --exclude=node007,node008,unimatrix
#SBATCH --output=/data/hohndor/pavs-popmodel/tool_restriction/log/%x-%A_%a.out
set -u
TOOL=$1; IDS=$2; TAG=$3; XMX=${4:-4g}
W=/data/hohndor/pavs-popmodel; T=/data/hohndor/pavs-tools; R=$W/tool_restriction; cd $W
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=$T/pylib  # pyarrow 17 for Exomiser PARQUET
K=$SLURM_ARRAY_TASK_ID; N=$SLURM_ARRAY_TASK_COUNT
mkdir -p $R/$TAG $R/log/$TAG
echo "host $(hostname) TMPDIR=$TMPDIR tool=$TOOL tag=$TAG chunk=$K/$N"
venv/bin/python tool_restriction_run.py --tool $TOOL --profiles $R/profiles.jsonl --ids $IDS --g2p hpo/genes_to_phenotype.txt \
  --chunks $N --chunk $K --work $TMPDIR/tr_${TOOL}_$K --out $TMPDIR/${TOOL}_chunk$K.csv --logdir $R/log/$TAG \
  --xmx $XMX --cpus $SLURM_CPUS_PER_TASK --universe-n ${UNIVERSE_N:-0} \
  --exo-jar $T/exomiser/exomiser-cli-15.1.0/exomiser-cli-15.1.0.jar --exo-props $R/exomiser_application.properties \
  --exo-default $R/exomiser_hiphive_default.yml --exo-human $R/exomiser_hiphive_human.yml \
  --lirical-jar $T/lirical/lirical-cli-2.4.2/lirical-cli-2.4.2.jar --lirical-data $T/lirical/data
rc=$?
cp $TMPDIR/${TOOL}_chunk$K.csv $R/$TAG/ && echo "TASKDONE $K rc=$rc"
