#!/bin/bash
#SBATCH --job-name=pavs-recover
#SBATCH --partition=debug
#SBATCH --nodelist=unimatrix
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=/home/hohndor/pavs-origin/recover/slurm-%j.out
set -u
R=/home/hohndor/pavs-origin/recover
W=${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}/work
mkdir -p "$W"
cd /home/hohndor/pavs-origin
/usr/bin/python3 -I store_origin_extract.py fetch-oa --openalex "$R/unknown_recover_openalex.json" --workdir "$W"
cp "$W"/fetch.jsonl "$R"/ && tar -cf "$R"/text.tar -C "$W" text
/usr/bin/python3 -I store_origin_extract.py llm --workdir "$W" --key-file /home/hohndor/.borg_llm_key --jobs 12
cp "$W"/origin.jsonl "$R"/
echo "DONE $(date)"
