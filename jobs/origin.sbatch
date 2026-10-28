#!/bin/bash
#SBATCH --job-name=pavs-origin
#SBATCH --partition=debug
#SBATCH --nodelist=unimatrix
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=/home/hohndor/pavs-origin/slurm-%j.out
set -u
W=${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}/work
mkdir -p "$W" /home/hohndor/pavs-origin
cd /home/hohndor/pavs-origin
# resume support: restore previous partial results if present
[ -f fetch.jsonl ] && cp fetch.jsonl "$W"/ && mkdir -p "$W/text" && { [ -f text.tar ] && tar -xf text.tar -C "$W"; }
[ -f origin.jsonl ] && cp origin.jsonl "$W"/
/usr/bin/python3 -I store_origin_extract.py fetch --pmids store_pmids_union.txt --workdir "$W"
cp "$W"/fetch.jsonl . && tar -cf text.tar -C "$W" text
/usr/bin/python3 -I store_origin_extract.py llm --workdir "$W" --key-file /home/hohndor/.borg_llm_key --jobs 12
cp "$W"/origin.jsonl .
echo "DONE $(date)"
