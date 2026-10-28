#!/bin/bash
#SBATCH --job-name=abs-extract
#SBATCH --partition=debug
#SBATCH --nodelist=unimatrix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --output=/home/hohndor/pavs-origin/abs/slurm-%j.out
cd /home/hohndor/pavs-origin
/usr/bin/python3 -I qwen_tasks.py abstract --cands abs/candidates.jsonl --out abs/abstracts_extract.jsonl --key-file /home/hohndor/.borg_llm_key --jobs 12
echo "DONE $(date)"
