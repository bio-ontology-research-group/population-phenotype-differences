#!/bin/bash
#SBATCH --job-name=hpo-select
#SBATCH --partition=debug
#SBATCH --nodelist=unimatrix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --output=/home/hohndor/pavs-origin/map2/slurm-select-%j.out
cd /home/hohndor/pavs-origin
/usr/bin/python3 -I qwen_tasks.py select --cands map2/all_cands.jsonl --out map2/all_select.jsonl --key-file /home/hohndor/.borg_llm_key --jobs 16 --batch 10
echo "DONE $(date)"
