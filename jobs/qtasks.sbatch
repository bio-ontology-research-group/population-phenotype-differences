#!/bin/bash
#SBATCH --job-name=qwen-tasks
#SBATCH --partition=debug
#SBATCH --nodelist=unimatrix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --output=/home/hohndor/pavs-origin/pmcp/slurm-qtasks-%j.out
cd /home/hohndor/pavs-origin
/usr/bin/python3 -I qwen_tasks.py termclass --terms pmcp/termclass_terms.tsv --out pmcp/termclass.jsonl --key-file /home/hohndor/.borg_llm_key --jobs 8
/usr/bin/python3 -I qwen_tasks.py setting --cands pmcp/genetic_candidates.json --patients pmcp/patients.jsonl --affil pmcp/pmcp_pmid_country.tsv --out pmcp/setting.jsonl --key-file /home/hohndor/.borg_llm_key --jobs 12
echo "DONE $(date)"
