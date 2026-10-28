#!/bin/bash
#SBATCH --job-name=tdm-origin
#SBATCH --partition=debug
#SBATCH --nodelist=unimatrix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/home/hohndor/pavs-origin/tdm_store/slurm-%j.out
cd /home/hohndor/pavs-origin
/usr/bin/python3 -I store_origin_extract.py llm --workdir tdm_store --key-file /home/hohndor/.borg_llm_key --jobs 6
echo "DONE $(date)"
