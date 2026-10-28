#!/bin/bash
#SBATCH --job-name=storeabs
#SBATCH --partition=debug
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=storeabs-%j.out
# Qwen abstract extraction on the abstracts of Phenopacket Store publications (extraction benchmark).
cd /data/hohndor/pavs-popmodel/storeabs
STAGE=${STAGE:-abstract}
if [ "$STAGE" = abstract ]; then
  ../venv/bin/python qwen_tasks.py abstract --cands cands.jsonl --out abs_store.jsonl --key-file ~/.borg_llm_key --jobs 12
else
  ../venv/bin/python qwen_tasks.py select --cands map_cands.jsonl --out sel_store.jsonl --key-file ~/.borg_llm_key --jobs 12
fi
echo DONE
