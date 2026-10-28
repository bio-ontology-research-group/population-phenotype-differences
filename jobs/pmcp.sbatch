#!/bin/bash
#SBATCH --job-name=pmcp-extract
#SBATCH --partition=debug
#SBATCH --nodelist=unimatrix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/hohndor/pavs-origin/pmcp/slurm-%j.out
set -u
P=/home/hohndor/pavs-origin/pmcp
W=${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}/work
mkdir -p "$W"
cd /home/hohndor/pavs-origin
# stage 1: per-patient triage + phenotype extraction (resumable; output kept on /home)
/usr/bin/python3 -I pmcp_extract.py patient --cands $P/genetic_candidates.json --out $P/patients.jsonl --key-file /home/hohndor/.borg_llm_key --jobs 12
# stage 2: article-level patient origin for articles with >=1 germline Mendelian patient
/usr/bin/python3 -I -c "
import json
p=set()
for l in open('$P/patients.jsonl'):
    d=json.loads(l)
    if d.get('ok') and (d['result'] or {}).get('mendelian_germline'): p.add(d['pmid'])
open('$P/mendelian_pmids.txt','w').write('\n'.join(sorted(p))+'\n'); print('mendelian articles', len(p))
"
[ -f $P/fetch.jsonl ] && cp $P/fetch.jsonl "$W"/ && mkdir -p "$W/text" && { [ -f $P/text.tar ] && tar -xf $P/text.tar -C "$W"; }
/usr/bin/python3 -I store_origin_extract.py fetch --pmids $P/mendelian_pmids.txt --workdir "$W"
cp "$W"/fetch.jsonl $P/ && tar -cf $P/text.tar -C "$W" text
[ -f $P/origin.jsonl ] && cp $P/origin.jsonl "$W"/
/usr/bin/python3 -I store_origin_extract.py llm --workdir "$W" --key-file /home/hohndor/.borg_llm_key --jobs 12
cp "$W"/origin.jsonl $P/
echo "DONE $(date)"
