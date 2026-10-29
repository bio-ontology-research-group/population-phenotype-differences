#!/usr/bin/env python3
"""Propose source-linked assessment denominators for subsequent source review.

Private model output is an extraction aid, never a validated annotation.
Uses the same logged request/finish checks as the paired-source experiment.
"""
import argparse, concurrent.futures, csv, hashlib, json
from pathlib import Path
import paired_source_extract as extractor

PROMPT = '''Read the supplied research article for the target gene. Identify PRIMARY human cohorts with a genetically confirmed disorder involving that gene. Distinguish primary patients from literature reviews, unaffected relatives, cell experiments and mixed-gene diagnostic cohorts. Keep genetically unresolved patients separate. Never assume every member of a cohort was assessed for every phenotype. Unmentioned findings are unknown, not absent. Do not infer ancestry from recruitment location.

For KMT2D, recognize the historical name MLL2 used in Kabuki syndrome reports; KDM6A is a different gene.

Return JSON with these fields:
patients: [] (always empty for this task)
cohorts: [{cohort_id, primary_or_review, diagnosis, genetic_confirmation, n_patients, n_families, population_as_stated, care_setting_as_stated, recruitment, age_range, previously_reported_patients, evidence}]
cohort_frequency_statements: [{cohort_id, phenotype_text, hpo_label, affected, assessed, denominator_basis, assessed_patient_ids, affected_patient_ids, source_location, evidence, ambiguity}]
review_flags: [text]

Every evidence field must be a verbatim quote of at most 30 words. Only supply integer affected and assessed values if both are explicitly recoverable from the prose or a table of present/absent/unknown results. A table listing positive clinical features per patient does NOT show that unlisted features were assessed. Explicit all-patient findings may have the cohort denominator. A normal test does not imply every individual phenotype measured by that test was excluded. Keep tissue-specific findings and thresholds separate. Flag disagreements between abstract, table and body; retain conflicting statements separately with an ambiguity flag. Do not combine primary and reviewed cases or pool studies. Flag whether individual ages, ancestry and care can be linked to assessed patients. If only positive feature lists without assessed denominators are available, give an empty frequency list and explain in review_flags. In mixed-gene studies include only denominators explicitly recoverable for the target gene subgroup. Do not use measurements from experimental cell lines as clinical phenotype frequencies.
Specific counting rules: a single-patient case report CAN contribute 1/1 for an explicitly present finding or 0/1 for an explicitly assessed absence, with denominator_basis=single_case. For a table with patient columns and separate phenotype rows, count + as affected and assessed, - as assessed absence, and blank/NA/unknown/not-done as unassessed; the denominator is the number of + and - entries in THAT row. Do not reject such a row because other phenotypes were unassessed. A clinical-list cell containing only positive features does not establish absence of other features. For example, if two patients have one + and one - for cleft palate, report 1/2; if one has + and the other is unknown, report 1/1. Keep this distinction explicit in denominator_basis. These are conditional frequencies among described and assessed patients, not population prevalence estimates.
Target gene: {gene}
TEXT:
{text}
'''

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--raw',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--key-file',type=Path,required=True)
    ap.add_argument('--endpoint',required=True)
    ap.add_argument('--model',default='qwen3.8-27b')
    ap.add_argument('--jobs',type=int,default=2)
    a=ap.parse_args();a.out.parent.mkdir(parents=True,exist_ok=True)
    extractor.PROMPT=PROMPT
    done={}
    if a.out.exists():
        for l in a.out.open():
            d=json.loads(l)
            if d.get('ok'):done[d['key']]=d
    ph=hashlib.sha256((extractor.SYSTEM+'\n'+PROMPT).encode()).hexdigest()
    items=[]
    for r in csv.DictReader(a.manifest.open()):
        if r['status']!='ready':continue
        path=a.raw/f"{r['pmid']}.fulltext.txt"
        item=dict(key=f"{r['pmid']}:{r['gene']}:frequency",pmid=r['pmid'],gene=r['gene'],source='fulltext',path=str(path))
        if item['key'] in done:
            old=done[item['key']]
            if old['source_sha256']!=hashlib.sha256(path.read_bytes()).hexdigest() or old['prompt_sha256']!=ph or old['model']!=a.model:
                raise ValueError('Incompatible saved extraction')
        else:items.append(item)
    key=a.key_file.expanduser().read_text().strip()
    with a.out.open('a') as out,concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        futures=[pool.submit(extractor.execute,it,a,key) for it in items]
        for f in concurrent.futures.as_completed(futures):
            r=f.result()
            if r.get('ok') and not all(isinstance(r['result'].get(k),list) for k in ['cohorts','cohort_frequency_statements','review_flags']):
                r.update(ok=False,error='invalid_frequency_schema')
            out.write(json.dumps(r)+'\n');out.flush()
            print(r['key'],r['ok'],r.get('error',''),flush=True)

if __name__=='__main__':main()
