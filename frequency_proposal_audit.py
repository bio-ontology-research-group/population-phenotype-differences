#!/usr/bin/env python3
"""Audit model frequency proposals without treating schema checks as validation."""
import argparse,collections,hashlib,json,re
from pathlib import Path
from source_robustness import write

def norm(s):return re.sub(r'\s+',' ',str(s or '')).strip().casefold()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--extract',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    docs={}
    for l in a.extract.open():
        d=json.loads(l)
        if d.get('ok'):docs[d['key']]=d
    papers=[];proposals=[]
    for d in docs.values():
        source=norm((a.raw/f"{d['pmid']}.fulltext.txt").read_text())
        result=d['result'];cohorts={c.get('cohort_id'):c for c in result['cohorts']}
        papers.append(dict(pmid=d['pmid'],gene=d['gene'],proposed_cohorts=len(cohorts),
                           proposed_frequencies=len(result['cohort_frequency_statements']),
                           review_flags=len(result['review_flags']),source_verified=False))
        for i,f in enumerate(result['cohort_frequency_statements']):
            k,n=f.get('affected'),f.get('assessed');cohort=cohorts.get(f.get('cohort_id'),{})
            valid=type(k)==int and type(n)==int and n>0 and 0<=k<=n
            quote=norm(f.get('evidence'));quote_found=bool(quote) and quote in source
            proposals.append(dict(pmid=d['pmid'],gene=d['gene'],proposal_index=i,
                cohort_id=f.get('cohort_id'),phenotype=f.get('phenotype_text'),proposed_hpo_label=f.get('hpo_label'),
                affected=k,assessed=n,valid_integer_counts=valid,cohort_resolved=bool(cohort),
                primary_or_review=cohort.get('primary_or_review'),normalized_quote_found=quote_found,
                ambiguity_flag=bool(f.get('ambiguity')),status='requires_source_review'))
    write(a.out/'papers.csv',papers);write(a.out/'proposals.csv',proposals)
    summary=dict(successful_sources=len(docs),proposed_statements=len(proposals),
        count_schema_pass=sum(p['valid_integer_counts'] for p in proposals),
        resolved_cohort=sum(p['cohort_resolved'] for p in proposals),
        normalized_quote_found=sum(p['normalized_quote_found'] for p in proposals),
        primary_or_review=dict(collections.Counter(p['primary_or_review'] for p in proposals)),
        model=sorted({d['model'] for d in docs.values()}),prompt_sha256=sorted({d['prompt_sha256'] for d in docs.values()}),
        extraction_sha256=hashlib.sha256(a.extract.read_bytes()).hexdigest(),
        interpretation='Counts describe model proposals. Schema/quote checks do not validate genotype, disorder, cohort identity, assessment denominator, or independence. No model-only proposal is an accepted frequency annotation.')
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
