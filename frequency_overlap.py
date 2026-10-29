#!/usr/bin/env python3
"""Audit availability of comparable study-specific count annotations in HPOA.

Matching a reference and phenotype is not evidence of independent cohorts.
No population assignment or inferential test is made from these annotations.
"""
import argparse,collections,csv,json,re
from pathlib import Path
from hpoa_provenance_audit import frequency
from source_robustness import write

FIELDS=['database_id','hpo_id','onset','sex','modifier','aspect','qualifier']

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--hpoa',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    groups=collections.defaultdict(list);coarse=collections.defaultdict(set)
    with a.hpoa.open() as f:
        for r in csv.DictReader((l for l in f if not l.startswith('#')),delimiter='\t'):
            typ,k,n=frequency(r['frequency'])
            if typ!='count' or not re.fullmatch(r'PMID:\d+',r['reference']):continue
            if r['aspect']!='P' or r['qualifier']:continue
            row={key:r[key] for key in FIELDS}
            row.update(reference=r['reference'],affected=k,assessed=n,disease_name=r['disease_name'])
            groups[tuple(r[key] for key in FIELDS)].append(row)
            coarse[r['database_id'],r['hpo_id']].add(r['reference'])
    rows=[];studies=[]
    for key,items in sorted(groups.items()):
        refs={r['reference'] for r in items}
        if len(refs)<2:continue
        row=dict(zip(FIELDS,key));row.update(n_references=len(refs),n_annotation_rows=len(items),
            min_annotation_frequency=min(r['affected']/r['assessed'] for r in items),
            max_annotation_frequency=max(r['affected']/r['assessed'] for r in items),
            references=';'.join(sorted(refs)),independent_cohorts_verified=False,population_provenance_verified=False)
        rows.append(row);studies.extend(items)
    write(a.out/'comparable_annotation_groups.csv',rows);write(a.out/'study_annotations.csv',studies)
    target=['OMIM:256500','OMIM:256810','OMIM:147920','OMIM:613826']
    summary=dict(disease_phenotype_pairs_with_multiple_single_pmid_counts=sum(len(v)>=2 for v in coarse.values()),
                 comparable_groups_after_onset_sex_modifier_matching=len(rows),
                 selected_disorders={d:sum(r['database_id']==d for r in rows) for d in target},
                 scope='Reference-specific HPOA count annotations. Source independence, denominator validity and population identity require source verification.')
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
