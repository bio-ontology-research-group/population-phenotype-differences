#!/usr/bin/env python3
"""Freeze a disease-focused frequency-curation set and retrieve supporting reports.

Diseases chosen for reference support, clinical diversity, and availability of
separate evaluation data, before comparing frequency-aware ranking performance.
"""
import argparse,collections,concurrent.futures,csv,hashlib,json,shutil
from pathlib import Path
from paired_source_inputs import fetch
from source_robustness import write

DISEASES={'OMIM:256500':('SPINK5','Netherton syndrome'),
          'OMIM:256810':('MPV17','Mitochondrial DNA depletion syndrome 6'),
          'OMIM:147920':('KMT2D','Kabuki syndrome 1'),
          'OMIM:613826':('RPGRIP1','Leber congenital amaurosis 6')}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.raw.mkdir(parents=True,exist_ok=True);a.out.mkdir(parents=True,exist_ok=True)
    annotations=list(csv.DictReader((a.root/'results/v5/hpoa_audit/frequency_provenance.csv').open()))
    meta=list(csv.DictReader((a.root/'results/v5/hpoa_audit/case_metadata.csv').open()))
    clinical=[r for r in meta if r['gene'] in {v[0] for v in DISEASES.values()} and r['corpus'] in ('pavs_clin','pavs_us')]
    heldout_pmids={r['pmid'] for r in clinical}
    groups=dict(csv.reader((a.root/'results/groups_lastauthor.tsv').open(),delimiter='\t'))
    heldout_groups={groups.get(p,'P:'+p) for p in heldout_pmids}
    candidates=[];excluded=[]
    for disease,(gene,name) in DISEASES.items():
        rr=[r for r in annotations if r['disease_id']==disease and r['frequency_format']=='count' and r['aspect']=='P']
        refs=collections.Counter(p for r in rr for p in r['pmids'].split(';') if p)
        for pmid,n in sorted(refs.items()):
            row=dict(pmid=pmid,gene=gene,disease_id=disease,disease_name=name,n_hpoa_count_rows=n,
                     reporting_group=groups.get(pmid,'P:'+pmid))
            if pmid in heldout_pmids or row['reporting_group'] in heldout_groups:
                excluded.append(dict(row,exclusion='clinical_evaluation_publication_or_reporting_group'))
            else:candidates.append(row)
    write(a.out/'curation_selection.csv',candidates);write(a.out/'excluded_references.csv',excluded)
    evalrows=[{k:r[k] for k in ['id','pmid','gene','corpus','region','disease_ids']} for r in clinical]
    write(a.out/'frozen_clinical_evaluation.csv',evalrows)
    plan={'selected_disorders':DISEASES,'selection_rationale':'Multiple source publications with counts; varied clinical phenotypes; clinical test cases available for three genes. SPINK5 is a frequency-provenance case study with no eligible clinical queries in the current collection.',
          'candidate_publication_disorder_pairs':len(candidates),'unique_candidate_publications':len({r['pmid'] for r in candidates}),
          'clinical_cases':dict(collections.Counter(r['gene']+':'+r['corpus'] for r in clinical)),
          'evaluation_rule':'All existing clinical PAVS queries for selected genes; their publications and last-author groups excluded from new curation. Eligibility, disease subtype and residual patient overlap require checks before evaluation.',
          'hpoa_evaluation_leakage_rule':'For baseline and updated frequency references, remove source contributions from evaluation publications/reporting groups; multi-reference aggregates involving these sources cannot be disaggregated and must be excluded.',
          'status':'frozen_before_frequency_model_evaluation'}
    (a.out/'design.json').write_text(json.dumps(plan,indent=2)+'\n')
    rows=[];cache=a.root/'data/v5_private/paired_sources'
    for r in candidates:
        for suffix in ['metadata.json','xml']:
            src=cache/f"{r['pmid']}.{suffix}";dst=a.raw/src.name
            if src.exists() and not dst.exists():shutil.copy2(src,dst)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(fetch,r,a.raw) for r in candidates]
        for f in concurrent.futures.as_completed(futures):
            r=f.result();rows.append(r);print(r['disease_id'],r['pmid'],r['status'],flush=True)
    fields=sorted({k for r in rows for k in r})
    with (a.out/'manifest.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(sorted(rows,key=lambda r:(r['disease_id'],r['pmid'])))
    (a.out/'availability.json').write_text(json.dumps(dict(collections.Counter(r['status'] for r in rows)),indent=2)+'\n')

if __name__=='__main__':main()
