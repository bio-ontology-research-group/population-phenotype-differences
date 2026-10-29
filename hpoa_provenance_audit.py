#!/usr/bin/env python3
"""Audit frequency formats and evidence-linked provenance; do not infer penetrance.

Case metadata are machine-derived/curated source data, not human validation.
Publication or disease overlap identifies review candidates, not verified cohort
membership. No country of an author is used as patient ancestry.
"""
import argparse
import collections
import csv
import hashlib
import json
import re
from pathlib import Path


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader((l for l in f if not l.startswith('#')),delimiter='\t'))


def write(path,rows):
    if not rows: return
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def frequency(value):
    if value in ('','-'): return 'missing',None,None
    m=re.fullmatch(r'(\d+)/(\d+)',value)
    if m:
        n,d=map(int,m.groups())
        return ('count' if 0<=n<=d and d>0 else 'invalid_count'),n,d
    if re.fullmatch(r'\d+(\.\d+)?%',value):
        return ('percentage' if 0<=float(value[:-1])<=100 else 'invalid_percentage'),None,None
    if value in {f'HP:004028{i}' for i in range(6)}: return 'category',None,None
    return 'unrecognized',None,None


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    base=a.root/'results/phenopackets'
    index=read_tsv(base/'index.tsv')
    fine={r['id']:r for f in (a.root/'results/v4/fine2').glob('*.tsv') for r in read_tsv(f)}
    cases=[]
    for r in index:
        path=base/r['file'];d=json.loads(path.read_text())
        ext={x['name']:x.get('value',{}) for x in d.get('subject',{}).get('extensions',[])}
        def field(name):
            v=ext.get(name,{})
            return v.get('text') or v.get('label') or ''
        disease_ids=set();alleles=set();expressions=set();variantclasses=set()
        for it in d.get('interpretations',[]):
            diag=it.get('diagnosis',{})
            if diag.get('disease',{}).get('id'): disease_ids.add(diag['disease']['id'])
            for gi in diag.get('genomicInterpretations',[]):
                vd=gi.get('variantInterpretation',{}).get('variationDescriptor',{})
                if vd.get('allelicState',{}).get('label'): alleles.add(vd['allelicState']['label'])
                for e in vd.get('expressions',[]):
                    if e.get('value'): expressions.add(e.get('syntax','')+':'+e['value'])
                if vd.get('structuralType',{}).get('label'): variantclasses.add(vd['structuralType']['label'])
        for disease in d.get('diseases',[]):
            if disease.get('term',{}).get('id'): disease_ids.add(disease['term']['id'])
        cases.append(dict(id=r['id'],pmid=r['pmid'],gene=r['gene'],corpus=r['corpus'],
                          region=fine.get(r['id'],{}).get('region',r['region']),
                          ancestry=field('ancestry'),origin_country=field('geographicOrigin'),
                          origin_level=field('originLevel'),region_source=field('populationRegionSource'),
                          care_country=field('placeOfCare') or r.get('place_of_care',''),
                          care_basis=field('placeOfCareBasis'),
                          disease_ids=';'.join(sorted(disease_ids)),
                          allelic_states=';'.join(sorted(alleles)),zygosity_reported=field('zygosityAsReported'),
                          variant_expressions=';'.join(sorted(expressions)),variant_classes=';'.join(sorted(variantclasses))))
    write(a.out/'case_metadata.csv',cases)
    counts=[]
    for corpus in sorted({r['corpus'] for r in cases}):
        cs=[r for r in cases if r['corpus']==corpus]
        counts.append(dict(corpus=corpus,n_cases=len(cs),n_publications=len({r['pmid'] for r in cs}),
                           ancestry=sum(bool(r['ancestry']) for r in cs),origin_country=sum(bool(r['origin_country']) for r in cs),
                           care_country=sum(bool(r['care_country']) for r in cs),disease_id=sum(bool(r['disease_ids']) for r in cs),
                           allelic_state=sum(bool(r['allelic_states']) for r in cs),
                           variant_expression=sum(bool(r['variant_expressions']) for r in cs),
                           variant_class=sum(bool(r['variant_classes']) for r in cs)))
    write(a.out/'genetic_metadata_coverage.csv',counts)
    bypub=collections.defaultdict(list)
    for r in cases: bypub[r['pmid']].append(r)
    hpo=a.root/'data/hpo/phenotype.hpoa'
    annotations=read_tsv(hpo)
    dg=collections.defaultdict(set)
    for r in read_tsv(a.root/'data/hpo/genes_to_phenotype.txt'):
        dg[r['disease_id']].add(r['gene_symbol'])
    rows=[];formats=collections.Counter();diseases=collections.defaultdict(list)
    for i,r in enumerate(annotations):
        f,n,den=frequency(r['frequency']);formats[f]+=1
        refs=re.findall(r'PMID:(\d+)',r['reference'])
        references=[v.strip() for v in r['reference'].split(';') if v.strip()]
        scope=('single_pmid_only' if len(references)==1 and len(refs)==1 else 'multiple_or_mixed_references_unresolved' if refs else 'no_pmid')
        linked=[c for ref in refs for c in bypub.get(ref,[])]
        matched=[c for c in linked if r['database_id'] in c['disease_ids'].split(';')]
        status='disease_and_publication_overlap' if matched else 'publication_overlap_only' if linked else 'unlinked_pmid' if refs else 'no_pmid_reference'
        support=matched or linked
        def unique(k): return ';'.join(sorted({c[k] for c in support if c[k]}))
        row=dict(annotation_row=i+1,disease_id=r['database_id'],disease_name=r['disease_name'],
                 genes=';'.join(sorted(dg[r['database_id']])),hpo_id=r['hpo_id'],qualifier=r['qualifier'],
                 aspect=r['aspect'],evidence=r['evidence'],onset=r['onset'],sex=r['sex'],
                 frequency=r['frequency'],frequency_format=f,affected=n,assessed=den,
                 reference=r['reference'],pmids=';'.join(refs),n_pmid_references=len(refs),
                 frequency_scope=scope,provenance_status=status,
                 linked_cases=len(support),linked_ancestries=unique('ancestry'),linked_origins=unique('origin_country'),
                 linked_care_countries=unique('care_country'),linked_care_bases=unique('care_basis'),
                 linked_case_regions=unique('region'),cohort_provenance_verified=False)
        rows.append(row);diseases[r['database_id']].append(row)
    write(a.out/'frequency_provenance.csv',rows)
    shortlist=[]
    for disease,rr in diseases.items():
        count=[r for r in rr if r['frequency_format']=='count' and r['aspect']=='P' and not r['qualifier']]
        pubs={p for r in count for p in r['pmids'].split(';') if p}
        linked=[r for r in count if r['provenance_status']=='disease_and_publication_overlap']
        if len(pubs)<2: continue
        shortlist.append(dict(disease_id=disease,disease_name=rr[0]['disease_name'],genes=rr[0]['genes'],
                              count_annotations=len(count),count_publications=len(pubs),
                              single_pmid_count_annotations=sum(r['frequency_scope']=='single_pmid_only' for r in count),
                              single_pmid_count_publications=len({r['pmids'] for r in count if r['frequency_scope']=='single_pmid_only'}),
                              disease_matched_annotations=len(linked),
                              disease_matched_publications=len({p for r in linked for p in r['pmids'].split(';') if p}),
                              max_term_denominator=max((r['assessed'] for r in count),default=0),
                              known_case_regions=';'.join(sorted({p for r in linked for p in r['linked_case_regions'].split(';') if p and p not in ('UNKNOWN','MIXED')}))))
    shortlist.sort(key=lambda r:(-r['disease_matched_publications'],-r['count_publications'],r['disease_id']))
    write(a.out/'disorder_candidates.csv',shortlist)
    summary={'hpoa_sha256':hashlib.sha256(hpo.read_bytes()).hexdigest(),'hpoa_version':next(l.strip() for l in hpo.open() if l.startswith('#version')),
             'annotations':len(rows),'frequency_formats':dict(formats),
             'provenance_status':dict(collections.Counter(r['provenance_status'] for r in rows)),
             'count_annotation_provenance':dict(collections.Counter(r['provenance_status'] for r in rows if r['frequency_format']=='count')),
             'count_frequency_scope':dict(collections.Counter(r['frequency_scope'] for r in rows if r['frequency_format']=='count')),
             'unique_count_pmids':len({p for r in rows if r['frequency_format']=='count' for p in r['pmids'].split(';') if p}),
             'disorders_with_at_least_two_count_publications':len(shortlist),
             'interpretation':'Overlap nominates source cohorts for review; it does not assign case ancestry to an HPOA denominator. No cohort provenance has yet been human-verified.'}
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
