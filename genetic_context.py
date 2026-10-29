#!/usr/bin/env python3
"""Recover stated disease inheritance and descriptive HGVS protein consequences.

Consequences are syntax-based annotations, not predictions of pathogenicity or
molecular disease mechanism. Missing information stays unknown.
"""
import argparse,collections,csv,json,re
from pathlib import Path
from source_robustness import write
from hpo_map2 import load_obo
from phenotype_specificity import ancestors_fn


def protein_class(expression):
    if not expression.startswith('hgvs.p:') or ':p.' not in expression:return None
    p=expression.split(':p.',1)[1].strip('()')
    if not p or '?' in p:return None
    if 'fs' in p:return 'frameshift'
    if 'delins' in p:return 'protein_delins'
    if 'del' in p:return 'protein_deletion'
    if 'dup' in p:return 'protein_duplication'
    if 'ins' in p:return 'protein_insertion'
    if 'ext' in p:return 'protein_extension'
    if '=' in p:return 'synonymous'
    if re.fullmatch(r'[A-Z][a-z]{2}\d+(Ter|\*)',p):return 'stop_gained'
    if re.fullmatch(r'(Ter|\*)\d+[A-Z][a-z]{2}',p):return 'stop_lost'
    if re.fullmatch(r'[A-Z][a-z]{2}\d+[A-Z][a-z]{2}',p):return 'amino_acid_substitution'
    return None


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--metadata',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    T=load_obo(str(a.root/'data/hpo/hp.obo'));anc=ancestors_fn(T)
    inheritance=collections.defaultdict(set)
    with (a.root/'data/hpo/phenotype.hpoa').open() as f:
        for r in csv.DictReader((l for l in f if not l.startswith('#')),delimiter='\t'):
            if not r['qualifier'] and 'HP:0000005' in anc(r['hpo_id']) and r['hpo_id']!='HP:0000005':
                inheritance[r['database_id']].add(r['hpo_id'])
    rows=list(csv.DictReader(a.metadata.open()))
    for r in rows:
        classes={protein_class(v) for v in r['variant_expressions'].split(';')}-{None}
        # Require every available protein expression to be classified for matching.
        proteins=[v for v in r['variant_expressions'].split(';') if v.startswith('hgvs.p:')]
        r['protein_consequences']=';'.join(sorted(classes))
        r['protein_class_complete']='yes' if proteins and all(protein_class(v) for v in proteins) else 'no'
        inherit=set().union(*(inheritance[d] for d in r['disease_ids'].split(';')))
        r['disease_inheritance_ids']=';'.join(sorted(inherit))
        r['disease_inheritance_labels']=';'.join(T.get(i,{}).get('name',i) for i in sorted(inherit))
    write(a.out/'case_metadata.csv',rows)
    coverage=[]
    for corpus in sorted({r['corpus'] for r in rows}):
        cs=[r for r in rows if r['corpus']==corpus]
        coverage.append(dict(corpus=corpus,n_cases=len(cs),disease_inheritance=sum(bool(r['disease_inheritance_ids']) for r in cs),
                             protein_consequence=sum(bool(r['protein_consequences']) for r in cs),
                             complete_protein_class=sum(r['protein_class_complete']=='yes' for r in cs)))
    write(a.out/'coverage.csv',coverage)
    (a.out/'interpretation.json').write_text(json.dumps({'inheritance':'Inheritance annotations linked to the stated disease identifier; not a patient-level segregation inference.',
        'protein_classes':'Descriptive HGVS protein-consequence syntax; no inference about pathogenicity, NMD, gain/loss of function, or penetrance.',
        'classes':dict(collections.Counter(r['protein_consequences'] for r in rows if r['protein_consequences']))},indent=2)+'\n')

if __name__=='__main__':main()
