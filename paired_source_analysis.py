#!/usr/bin/env python3
"""Describe paired publication-level extraction outputs; human accuracy is separate.

The one-individual sensitivity is selected by output cardinality, so it is not a
human-validated patient match. Raw outputs are never edited to remove errors.
"""
import argparse,collections,csv,hashlib,json
from pathlib import Path
import numpy as np
from source_robustness import BED,TEST,POPS,write
from symptom_sign import bh
from hpo_map2 import load_obo
from phenotype_specificity import ancestors_fn

def successful(path,key):
    ans={}
    for l in path.open():
        d=json.loads(l)
        if d.get('ok'):ans[d[key]]=d
    return ans

def metrics(terms,classes):
    nt=sum(classes.get(t) in TEST for t in terms)
    nb=sum(classes.get(t) in BED for t in terms)
    return dict(n_terms=len(terms),n_test=nt,n_bed=nb,
                n_mixed=sum(classes.get(t)=='mixed' for t in terms),
                n_unclassified=sum(t not in classes for t in terms),
                test_share=nt/(nt+nb) if nt+nb else np.nan)

def paired_intervals(values,rng,boot):
    v=np.asarray(values,float);v=v[np.isfinite(v)]
    if not len(v):return len(v),np.nan,np.nan,np.nan
    if len(v)==1:return 1,float(v[0]),np.nan,np.nan
    sampled=v[rng.integers(0,len(v),size=(boot,len(v)))].mean(axis=1)
    lo,hi=np.percentile(sampled,[2.5,97.5])
    return len(v),float(v.mean()),float(lo),float(hi)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--termclass-extra',type=Path)
    ap.add_argument('--boot',type=int,default=5000)
    ap.add_argument('--seed',type=int,default=20261031)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    private=a.root/'data/v5_private';mapped=private/'paired_mapping_final'
    docs=successful(private/'paired_extraction.jsonl','key')
    selected=successful(mapped/'selected.jsonl','key')
    phrases={json.loads(l)['key'] for l in (mapped/'phrases.jsonl').open()}
    if phrases-set(selected):raise ValueError(f'{len(phrases-set(selected))} phrases still require mapping')
    classes={d['id']:d['class'] for d in map(json.loads,(a.root/'results/termclass.jsonl').open()) if d.get('ok')}
    if a.termclass_extra:
        for d in map(json.loads,a.termclass_extra.open()):
            if d.get('ok') and d['class'] in BED|TEST|{'mixed'}:
                if d['id'] in classes and classes[d['id']]!=d['class']:raise ValueError('Changed frozen class')
                classes[d['id']]=d['class']
    T=load_obo(str(a.root/'data/hpo/hp.obo'));anc=ancestors_fn(T)
    mentions=[json.loads(l) for l in (mapped/'mentions.jsonl').open()]
    termsets=collections.defaultdict(set);n_mentions=collections.Counter();nulls=collections.Counter()
    for m in mentions:
        if m['status']!='present':continue
        n_mentions[m['document_key']]+=1
        h=selected[m['phrase_key']]['hpo']
        if h:termsets[m['document_key']].add(h)
        else:nulls[m['document_key']]+=1
    manifest={r['pmid']:r for r in csv.DictReader((a.root/'results/v5/paired_sources/manifest.csv').open()) if r['status']=='ready'}
    rows=[]
    for pmid,r in sorted(manifest.items()):
        ds={s:docs[f"{pmid}:{r['gene']}:{s}"] for s in ['abstract','fulltext']}
        row=dict(pmid=pmid,gene=r['gene'],region=r['region'])
        one=all(len(d['result']['patients'])==1 and d['result']['patients'][0].get('unit')=='individual' for d in ds.values())
        row['single_individual_both']=int(one)
        for s,d in ds.items():
            ts=termsets[d['key']]
            row.update({f'{s}_{k}':v for k,v in metrics(ts,classes).items()})
            row[f'{s}_patient_entries']=len(d['result']['patients'])
            row[f'{s}_aggregate_entries']=sum(p.get('unit')=='aggregate' for p in d['result']['patients'])
            row[f'{s}_present_mentions']=n_mentions[d['key']]
            row[f'{s}_null_mappings']=nulls[d['key']]
        ats=termsets[ds['abstract']['key']];fts=termsets[ds['fulltext']['key']]
        row['exact_jaccard']=len(ats&fts)/len(ats|fts) if ats|fts else np.nan
        row['abstract_terms_retained_exact']=len(ats&fts)/len(ats) if ats else np.nan
        fc=set().union(*(anc(t) for t in fts))|fts
        row['abstract_terms_retained_ancestor']=len(ats&fc)/len(ats) if ats else np.nan
        for k in ['n_terms','n_test','n_bed','test_share']:
            row[f'delta_{k}']=row[f'fulltext_{k}']-row[f'abstract_{k}']
        rows.append(row)
    write(a.out/'publication_pairs.csv',rows)
    rng=np.random.default_rng(a.seed);summary=[];interactions=[]
    outcomes=['delta_n_terms','delta_n_test','delta_n_bed','delta_test_share',
              'exact_jaccard','abstract_terms_retained_exact','abstract_terms_retained_ancestor']
    for subset in ['all_publications','one_individual_each_output']:
        pool=[r for r in rows if subset=='all_publications' or r['single_individual_both']]
        for region in ['ALL']+POPS:
            sub=[r for r in pool if region=='ALL' or r['region']==region]
            for outcome in outcomes:
                n,mean,lo,hi=paired_intervals([r[outcome] for r in sub],rng,a.boot)
                summary.append(dict(subset=subset,region=region,outcome=outcome,publications=n,mean=mean,lo=lo,hi=hi))
        # Exploratory regional source interaction. Publication-label permutation;
        # require five finite paired observations per group, no threshold relaxation.
        for region in POPS[1:]:
            for outcome in outcomes[:4]:
                r1=np.array([r[outcome] for r in pool if r['region']==region and np.isfinite(r[outcome])])
                r0=np.array([r[outcome] for r in pool if r['region']=='EUR' and np.isfinite(r[outcome])])
                if min(len(r1),len(r0))<5:continue
                observed=r1.mean()-r0.mean();v=np.r_[r1,r0];null=[]
                for _ in range(4999):
                    z=rng.permutation(v);null.append(z[:len(r1)].mean()-z[len(r1):].mean())
                p=(1+sum(abs(x)>=abs(observed) for x in null))/5000
                interactions.append(dict(subset=subset,region=region,outcome=outcome,n_region=len(r1),n_EUR=len(r0),
                                         difference=observed,p=p))
    for subset in {r['subset'] for r in interactions}:
        sub=[r for r in interactions if r['subset']==subset]
        for r,q in zip(sub,bh([r['p'] for r in sub])):r['q']=float(q)
    write(a.out/'paired_summary.csv',summary);write(a.out/'region_interactions.csv',interactions)
    missing=sorted(set().union(*termsets.values())-set(classes))
    write(a.out/'unclassified_terms.csv',[dict(id=t,name=T.get(t,{}).get('name','')) for t in missing])
    if not missing and (a.out/'unclassified_terms.csv').exists():
        (a.out/'unclassified_terms.csv').unlink()
    qc=dict(paired_publications=len(rows),single_individual_both=sum(r['single_individual_both'] for r in rows),
            extraction_documents=len(docs),phrase_mappings=len(phrases),unclassified_terms=len(missing),
            bootstrap_replicates=a.boot,seed=a.seed,human_validation='pending',
            estimand='Difference between extracted publication-level phenotype unions. The single-individual subset is output-based, with identity validation pending.',
            interaction_scope='Unadjusted exploratory differences of within-publication source effects across regions; gene and recruitment distributions differ.',
            hypotheses='BH across all eligible regional interaction tests separately per subset',
            hashes={str(p.relative_to(a.root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [private/'paired_extraction.jsonl',mapped/'selected.jsonl',a.root/'results/termclass.jsonl']})
    (a.out/'summary.json').write_text(json.dumps(qc,indent=2)+'\n')
    print(json.dumps({k:v for k,v in qc.items() if k!='hashes'},indent=2))

if __name__=='__main__':main()
