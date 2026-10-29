#!/usr/bin/env python3
"""Exploratory gene-specific phenotype contrasts with overlap/metadata audits.

Regional groups are observed labels, not verified genetic ancestry. The primary
family is all eligible gene x region comparisons, matched by corpus; additional
versions condition on author region, care country with stated origin, or disease
and allelic state. HPO tests use an is_a ancestor closure below organ-system level.
"""
import argparse,collections,csv,json
from pathlib import Path
import numpy as np
from source_robustness import load_cases,write,POPS
from population_model import make_plan,draw
from symptom_sign import bh
from hpo_map2 import load_obo
from phenotype_specificity import ancestors_fn


def eligible(cs,pop,key,min_cases,min_pubs):
    by=collections.defaultdict(set)
    for c in cs:by[key(c)].add(c['region'])
    kept=[c for c in cs if by[key(c)]=={pop,'EUR'}]
    counts={r:sum(c['region']==r for c in kept) for r in (pop,'EUR')}
    pubs={r:len({c['pmid'] for c in kept if c['region']==r}) for r in (pop,'EUR')}
    status='eligible' if min(counts.values())>=min_cases and min(pubs.values())>=min_pubs else 'insufficient_overlap_or_support'
    return kept,counts,pubs,status


def permutation(cs,pop,key,Y,rng,perms):
    keys={k:i for i,k in enumerate(sorted({key(c) for c in cs}))}
    sid=np.array([keys[key(c)] for c in cs]);ns=len(keys)
    labels=np.array([c['region'] for c in cs]);n=np.bincount(sid,minlength=ns)
    plan=make_plan(cs,range(len(cs)),lambda i:key(cs[i]),lambda i:cs[i]['group'],labels)
    def stat(lab):
        one=(lab==pop).astype(float);n1=np.bincount(sid,weights=one,minlength=ns);n0=n-n1
        w=n1*n0/n;norm=w.sum()
        if norm==0:return np.zeros(Y.shape[1])
        coeff=(one*w[sid]/np.maximum(n1[sid],1)-(1-one)*w[sid]/np.maximum(n0[sid],1))/norm
        return coeff@Y
    observed=stat(labels)
    null=np.array([stat(draw(plan,labels,rng)) for _ in range(perms)])
    center=null.mean(axis=0);ex=observed-center
    p=(1+(np.abs(null-center)>=np.abs(ex)).sum(axis=0))/(perms+1)
    lo,hi=np.percentile(null-center,[2.5,97.5],axis=0)
    return observed,ex,lo,hi,p


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--metadata',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--perms',type=int,default=4999);ap.add_argument('--min-cases',type=int,default=10)
    ap.add_argument('--min-publications',type=int,default=3);ap.add_argument('--seed',type=int,default=20261030)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    meta={r['id']:r for r in csv.DictReader(a.metadata.open())}
    cases=load_cases(a.root)
    for c in cases:
        c.update({k:v for k,v in meta[c['id']].items() if k not in ('region','corpus','gene','id','pmid')})
    # Test share has a defined denominator for all cases used below.
    cases=[c for c in cases if np.isfinite(c['share'])]
    T=load_obo(str(a.root/'data/hpo/hp.obo'));anc=ancestors_fn(T)
    root='HP:0000118';organs={t for t,v in T.items() if root in v['is_a']}
    for c in cases:
        closure=set().union(*(anc(t) for t in c['terms'].split(';') if t in T))
        c['closure']={t for t in closure if root in anc(t) and t!=root and t not in organs}
    by=collections.defaultdict(list)
    for c in cases:by[c['gene']].append(c)
    versions={
      'corpus':(lambda c:True,lambda c:(c['corpus'],)),
      'author_region':(lambda c:c['author'] not in ('UNKNOWN',''),lambda c:(c['corpus'],c['author'])),
      'care_stated_origin':(lambda c:bool(c['care_country'] and (c['ancestry'] or c['origin_country'])) and c['region_source']!='care',lambda c:(c['corpus'],c['care_country'])),
      'disease_allelic_state':(lambda c:bool(c['disease_ids'] and c['allelic_states']) and 'unknown' not in c['allelic_states'].lower(),lambda c:(c['corpus'],c['disease_ids'],c['allelic_states']))}
    support=[];out=[];rng=np.random.default_rng(a.seed)
    for version,(select,key) in versions.items():
        for gene,allcs in sorted(by.items()):
            for pop in POPS[1:]:
                cs=[c for c in allcs if c['region'] in (pop,'EUR') and select(c)]
                if not any(c['region']==pop for c in cs) or not any(c['region']=='EUR' for c in cs):continue
                cs,n,pubs,status=eligible(cs,pop,key,a.min_cases,a.min_publications)
                support.append(dict(version=version,gene=gene,population=pop,status=status,n_population=n[pop],n_eur=n['EUR'],publications_population=pubs[pop],publications_eur=pubs['EUR']))
                if status!='eligible':continue
                count=collections.Counter(t for c in cs for t in c['closure'])
                termpub=collections.defaultdict(set)
                for c in cs:
                    for t in c['closure']:termpub[t].add(c['pmid'])
                terms=sorted(t for t,nobs in count.items() if 5<=nobs<=len(cs)-5 and len(termpub[t])>=3)
                features=['test_share_pp','log1p_bed','log1p_test']+terms
                Y=np.array([[100*c['share'],np.log1p(c['n_bed']),np.log1p(c['n_test'])]+[100*(t in c['closure']) for t in terms] for c in cs])
                observed,ex,lo,hi,p=permutation(cs,pop,key,Y,rng,a.perms)
                for j,feature in enumerate(features):
                    out.append(dict(version=version,gene=gene,population=pop,outcome=feature,label=T.get(feature,{}).get('name',feature),
                                    family='term' if j>=3 else 'composition',observed=observed[j],estimate=ex[j],null_lo=lo[j],null_hi=hi[j],p=p[j],
                                    n_population=n[pop],n_eur=n['EUR'],publications_population=pubs[pop],publications_eur=pubs['EUR'],
                                    shared_publications=len({c['pmid'] for c in cs if c['region']==pop}&{c['pmid'] for c in cs if c['region']=='EUR'}),
                                    n_strata=len({key(c) for c in cs})))
                print(version,gene,pop,n[pop],n['EUR'],len(terms),'terms',flush=True)
        # Write each completed version so progress is recoverable.
        for fam in ['composition','term']:
            sub=[r for r in out if r['version']==version and r['family']==fam]
            for r,q in zip(sub,bh([r['p'] for r in sub])):r['q']=float(q);r['bh_family_size']=len(sub)
        write(a.out/'support.csv',support);write(a.out/'contrasts.csv',out)
    summary={'permutations':a.perms,'seed':a.seed,'min_cases_per_group':a.min_cases,'min_publications_per_group':a.min_publications,
             'eligible_comparisons':dict(collections.Counter(r['version'] for r in support if r['status']=='eligible')),
             'significant_tests':dict(collections.Counter(r['version']+':'+r['family'] for r in out if r['q']<.05)),
             'closure':'asserted is_a from frozen OBO, excludes root and direct organ-system classes',
             'status':'complete','scope':'Exploratory associations; care-based regional labels excluded from the care-conditioned version. Metadata are not newly human-validated.'}
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__=='__main__':main()
