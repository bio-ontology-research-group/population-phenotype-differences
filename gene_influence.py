#!/usr/bin/env python3
"""Gene contribution and exact-variant overlap audits for regional comparisons."""
import argparse,collections,csv,json
from pathlib import Path
import numpy as np
from source_robustness import load_cases,write,POPS
from gene_regional_analysis import eligible,permutation
from symptom_sign import bh


def signature(s):
    fields=s.split(';') if s else []
    for syntax in ['hgvs.g:','hgvs.c:']:
        selected=sorted(set(v for v in fields if v.startswith(syntax)))
        if selected:return '|'.join(selected)
    return ''


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--metadata',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    cases=[r for r in load_cases(a.root) if np.isfinite(r['share'])]
    influence=[]
    for pop in POPS[1:]:
        strata=collections.defaultdict(list)
        for c in cases:
            if c['region'] in (pop,'EUR'):strata[c['gene'],c['corpus']].append(c)
        parts=collections.defaultdict(lambda:[0.,0.,0,0])
        for (gene,corpus),cs in strata.items():
            one=[c for c in cs if c['region']==pop];zero=[c for c in cs if c['region']=='EUR']
            if not one or not zero:continue
            w=len(one)*len(zero)/(len(one)+len(zero))
            d=100*(np.mean([c['share'] for c in one])-np.mean([c['share'] for c in zero]))
            parts[gene][0]+=w;parts[gene][1]+=w*d;parts[gene][2]+=len(one);parts[gene][3]+=len(zero)
        totalw=sum(x[0] for x in parts.values());num=sum(x[1] for x in parts.values());full=num/totalw
        for gene,(w,n,n1,n0) in parts.items():
            leave=(num-n)/(totalw-w) if w<totalw else np.nan
            influence.append(dict(population=pop,gene=gene,n_population=n1,n_eur=n0,weight_share=w/totalw,
                                  gene_difference_pp=n/w,contribution_pp=n/totalw,full_difference_pp=full,
                                  leave_gene_out_pp=leave,change_on_removal_pp=leave-full))
    write(a.out/'leave_gene_out.csv',influence)
    meta={r['id']:r for r in csv.DictReader(a.metadata.open())}
    variants=collections.defaultdict(list)
    for c in cases:
        m=meta[c['id']];sig=signature(m['variant_expressions'])
        if sig and m['allelic_states'] and 'unknown' not in m['allelic_states'].lower():
            c['variant_signature']=sig;c['allelic_states']=m['allelic_states'];variants[c['gene'],sig,m['allelic_states']].append(c)
    overlaps=[];contrast_rows=[];rng=np.random.default_rng(20261031)
    for (gene,sig,allelic),cs in variants.items():
        if len({c['region'] for c in cs})<2:continue
        overlaps.append(dict(gene=gene,variant_signature=sig,allelic_state=allelic,n_cases=len(cs),n_publications=len({c['pmid'] for c in cs}),
                             region_counts=json.dumps(dict(collections.Counter(c['region'] for c in cs)),sort_keys=True)))
    # Preserve the same minimum support as gene analyses; do not relax after seeing overlap.
    by=collections.defaultdict(list)
    for cs in variants.values():
        for c in cs:by[c['gene']].append(c)
    support=[]
    for gene,cs in sorted(by.items()):
        for pop in POPS[1:]:
            pool=[c for c in cs if c['region'] in (pop,'EUR')]
            key=lambda c:(c['corpus'],c['variant_signature'],c['allelic_states'])
            kept,n,pubs,status=eligible(pool,pop,key,10,3)
            support.append(dict(gene=gene,population=pop,status=status,n_population=n[pop],n_eur=n['EUR'],publications_population=pubs[pop],publications_eur=pubs['EUR']))
            if status=='eligible':
                Y=np.array([[100*c['share']] for c in kept])
                obs,ex,lo,hi,p=permutation(kept,pop,key,Y,rng,4999)
                contrast_rows.append(dict(gene=gene,population=pop,observed=obs[0],estimate=ex[0],null_lo=lo[0],null_hi=hi[0],p=p[0],n_population=n[pop],n_eur=n['EUR']))
    for r,q in zip(contrast_rows,bh([r['p'] for r in contrast_rows])):r['q']=q
    write(a.out/'shared_variants.csv',overlaps);write(a.out/'variant_support.csv',support);write(a.out/'variant_contrasts.csv',contrast_rows)
    summary={'shared_variant_allelic_profiles':len(overlaps),'eligible_variant_matched_comparisons':len(contrast_rows),
             'variant_matching':'Exact hgvs.g expressions preferred, otherwise exact hgvs.c; same corpus and allelic state; no equivalence normalization across transcripts/builds.',
             'influence':'Uncentered Mantel-Haenszel contrasts with cases lacking classified findings excluded. Descriptive deletion analysis; no new p-values.',
             'max_gene_influence':{p:max((dict(gene=r['gene'],change_pp=r['change_on_removal_pp']) for r in influence if r['population']==p),key=lambda r:abs(r['change_pp'])) for p in POPS[1:]}}
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
