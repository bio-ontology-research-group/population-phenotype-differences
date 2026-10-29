#!/usr/bin/env python3
"""Exploratory, leakage-filtered disease-profile frequency benchmark.

The scoring specification was frozen in docs/frequency_evaluation_plan.md before
evaluation. Source-checked annotations are not independently clinician-validated.
"""
import argparse
import collections
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from rank_engine import Engine
from source_robustness import write

CATEGORIES = {'HP:0040280': (1., 1.), 'HP:0040281': (.8, .99),
              'HP:0040282': (.3, .79), 'HP:0040283': (.05, .29),
              'HP:0040284': (.01, .04), 'HP:0040285': (0., 0.)}


def read(path, delimiter=','):
    with path.open() as f:
        return list(csv.DictReader((s for s in f if not s.startswith('#')), delimiter=delimiter))


def excluded_reference(reference, heldout, groups, heldout_groups):
    pmids = set(re.findall(r'PMID:(\d+)', reference))
    return bool(pmids & heldout or any(groups.get(p, 'P:' + p) in heldout_groups for p in pmids))


def frequency_value(raw, endpoint='mid'):
    if raw in CATEGORIES:
        lo, hi = CATEGORIES[raw]
        return {'low': lo, 'high': hi, 'mid': (lo + hi) / 2}[endpoint]
    if re.fullmatch(r'\d+/\d+', raw):
        k, n = map(int, raw.split('/'))
        return (k + .5) / (n + 1) if 0 <= k <= n and n else None
    if raw.endswith('%'):
        try:
            v = float(raw[:-1]) / 100
            return v if 0 <= v <= 1 else None
        except ValueError:
            return None
    return None


def resolved_weight(values, missing=1., endpoint='mid'):
    # Missing rows do not override an explicit value, but conflicting explicit
    # count/category/percentage assertions remain unresolved.
    valid = {v for v in values if frequency_value(v) is not None}
    return frequency_value(next(iter(valid)), endpoint) if len(valid) == 1 else missing


def normalize(matrix):
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise ValueError('A disease profile has no positive frequency weight')
    return sp.diags(1. / totals) @ matrix


def rank_interval(scores, index):
    score = scores[index]
    tied = np.isclose(scores, score, rtol=0, atol=1e-10)
    best = int(np.sum((scores > score) & ~tied)) + 1
    worst = best + int(tied.sum()) - 1
    return best, worst, (best + worst) / 2.


def curated_counts(annotations, assessments, metadata, family=False, truncating=False, omit=None):
    meta = {(r['cohort_id'], r['patient_label']): r for r in metadata}
    cohort_disease = {r['cohort_id']: r['disease_id'] for r in annotations}
    members = collections.Counter(r['family_id'] for r in metadata)
    counts = collections.defaultdict(lambda: [0., 0.])
    for r in assessments:
        if r['pmid'] == omit or r['status'] not in ('present', 'absent'):
            continue
        if truncating and r['pmid'] == '25896430' and r['patient_label'] not in ('Patient 4', 'Patient 6', 'Patient 8'):
            continue
        m = meta[r['cohort_id'], r['patient_label']]
        w = 1. / members[m['family_id']] if family else 1.
        pair = cohort_disease[r['cohort_id']], r['hpo_id']
        counts[pair][0] += w * (r['status'] == 'present')
        counts[pair][1] += w
    return dict(counts)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args(); root = a.root; a.out.mkdir(parents=True, exist_ok=True)
    v5 = root / 'results/v5'; curation = v5 / 'frequency_curation_combined'
    ann = read(curation / 'frequency_annotations.csv')
    assessments = read(curation / 'patient_assessments.csv')
    metadata = read(curation / 'patient_metadata.csv')
    queries = read(v5 / 'clinical_evaluation_audit/query_eligibility.csv')
    groups = dict(csv.reader((root / 'results/groups_lastauthor.tsv').open(), delimiter='\t'))
    heldout = {r['pmid'] for r in queries}
    heldout_groups = {groups.get(p, 'P:' + p) for p in heldout}
    overlap = []
    for pmid in sorted({r['pmid'] for r in ann}):
        source = [r for r in ann if r['pmid'] == pmid][0]
        gene = source['gene']
        # A distinct stated care country, cohort and causal-variant description
        # supports separation; text-level variant screening is not HGVS equivalence.
        variants = {v for r in metadata if r['pmid'] == pmid for v in r['variants_as_reported'].split(';')}
        target_variants = {r['source_variant_cdna'] for r in queries if r['gene'] == gene and r['eligible_primary'] == 'True'}
        hits = sorted(variants & target_variants)
        assert pmid not in heldout and groups.get(pmid, 'P:' + pmid) not in heldout_groups and not hits
        overlap.append(dict(pmid=pmid, gene=gene, care_country=source['care_country'],
                            reporting_group=groups.get(pmid, 'P:' + pmid), shared_query_publication=False,
                            shared_query_reporting_group=False, exact_reported_cdna_overlap=';'.join(hits),
                            disposition='no_overlap_identified_from_source_cohort_country_and_reported_variants',
                            limitation='No private patient identifiers or cross-transcript variant normalization; US within-gene identity ambiguity retained.'))
    write(a.out / 'curated_query_overlap.csv', overlap)
    eng = Engine(str(root / 'data/hpo/hp.obo'), str(root / 'data/hpo/genes_to_phenotype.txt'))
    disease_genes = collections.defaultdict(set)
    for r in read(root / 'data/hpo/genes_to_phenotype.txt', '\t'):
        if r['disease_id'].startswith('OMIM:'):
            disease_genes[r['disease_id']].add(r['gene_symbol'])
    removed_kdm6a = 'KDM6A' in disease_genes['OMIM:147920']
    disease_genes['OMIM:147920'] = {'KMT2D'}
    support = collections.defaultdict(set); frequencies = collections.defaultdict(set)
    audit = collections.Counter(); exclusions = []
    for r in read(root / 'data/hpo/phenotype.hpoa', '\t'):
        d, h = r['database_id'], eng.alt.get(r['hpo_id'], r['hpo_id'])
        if d not in disease_genes or r['aspect'] != 'P' or 'NOT' in r['qualifier'] or r['frequency'] == 'HP:0040285':
            continue
        if h not in eng.tid or eng.terms[h]['obsolete'] or 'HP:0000118' not in eng.anc_sets[h]:
            continue
        audit['eligible_before_reference_filter'] += 1
        if excluded_reference(r['reference'], heldout, groups, heldout_groups):
            audit['excluded_reference_rows'] += 1
            exclusions.append({k:r[k] for k in ('database_id','hpo_id','reference','frequency')})
            continue
        support[d].add(h)
        if not any(r[k] for k in ('onset', 'sex', 'modifier')):
            frequencies[d, h].add(r['frequency'])
        audit['retained_rows'] += 1
    write(a.out / 'excluded_reference_annotations.csv', exclusions)
    # IC is calculated on the union of filtered disease terms per candidate gene.
    g2t = collections.defaultdict(set)
    for d, terms in support.items():
        for gene in disease_genes[d]: g2t[gene].update(terms)
    candidate_genes = sorted(g2t); gid = {g:i for i,g in enumerate(candidate_genes)}
    eng.genes = candidate_genes; eng.gid = gid
    eng.build_gene_matrices(g2t); eng.set_ic('extrinsic')
    fixed_ic = eng.ic.copy()
    # Reuse the engine on disease profiles; keep the gene-derived IC fixed.
    diseases = sorted(support); eng.genes = diseases; eng.gid = {d:i for i,d in enumerate(diseases)}
    edges = [(eng.gid[d], gid[g]) for d in diseases for g in sorted(disease_genes[d])]
    di, gi = np.array(edges).T
    counts = curated_counts(ann, assessments, metadata)
    count_variants = {'standard': counts,
                      'family': curated_counts(ann, assessments, metadata, family=True),
                      'truncating': curated_counts(ann, assessments, metadata, truncating=True)}
    for pmid in sorted({r['pmid'] for r in ann}):
        count_variants['omit_' + pmid] = curated_counts(ann, assessments, metadata, omit=pmid)
    updates = []
    for (d,h),(k,n) in sorted(counts.items()):
        updates.append(dict(disease_id=d,hpo_id=h,affected=k,assessed=n,weight=(k+.5)/(n+1),
                            fixed_support=h in support[d],new_positive_association=h not in support[d] and k>0,
                            existing_values=';'.join(sorted(frequencies[d,h])),
                            original_weight=resolved_weight(frequencies[d,h])))
    write(a.out / 'frequency_updates.csv', updates)
    expanded = {d:set(ts) for d,ts in support.items()}
    for (d,h),(k,n) in counts.items():
        if k > 0: expanded[d].add(h)
    specifications = [dict(name='unweighted', weighted=False),dict(name='existing', weighted=True),
                      dict(name='updated', weighted=True,update='standard'),
                      dict(name='association_only', weighted=False,expand=True),
                      dict(name='association_and_frequency', weighted=True,expand=True,update='standard'),
                      dict(name='updated_family',weighted=True,update='family'),
                      dict(name='updated_truncating',weighted=True,update='truncating')]
    for endpoint in ['low','high']:
        for arm in ['existing','updated']:
            specifications.append(dict(name=arm+'_category_'+endpoint,weighted=True,endpoint=endpoint,
                                       update='standard' if arm=='updated' else None))
    for arm in ['existing','updated']:
        specifications.append(dict(name=arm+'_missing_half',weighted=True,missing=.5,
                                   update='standard' if arm=='updated' else None))
    for pmid in sorted({r['pmid'] for r in ann}):
        specifications.append(dict(name='updated_omit_'+pmid,weighted=True,update='omit_'+pmid))
    index = {r['id']:r for r in read(root/'results/phenopackets/index.tsv','\t')}
    queries = [r for r in queries if r['eligible_primary']=='True']
    query_terms = {}
    for r in queries:
        packet=json.loads((root/'results/phenopackets'/index[r['id']]['file']).read_text())
        query_terms[r['id']]=eng.map_terms([f['type']['id'] for f in packet['phenotypicFeatures'] if not f.get('excluded',False)])
        assert query_terms[r['id']] and r['gene'] in gid
    ranks=[]
    for spec in specifications:
        profiles=expanded if spec.get('expand') else support
        eng.build_gene_matrices(profiles)
        assert np.array_equal(eng.ic,fixed_ic)
        rr=[];cc=[];vv=[]
        for d in diseases:
            for h in sorted(profiles[d]):
                w=1.
                if spec['weighted']:
                    w=resolved_weight(frequencies[d,h],spec.get('missing',1.),spec.get('endpoint','mid'))
                    if spec.get('update') and (d,h) in count_variants[spec['update']]:
                        k,n=count_variants[spec['update']][d,h];w=(k+.5)/(n+1)
                rr.append(eng.gid[d]);cc.append(eng.tid[h]);vv.append(w)
        weighted=normalize(sp.csr_matrix((vv,(rr,cc)),shape=eng.GT.shape))
        for r in queries:
            dscores=eng.bma_all_genes(query_terms[r['id']],GTnorm=weighted)
            gscores=np.full(len(candidate_genes),-np.inf);np.maximum.at(gscores,gi,dscores[di])
            best,worst,mid=rank_interval(gscores,gid[r['gene']])
            target_profiles=[d for d in diseases if r['gene'] in disease_genes[d]]
            best_profiles=[d for d in target_profiles if np.isclose(dscores[eng.gid[d]],gscores[gid[r['gene']]],rtol=0,atol=1e-10)]
            ranks.append(dict(id=r['id'],pmid=r['pmid'],gene=r['gene'],corpus=r['corpus'],
                              stringent_identity=r['stringent_identity_sensitivity'],arm=spec['name'],
                              n_query_terms=len(query_terms[r['id']]),rank=mid,best_tied_rank=best,worst_tied_rank=worst,
                              reciprocal_rank=1/mid,top1=mid<=1,top10=mid<=10,
                              highest_scoring_target_profiles=';'.join(best_profiles),
                              target_gene_score=gscores[gid[r['gene']]],
                              intended_disorder_score=dscores[eng.gid[r['benchmark_disorder']]]))
        print(spec['name'],'complete',flush=True)
    write(a.out/'query_ranks.csv',ranks)
    lookup={(r['id'],r['arm']):r for r in ranks};metrics=[]
    for subset in ['all','stringent_identity','pavs_clin','pavs_us','KMT2D','MPV17','RPGRIP1']:
        eligible=[q for q in queries if subset=='all' or (subset=='stringent_identity' and q['stringent_identity_sensitivity']=='True') or q['corpus']==subset or q['gene']==subset]
        for spec in specifications:
            rs=[lookup[q['id'],spec['name']] for q in eligible]
            base=[lookup[q['id'],'existing'] for q in eligible]
            metrics.append(dict(subset=subset,arm=spec['name'],n=len(rs),n_publications=len({r['pmid'] for r in rs}),
                                mrr=np.mean([r['reciprocal_rank'] for r in rs]),top1=np.mean([r['top1'] for r in rs]),
                                top10=np.mean([r['top10'] for r in rs]),median_rank=np.median([r['rank'] for r in rs]),
                                improved_vs_existing=sum(r['rank']<b['rank'] for r,b in zip(rs,base)),
                                worsened_vs_existing=sum(r['rank']>b['rank'] for r,b in zip(rs,base))))
    write(a.out/'metrics.csv',metrics)
    # Publication-cluster uncertainty acknowledges that only three source
    # publications contribute queries; gene/country strata are even smaller.
    rng=np.random.default_rng(20261029);paired=[]
    pubs=sorted({q['pmid'] for q in queries});draws=rng.integers(len(pubs),size=(5000,len(pubs)))
    for arm,baseline in [('existing','unweighted'),('updated','existing'),('association_only','unweighted'),('association_and_frequency','association_only')]:
        for metric in ['reciprocal_rank','top10']:
            sums=[];sizes=[]
            for pub in pubs:
                qs=[q for q in queries if q['pmid']==pub]
                sums.append(sum(float(lookup[q['id'],arm][metric])-float(lookup[q['id'],baseline][metric]) for q in qs));sizes.append(len(qs))
            vals=np.array(sums)[draws].sum(axis=1)/np.array(sizes)[draws].sum(axis=1)
            lo,hi=np.quantile(vals,[.025,.975])
            paired.append(dict(arm=arm,baseline=baseline,metric=metric,change=sum(sums)/sum(sizes),
                               ci_low=lo,ci_high=hi,n=len(queries),n_publications=len(pubs),bootstrap=5000))
    write(a.out/'paired_changes.csv',paired)
    audit.update(n_genes=len(candidate_genes),n_disease_profiles=len(diseases),n_queries=len(queries),
                 curated_source_rows=len(ann),curated_unique_disease_terms=len(counts),
                 fixed_support_updates=sum(r['fixed_support'] for r in updates),
                 new_positive_associations=sum(r['new_positive_association'] for r in updates))
    summary=dict(counts=dict(audit),heldout_pmids=sorted(heldout),heldout_reporting_groups=sorted(heldout_groups),
                 kabuki1_kdm6a_mapping_removed=removed_kdm6a,specifications=specifications,
                 status='exploratory_source_checked_clinician_review_pending',
                 limitations=['Source-level variant/country checks found no curated-query overlap; private identity resolution unavailable.',
                              'Only three query publications and three reporting groups; cluster intervals are descriptive and unstable.',
                              'Curated cohorts have different ages and ascertainment; pooled weights are exploratory, with source-deletion sensitivities.',
                              'Within each curated disorder, population/care country has no contrasting group; population effects are not identifiable.',
                              'Novel support uses at least one explicit positive; all-zero annotations can change existing weights but never add positive associations.',
                              'This evaluates this weighted BMA rule, not all frequency-aware diagnostic methods.'],
                 hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                         [Path(__file__).resolve(),root/'docs/frequency_evaluation_plan.md',root/'data/hpo/phenotype.hpoa',root/'data/hpo/genes_to_phenotype.txt',
                          curation/'frequency_annotations.csv',curation/'patient_assessments.csv',v5/'clinical_evaluation_audit/query_eligibility.csv']})
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary['counts'],indent=2))


if __name__=='__main__':main()
