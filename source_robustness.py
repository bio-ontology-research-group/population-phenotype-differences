#!/usr/bin/env python3
"""Source-specific regional contrasts and paired bedside-restriction benchmarks.

Uses immutable v4 case/rank outputs. Null ranges from stratified block permutations
are explicitly distinct from publication-cluster bootstrap confidence intervals.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from population_model import group_of, load_groups
from symptom_sign import bh, contrast

POPS = ['EUR', 'MENA', 'SASIA', 'LATAM', 'AFR', 'EASIA']
CORPORA = ['store', 'pmcp', 'pubmed_abs']
TEST = {'laboratory', 'imaging', 'physiological_test', 'pathology'}
BED = {'clinical_exam', 'symptom_history'}


def read(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t' if str(path).endswith('.tsv') else ','))


def write(path, rows):
    if not rows:
        return
    with Path(path).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n'); w.writeheader(); w.writerows(rows)


def cluster_boot(values, clusters, rng, repetitions):
    """Resample clusters with replacement; retain case-weighted means."""
    values = np.asarray(values, float)
    _, ix = np.unique(clusters, return_inverse=True)
    n = np.bincount(ix)
    totals = np.array([np.bincount(ix, weights=values[:, j]) for j in range(values.shape[1])]).T
    if len(n) < 2:
        return np.full((2, values.shape[1]), np.nan), len(n)
    means = np.empty((repetitions, values.shape[1]))
    for b in range(repetitions):
        draw = rng.integers(0, len(n), size=len(n))
        means[b] = totals[draw].sum(axis=0) / n[draw].sum()
    return np.percentile(means, [2.5, 97.5], axis=0), len(n)


def load_cases(root):
    groups = load_groups(str(root / 'results/groups_lastauthor.tsv'))
    with (root / 'results/termclass.jsonl').open() as f:
        types = {d['id']: d['class'] for d in map(json.loads, f)}
    cases = []
    for name in ['cases_store.tsv', 'cases_pmcp_map2.tsv', 'cases_abs.tsv']:
        for r in read(root / 'results/v4/fine2' / name):
            if r['region'] not in POPS:
                continue
            ts = set(r['terms'].split(';')) - {''}
            r['n_test'] = sum(types.get(t) in TEST for t in ts)
            r['n_bed'] = sum(types.get(t) in BED for t in ts)
            r['n_terms'] = len(ts)
            r['share'] = r['n_test'] / (r['n_test'] + r['n_bed']) if r['n_test'] + r['n_bed'] else np.nan
            r['group'] = group_of(groups, r['pmid'])
            r['author'] = r.get('author_region') or 'UNKNOWN'
            cases.append(r)
    return cases


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--perms', type=int, default=2000)
    ap.add_argument('--boot', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=20261029)
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    cases = load_cases(a.root)
    assert len(cases) == len({r['id'] for r in cases}), 'Duplicate case IDs'
    meta = {r['id']: r for r in cases}
    rank_files = sorted((a.root / 'results/v4/gda').glob('gda_chunk*.csv'))
    ranks = [r for path in rank_files for r in read(path)]
    assert len(ranks) == len({r['id'] for r in ranks}), 'Duplicate ranks'
    for r in ranks:
        m = meta[r['id']]
        for k in ['pmid', 'group', 'author']:
            r[k] = m[k]
        for k in ['n_terms', 'n_test', 'n_bed']:
            r[k] = int(r[k])
        for k in ['rank_all', 'rank_bed', 'coverage']:
            r[k] = float(r[k]) if r[k] else np.nan
    denominators, paired = [], []
    rng = np.random.default_rng(a.seed)
    for corpus in CORPORA:
        for pop in POPS:
            cs = [r for r in cases if r['corpus'] == corpus and r['region'] == pop]
            rr = [r for r in ranks if r['corpus'] == corpus and r['region'] == pop]
            pr = [r for r in rr if r['n_bed'] > 0 and r['n_test'] > 0]
            denominators.append(dict(corpus=corpus, population=pop, cases=len(cs),
                                     publications=len({r['pmid'] for r in cs}),
                                     no_classified_findings=sum(np.isnan(r['share']) for r in cs),
                                     ranked=len(rr), paired=len(pr)))
            if not pr:
                continue
            vals = np.array([[np.log2(r['rank_bed']/r['rank_all']),
                              100*(r['rank_all'] <= 10), 100*(r['rank_bed'] <= 10),
                              100*((r['rank_bed'] <= 10)*1-(r['rank_all'] <= 10)*1)] for r in pr])
            means = vals.mean(axis=0)
            for unit in ['pmid', 'group']:
                ci, clusters = cluster_boot(vals, [r[unit] for r in pr], rng, a.boot)
                paired.append(dict(corpus=corpus, population=pop, resampling_unit=unit,
                                   n_cases=len(pr), n_clusters=clusters,
                                   n_with_other_terms=sum(r['n_terms'] > r['n_bed']+r['n_test'] for r in pr),
                                   log2_rank_loss=means[0], log2_lo=ci[0,0], log2_hi=ci[1,0],
                                   geometric_rank_ratio=2**means[0], ratio_lo=2**ci[0,0], ratio_hi=2**ci[1,0],
                                   top10_all=means[1], top10_bed=means[2], top10_change_pp=means[3],
                                   top10_change_lo=ci[0,3], top10_change_hi=ci[1,3]))
            print(corpus,pop,len(pr),f'rank ratio {2**means[0]:.3f}',flush=True)
    write(a.out/'denominators.csv',denominators)
    write(a.out/'paired_ranking.csv',paired)
    outputs = []
    tasks = [
        ('phenotype', [r for r in cases if np.isfinite(r['share'])],
         {'test_share_pp': lambda r: 100*r['share'], 'log1p_bed': lambda r: np.log1p(r['n_bed']),
          'log1p_test': lambda r: np.log1p(r['n_test'])}),
        ('gene_ranking', ranks,
         {'log2_rank_all': lambda r: np.log2(r['rank_all']), 'top10_all_pp': lambda r: 100*(r['rank_all']<=10),
          'coverage_pp': lambda r: 100*r['coverage']})]
    for family, data, outcomes in tasks:
        for version in ['total', 'author_adjusted']:
            key = (lambda r: (r['gene'],r['corpus'])) if version=='total' else (lambda r: (r['gene'],r['corpus'],r['author']))
            for corpus in CORPORA:
                sub = [r for r in data if r['corpus']==corpus]
                for pop in POPS[1:]:
                    result = contrast(sub,'region',pop,'EUR',key,outcomes,rng,a.perms)
                    if result is None:
                        continue
                    est, lo, hi, p, n, n1, ns = result
                    for j,outcome in enumerate(outcomes):
                        outputs.append(dict(family=family,version=version,corpus=corpus,population=pop,
                                            outcome=outcome,estimate=est[j],null_lo=lo[j],null_hi=hi[j],p=p[j],
                                            n_cases=n,n_population=n1,n_strata=ns))
                    print('contrasts',family,version,corpus,pop,flush=True)
    for family in ['phenotype','gene_ranking']:
        for version in ['total','author_adjusted']:
            sub=[r for r in outputs if r['family']==family and r['version']==version]
            for r,q in zip(sub,bh([r['p'] for r in sub])):
                r['q']=float(q);r['bh_family_size']=len(sub)
    write(a.out/'regional_contrasts.csv',outputs)
    inputs=[a.root/'results/termclass.jsonl',a.root/'results/groups_lastauthor.tsv',
            *sorted((a.root/'results/v4/fine2').glob('*.tsv')),*rank_files]
    manifest={'seed':a.seed,'permutations':a.perms,'bootstrap_replicates':a.boot,
              'inputs':{str(p.relative_to(a.root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
              'case_count':len(cases),'rank_count':len(ranks),'status':'complete'}
    (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    main()
