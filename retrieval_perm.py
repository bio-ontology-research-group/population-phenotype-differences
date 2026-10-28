#!/usr/bin/env python3
"""Permutation test for the population-matching interaction in nearest-case retrieval.

For each query case (region r1 or r2) the query's gene is represented by k cases from
the query's own region ("same") or from the other region ("other"); every other gene is
scored by its best-matching library case. d = log2(E[rank_other] / E[rank_same]), with the
expectation taken exactly over random k-subsets of each library.

Statistic (gene-matched interaction): over genes with queries from both regions,
mean_g [(mean d over r1 queries + mean d over r2 queries) / 2]. Positive = cases are
retrieved better by cases of their own region, after removing any main effect of one
region's cases being better prototypes for everybody.

Nulls permute region labels between publications (all cases of a publication move
together) within gene ("gene"), or within (gene, author region) ("gene+author", the
patient-region effect beyond the reporting authors' region).

Usage: retrieval_perm.py --cases T.tsv --hpo-dir D --cache S.npy [--authors A.tsv]
       --r1 MENA --r2 EUR --perms 1000 --out PREFIX
"""
import argparse, collections, csv
from math import comb
import numpy as np
import case_pairs as cp


def exp_rank(vals, k, srt):
    """Expected rank of the max of a random k-subset of vals against sorted competitor scores."""
    v = np.sort(vals)[::-1]; n = len(v)
    p = np.array([comb(n - m, k - 1) for m in range(1, n - k + 2)], float) / comb(n, k)
    r = 1 + len(srt) - np.searchsorted(srt, v[:n - k + 1], side="right")
    return float(p @ r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", action="append", required=True); ap.add_argument("--hpo-dir", required=True)
    ap.add_argument("--cache"); ap.add_argument("--authors"); ap.add_argument("--sim", choices=["z", "raw"], default="z")
    ap.add_argument("--r1", default="MENA"); ap.add_argument("--r2", default="EUR")
    ap.add_argument("--k", type=int, default=5); ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--agg", choices=["mean", "median"], default="mean", help="aggregate over genes")
    ap.add_argument("--query-corpus", help="comma list: only these corpora as queries")
    ap.add_argument("--lib-corpus", help="comma list: only these corpora as same-gene library cases")
    ap.add_argument("--groups", help="pmid -> reporting group TSV; permute labels between groups instead of publications")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sim, cases = cp.setup(a)
    GRP = {r[0]: r[1] for r in csv.reader(open(a.groups), delimiter="\t")} if a.groups else {}
    rng = np.random.default_rng(7)
    genes = sorted({c["gene"] for c in cases}); gidx = {g: i for i, g in enumerate(genes)}
    cg = np.array([gidx[c["gene"]] for c in cases]); pm = np.array([c["pmid"] for c in cases])
    use = [i for i, c in enumerate(cases) if c["region"] in (a.r1, a.r2)]
    qc = set(a.query_corpus.split(",")) if a.query_corpus else None
    lc = set(a.lib_corpus.split(",")) if a.lib_corpus else None
    bygene = collections.defaultdict(list)
    for i in use: bygene[cases[i]["gene"]].append(i)
    # per query: competitor scores (best per other gene) and own-gene candidate library
    Q = {}
    for g, idx in bygene.items():
        if len(idx) < 3: continue
        for q in idx:
            if qc and cases[q]["corpus"] not in qc: continue
            s = sim.sim_row(q)
            valid = np.array([not cp.related(p, pm[q]) for p in pm])
            best = np.full(len(genes), -np.inf); np.maximum.at(best, cg[valid], s[valid])
            best = np.delete(best, gidx[g]); best = np.sort(best[np.isfinite(best)])
            cand = [j for j in idx if valid[j] and (not lc or cases[j]["corpus"] in lc)]
            if cand: Q[q] = (best, np.array(cand), s[cand])
    print(f"queries with a library: {len(Q)}", flush=True)
    memo = {}

    def stat(lab):
        per = collections.defaultdict(lambda: {a.r1: [], a.r2: []})
        for q, (best, cand, sv) in Q.items():
            lq = lab[q]; lc = np.array([lab[j] for j in cand])
            same, oth = sv[lc == lq], sv[lc != lq]
            k = min(len(same), len(oth), a.k)
            if k < 1: continue
            key = (q, k, same.tobytes(), oth.tobytes())
            if key not in memo:
                memo[key] = np.log2(exp_rank(oth, k, best)) - np.log2(exp_rank(same, k, best))
            per[cases[q]["gene"]][lq].append(memo[key])
        agg = np.median if a.agg == "median" else np.mean
        both = [(np.mean(v[a.r1]) + np.mean(v[a.r2])) / 2 for v in per.values() if v[a.r1] and v[a.r2]]
        if not both:  # queries from one region only (e.g. clinical queries): gene-averaged d of those queries
            one = [np.mean(v[a.r1] or v[a.r2]) for v in per.values() if v[a.r1] or v[a.r2]]
            return (agg(one) if one else np.nan), len(one)
        return (agg(both) if both else np.nan), len(both)

    lab0 = {i: cases[i]["region"] for i in use}
    t0, ng = stat(lab0)
    print(f"{a.r1} vs {a.r2}: gene-matched interaction {t0:+.3f} over {ng} genes", flush=True)
    pubs = collections.defaultdict(list)  # stratum -> list of (pmid, region)
    for i in use:
        unit = GRP.get(cases[i]["pmid"], cases[i]["pmid"]) if GRP else cases[i]["pmid"]
        pubs[(cases[i]["gene"], cases[i].get("corpus", "")), cases[i]["author"], (unit, cases[i]["region"])].append(i)
    out = []
    for null in ("gene", "gene+author"):
        strata = collections.defaultdict(list)
        for (g, au, p), idx in pubs.items():
            strata[(g,) if null == "gene" else (g, au)].append(idx)
        T = []
        for b in range(a.perms):
            lab = dict(lab0)
            for groups in strata.values():
                if len(groups) < 2: continue
                regs = [lab0[idx[0]] for idx in groups]; rng.shuffle(regs)
                for idx, r in zip(groups, regs):
                    for i in idx: lab[i] = r
            T.append(stat(lab)[0])
        T = np.array(T); T = T[np.isfinite(T)]
        p = (1 + np.sum(T - T.mean() >= t0 - T.mean())) / (1 + len(T))
        print(f"  null {null:12s}: mean {T.mean():+.3f} 95% [{np.percentile(T,2.5):+.3f}, {np.percentile(T,97.5):+.3f}]; "
              f"excess {t0 - T.mean():+.3f}; one-sided p={p:.4f}", flush=True)
        out.append({"r1": a.r1, "r2": a.r2, "genes": ng, "stat": t0, "null": null, "null_mean": T.mean(),
                    "lo": np.percentile(T, 2.5), "hi": np.percentile(T, 97.5), "p_one_sided": p, "perms": len(T)})
    with open(a.out + "_retrieval_perm.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)


if __name__ == "__main__":
    main()
