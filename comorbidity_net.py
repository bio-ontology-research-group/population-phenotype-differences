#!/usr/bin/env python3
"""Co-reporting networks of organ systems by population.

Nodes are organ systems (children of HP:0000118 with findings in at least --min-cases cases), optionally split into
their test-based and bedside findings (--split). A case "has" a node when one of its HPO terms lies below the organ
system (and, with --split, has the matching finding class in termclass.jsonl).

Edge strength is the within-gene co-reporting correlation: node indicators are centred within each gene x corpus
stratum separately for each population, so gene composition cannot create an edge, and the correlation is computed
from the pooled centred indicators. For population P against EUR, the edge statistic is the difference of the
Fisher-transformed correlations, atanh(r_P) - atanh(r_EUR), computed on strata that contain both populations.
Population labels are permuted between reporting groups within strata (as in population_model.py terms); the
reported excess is observed minus permutation mean, p is two-sided around the permutation mean, and q is
Benjamini-Hochberg over the edges of one comparison.

  comorbidity_net.py --cases A.tsv [...] --obo hp.obo --termdb ontofunc_hp --groups groups.tsv --termclass T
                     --pops MENA,SASIA,LATAM,AFR,EASIA [--split] --perms 1000 --out P
writes P_edges.csv (one row per population and edge), P_global.csv (global test: sum of squared standardised
excesses against its permutation distribution) and P_nodes.csv.
"""
import argparse, collections, csv, json, os
import numpy as np

from population_model import load_groups, group_of, make_plan, draw

TEST = {"laboratory", "imaging", "physiological_test", "pathology"}
BED = {"clinical_exam", "symptom_history"}


def centred(X, gid, ng):
    """Subtract the mean of each column within groups gid; return the centred matrix and the group sizes."""
    n = np.bincount(gid, minlength=ng).astype(float)
    M = np.stack([np.bincount(gid, weights=X[:, j], minlength=ng) for j in range(X.shape[1])], 1) / np.maximum(n, 1)[:, None]
    return X - M[gid], n


def corr_by_label(X, S, ns, L):
    """Within-stratum, within-label centred correlation matrices for labels 1 (r1) and 0 (r2)."""
    out = []
    for v in (1, 0):
        m = L == v
        Xc, n = centred(X[m], S[m], ns)
        C = Xc.T @ Xc
        d = np.sqrt(np.maximum(np.diag(C), 1e-12))
        out.append(C / np.outer(d, d))
    return out


def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); q = np.empty(m); prev = 1.0
    for k in range(m - 1, -1, -1):
        prev = min(prev, p[o[k]] * m / (k + 1)); q[o[k]] = prev
    return q


def main():
    from rank_engine import parse_obo
    from phenotype_specificity import names, ancestors_fn, ROOT_PA
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", action="append", required=True); ap.add_argument("--obo", required=True)
    ap.add_argument("--termdb"); ap.add_argument("--groups", required=True); ap.add_argument("--termclass", required=True)
    ap.add_argument("--pops", default="MENA,SASIA,LATAM,AFR,EASIA"); ap.add_argument("--ref", default="EUR")
    ap.add_argument("--split", action="store_true"); ap.add_argument("--min-cases", type=int, default=200)
    ap.add_argument("--perms", type=int, default=1000); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(29)
    G = load_groups(a.groups)
    terms, alt = parse_obo(a.obo); nm = names(a.obo)
    anc = ancestors_fn(terms, os.path.join(a.termdb, "graph_path.txt") if a.termdb else None)
    tc = {x["id"]: x["class"] for x in map(json.loads, open(a.termclass))}
    BR = sorted(t for t in terms if ROOT_PA in terms[t]["is_a"] and not terms[t]["obsolete"])
    pops = a.pops.split(",")
    cases = []
    for path in a.cases:
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["region"] not in pops + [a.ref]: continue
            ts = [alt.get(t, t) for t in r["terms"].split(";")]
            ts = [t for t in ts if t in terms and not terms[t]["obsolete"]]
            if not ts: continue
            feat = set()
            for t in set(ts):
                for b in BR:
                    if b in anc(t):
                        if a.split:
                            c = tc.get(t)
                            if c in TEST: feat.add((b, "test"))
                            elif c in BED: feat.add((b, "bedside"))
                        else:
                            feat.add((b, ""))
            cases.append({"gene": r["gene"], "corpus": r.get("corpus", ""), "region": r["region"],
                          "group": group_of(G, r["pmid"]), "feat": feat})
    cnt = collections.Counter(f for c in cases for f in c["feat"])
    nodes = sorted(f for f, n in cnt.items() if n >= a.min_cases)
    col = {f: j for j, f in enumerate(nodes)}
    label = {f: nm.get(f[0], f[0]).replace("Abnormality of the ", "").replace("Abnormality of ", "") + (f" ({f[1]})" if f[1] else "")
             for f in nodes}
    print(f"{len(cases)} cases, {len(nodes)} nodes")
    with open(a.out + "_nodes.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["node", "hpo", "kind", "label", "n_cases"])
        for nd in nodes: w.writerow([col[nd], nd[0], nd[1], label[nd], cnt[nd]])
    iu = np.triu_indices(len(nodes), 1)
    rows, glob = [], []
    # reference network: EUR, within-gene centred correlation over all EUR cases
    ref = [c for c in cases if c["region"] == a.ref]
    sid = {s: k for k, s in enumerate(sorted({(c["gene"], c["corpus"]) for c in ref}))}
    X = np.zeros((len(ref), len(nodes)), np.float32)
    for i, c in enumerate(ref):
        for f in c["feat"]:
            if f in col: X[i, col[f]] = 1
    Xc, _ = centred(X, np.array([sid[(c["gene"], c["corpus"])] for c in ref]), len(sid))
    C = Xc.T @ Xc; d = np.sqrt(np.maximum(np.diag(C), 1e-12)); Rref = C / np.outer(d, d)
    for P in pops:
        cs = [c for c in cases if c["region"] in (P, a.ref)]
        sreg = collections.defaultdict(set)
        for c in cs: sreg[(c["gene"], c["corpus"])].add(c["region"])
        cs = [c for c in cs if len(sreg[(c["gene"], c["corpus"])]) == 2]
        sid = {s: k for k, s in enumerate(sorted({(c["gene"], c["corpus"]) for c in cs}))}
        S = np.array([sid[(c["gene"], c["corpus"])] for c in cs]); ns = len(sid)
        X = np.zeros((len(cs), len(nodes)), np.float32)
        for i, c in enumerate(cs):
            for f in c["feat"]:
                if f in col: X[i, col[f]] = 1
        lab0 = np.array([c["region"] for c in cs])
        plan = make_plan(cs, range(len(cs)), lambda i: (cs[i]["gene"], cs[i]["corpus"]), lambda i: cs[i]["group"], lab0)

        def stat(L):
            R1, R2 = corr_by_label(X, S, ns, (L == P).astype(int))
            return (np.arctanh(np.clip(R1, -0.999, 0.999)) - np.arctanh(np.clip(R2, -0.999, 0.999)))[iu], R1[iu], R2[iu]
        t0, r1, r2 = stat(lab0)
        null = np.empty((a.perms, len(t0)))
        for b in range(a.perms): null[b] = stat(draw(plan, lab0, rng))[0]
        mu = null.mean(0)
        p = (1 + (np.abs(null - mu) >= np.abs(t0 - mu)).sum(0)) / (1 + a.perms)
        q = bh(p)
        sd = null.std(0) + 1e-12
        g0 = (((t0 - mu) / sd) ** 2).sum(); gn = (((null - mu) / sd) ** 2).sum(1)
        pg = (1 + (gn >= g0).sum()) / (1 + a.perms)
        glob.append({"pop": P, "global_stat": g0, "null_mean": gn.mean(), "p_global": pg, "edges_q05": int((q < 0.05).sum()),
                     "n_cases": len(cs), "n_pop": int((lab0 == P).sum()), "strata": ns})
        print(f"{P}: global p={pg:.4f}; {len(cs)} cases ({(lab0 == P).sum()} {P}), {ns} strata; edges q<0.05: {(q < 0.05).sum()} "
              f"(stronger {((q < 0.05) & (t0 > mu)).sum()}, weaker {((q < 0.05) & (t0 < mu)).sum()})")
        for e, (i, j) in enumerate(zip(*iu)):
            rows.append({"pop": P, "i": i, "j": j, "node_i": label[nodes[i]], "node_j": label[nodes[j]],
                         "r_ref_all": Rref[i, j], "r_pop": r1[e], "r_ref": r2[e], "diff": t0[e], "null_mean": mu[e],
                         "excess": t0[e] - mu[e], "p": p[e], "q": q[e], "n_cases": len(cs), "n_pop": int((lab0 == P).sum())})
    with open(a.out + "_edges.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(a.out + "_global.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(glob[0])); w.writeheader(); w.writerows(glob)
    print("written", a.out + "_edges.csv")


if __name__ == "__main__":
    main()
