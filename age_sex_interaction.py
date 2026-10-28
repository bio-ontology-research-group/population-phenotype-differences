#!/usr/bin/env python3
"""Does the population difference in the share of test-based findings differ by sex or by age band?

For population P against EUR, the Mantel-Haenszel difference in the per-case test share is computed within
gene x corpus x sex (or x age band) strata, separately for each sex (age band). Population labels are permuted
between reporting groups within these strata (one draw serves all subsets), so the difference between two subset
estimates has a joint permutation distribution. Reported: each subset estimate and the interaction contrasts
male minus female, and each age band minus adults, with permutation p (centred) and BH q over all contrasts.

  age_sex_interaction.py --cases A.tsv [...] --termclass T --groups G --perms 2000 --out P   -> P.csv
"""
import argparse, collections, csv, json
import numpy as np

from population_model import load_groups, group_of, make_plan, draw

POPS = ["MENA", "SASIA", "LATAM", "AFR", "EASIA"]
TEST = {"laboratory", "imaging", "physiological_test", "pathology"}
BED = {"clinical_exam", "symptom_history"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", action="append", required=True); ap.add_argument("--termclass", required=True)
    ap.add_argument("--groups", required=True); ap.add_argument("--perms", type=int, default=2000); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(41)
    G = load_groups(a.groups)
    tc = {x["id"]: x["class"] for x in map(json.loads, open(a.termclass))}
    cases = []
    for path in a.cases:
        for r in csv.DictReader(open(path), delimiter="\t"):
            ts = set(r["terms"].split(";")); nt = sum(tc.get(t) in TEST for t in ts); nb = sum(tc.get(t) in BED for t in ts)
            if nt + nb == 0: continue
            cases.append({"gene": r["gene"], "corpus": r["corpus"], "region": r["region"], "group": group_of(G, r["pmid"]),
                          "sex": r.get("sex", ""), "age": r.get("age_band", ""), "share": nt / (nt + nb)})
    rows = []
    for var, levels, ref in (("sex", ["male", "female"], "female"), ("age", ["infant", "child", "adolescent", "adult"], "adult")):
        for P in POPS:
            cs = [c for c in cases if c["region"] in (P, "EUR") and c[var] in levels]
            key = lambda c: (c["gene"], c["corpus"], c[var])
            sreg = collections.defaultdict(set)
            for c in cs: sreg[key(c)].add(c["region"])
            cs = [c for c in cs if len(sreg[key(c)]) == 2]
            sid = {s: k for k, s in enumerate(sorted({key(c) for c in cs}))}; S = np.array([sid[key(c)] for c in cs]); ns = len(sid)
            lev = np.array([levels.index(c[var]) for c in cs]); y = np.array([c["share"] for c in cs])
            lab0 = np.array([c["region"] for c in cs])
            plan = make_plan(cs, range(len(cs)), lambda i: key(cs[i]), lambda i: cs[i]["group"], lab0)

            def stat(L):
                L1 = (L == P).astype(float); L2 = 1 - L1
                n1 = np.bincount(S, weights=L1, minlength=ns); n2 = np.bincount(S, weights=L2, minlength=ns)
                w = np.where((n1 > 0) & (n2 > 0), n1 * n2 / np.maximum(n1 + n2, 1), 0)
                c = L1 * w[S] / np.maximum(n1[S], 1) - L2 * w[S] / np.maximum(n2[S], 1)
                out = []
                for k in range(len(levels)):
                    m = lev == k; out.append((c[m] @ y[m]) / max(w[np.unique(S[m])].sum(), 1e-12))
                return np.array(out)
            t0 = stat(lab0); null = np.array([stat(draw(plan, lab0, rng)) for _ in range(a.perms)]); mu = null.mean(0)
            r_ = levels.index(ref)
            for k, lv in enumerate(levels):
                p = (1 + np.sum(np.abs(null[:, k] - mu[k]) >= abs(t0[k] - mu[k]))) / (1 + a.perms)
                rows.append({"variable": var, "population": P, "contrast": lv, "estimate_points": 100 * (t0[k] - mu[k]), "p": p,
                             "n_pop": int(((lab0 == P) & (lev == k)).sum()), "kind": "subset"})
                if k == r_: continue
                d0 = (t0[k] - mu[k]) - (t0[r_] - mu[r_]); dn = (null[:, k] - mu[k]) - (null[:, r_] - mu[r_])
                p = (1 + np.sum(np.abs(dn - dn.mean()) >= abs(d0 - dn.mean()))) / (1 + a.perms)
                rows.append({"variable": var, "population": P, "contrast": f"{lv} minus {ref}", "estimate_points": 100 * d0, "p": p,
                             "n_pop": int((lab0 == P).sum()), "kind": "interaction"})
    inter = [r for r in rows if r["kind"] == "interaction"]
    p = np.array([r["p"] for r in inter]); o = np.argsort(p); m = len(p); prev = 1.0
    for k in range(m - 1, -1, -1):
        prev = min(prev, p[o[k]] * m / (k + 1)); inter[o[k]]["q"] = prev
    for r in rows: r.setdefault("q", "")
    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(f"{r['variable']:4s} {r['population']:6s} {r['contrast']:22s} {r['estimate_points']:+6.1f} p={r['p']:.4f} q={r['q'] if r['q'] == '' else round(r['q'], 4)} n={r['n_pop']}")


if __name__ == "__main__":
    main()
