#!/usr/bin/env python3
"""Symptoms versus signs across populations.

Each observed HPO term of a case is a sign (termclass clinical_exam: seen or measured at the bedside) or a symptom
(termclass symptom_history: reported by the patient or family, or taken from the history); test-based and mixed
terms are ignored. Per case: log(1+ signs), log(1+ symptoms) and, for cases with at least one bedside finding, the
share of symptoms among bedside findings. Population vs EUR differences use Mantel-Haenszel weights within strata and
permutation of labels between reporting groups (as population_model.py terms): total (gene x corpus), patient-carried
(+ author region), age-adjusted (+ age band) and sex-adjusted (+ sex). The setting contrasts compare ancestry x
place-of-care groups (setting column) within gene x corpus. BH q within each version (and within the setting family).

  symptom_sign.py --cases A.tsv [...] --setting S.tsv [...] --termclass T --groups G --perms 2000 --out P   -> P.csv
"""
import argparse, collections, csv, json
import numpy as np

from population_model import load_groups, group_of, make_plan, draw

POPS = ["MENA", "SASIA", "LATAM", "AFR", "EASIA"]
SETTING = [("MENAanc_MENAcare", "MENAanc_EURcare"), ("MENAanc_MENAcare", "EURorlocal_EURcare"),
           ("MENAanc_EURcare", "EURorlocal_EURcare"), ("SASIAanc_SASIAcare", "EURorlocal_EURcare"),
           ("AFRanc_AFRcare", "EURorlocal_EURcare"), ("LATAManc_LATAMcare", "EURorlocal_EURcare"),
           ("EASIAanc_EASIAcare", "EURorlocal_EURcare")]


def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); q = np.empty(m); prev = 1.0
    for k in range(m - 1, -1, -1):
        prev = min(prev, p[o[k]] * m / (k + 1)); q[o[k]] = prev
    return q


def contrast(cs, labcol, r1, r2, key, Y, rng, perms):
    cs = [c for c in cs if c[labcol] in (r1, r2)]
    sreg = collections.defaultdict(set)
    for c in cs: sreg[key(c)].add(c[labcol])
    idx = [i for i, c in enumerate(cs) if len(sreg[key(c)]) == 2]
    cs = [cs[i] for i in idx]
    if not cs: return None
    sid = {s: k for k, s in enumerate(sorted({key(c) for c in cs}))}; S = np.array([sid[key(c)] for c in cs]); ns = len(sid)
    Yc = np.array([[f(c) for f in Y.values()] for c in cs], float); ok = ~np.isnan(Yc)
    lab0 = np.array([c[labcol] for c in cs])
    plan = make_plan(cs, range(len(cs)), lambda i: key(cs[i]), lambda i: cs[i]["group"], lab0)

    def stat(L):
        res = []
        for j in range(Yc.shape[1]):
            m = ok[:, j]; L1 = ((L == r1) & m).astype(float); L2 = ((L == r2) & m).astype(float)
            n1 = np.bincount(S, weights=L1, minlength=ns); n2 = np.bincount(S, weights=L2, minlength=ns)
            w = np.where((n1 > 0) & (n2 > 0), n1 * n2 / np.maximum(n1 + n2, 1), 0)
            c = (L1 * w[S] / np.maximum(n1[S], 1) - L2 * w[S] / np.maximum(n2[S], 1)) / max(w.sum(), 1e-12)
            res.append(c @ np.where(m, Yc[:, j], 0))
        return np.array(res)
    t0 = stat(lab0); null = np.array([stat(draw(plan, lab0, rng)) for _ in range(perms)]); mu = null.mean(0)
    p = (1 + (np.abs(null - mu) >= np.abs(t0 - mu)).sum(0)) / (1 + perms)
    lo, hi = np.percentile(null - mu, [2.5, 97.5], axis=0)
    return t0 - mu, lo, hi, p, len(cs), int((lab0 == r1).sum()), ns


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", action="append", required=True); ap.add_argument("--termclass", required=True)
    ap.add_argument("--setting", action="append", default=[], help="setting2 tables: their region column holds the setting label")
    ap.add_argument("--groups", required=True); ap.add_argument("--perms", type=int, default=2000); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(37)
    G = load_groups(a.groups)
    tc = {x["id"]: x["class"] for x in map(json.loads, open(a.termclass))}
    setl = {}
    for path in a.setting:
        for r in csv.DictReader(open(path), delimiter="\t"): setl[r["id"]] = r["region"]
    cases = []
    for path in a.cases:
        for r in csv.DictReader(open(path), delimiter="\t"):
            ts = set(r["terms"].split(";"))
            ns_ = sum(tc.get(t) == "clinical_exam" for t in ts); ny = sum(tc.get(t) == "symptom_history" for t in ts)
            cases.append({"gene": r["gene"], "corpus": r["corpus"], "region": r["region"], "setting": setl.get(r["id"], "NA"),
                          "author": r.get("author_region") or "UNKNOWN", "age": r.get("age_band", ""), "sex": r.get("sex", ""),
                          "group": group_of(G, r["pmid"]), "signs": ns_, "symptoms": ny})
    Y = {"log(1+ signs)": lambda c: np.log1p(c["signs"]), "log(1+ symptoms)": lambda c: np.log1p(c["symptoms"]),
         "share of symptoms among bedside findings": lambda c: c["symptoms"] / (c["signs"] + c["symptoms"]) if c["signs"] + c["symptoms"] else np.nan}
    for r in ("EUR",) + tuple(POPS):
        sub = [c for c in cases if c["region"] == r and c["signs"] + c["symptoms"] > 0]
        print(f"{r:6s} cases {len(sub):6d}  mean signs {np.mean([c['signs'] for c in sub]):.2f}  symptoms {np.mean([c['symptoms'] for c in sub]):.2f}  "
              f"symptom share {np.mean([c['symptoms'] / (c['signs'] + c['symptoms']) for c in sub]):.3f}")
    versions = {"total": lambda c: (c["gene"], c["corpus"]), "patient-carried": lambda c: (c["gene"], c["corpus"], c["author"]),
                "age-adjusted": lambda c: (c["gene"], c["corpus"], c["age"]), "sex-adjusted": lambda c: (c["gene"], c["corpus"], c["sex"])}
    rows = []
    for v, key in versions.items():
        pool = cases if v in ("total", "patient-carried") else [c for c in cases if c["age" if v == "age-adjusted" else "sex"] in
                                                                  (("infant", "child", "adolescent", "adult") if v == "age-adjusted" else ("male", "female"))]
        for P in POPS:
            res = contrast(pool, "region", P, "EUR", key, Y, rng, a.perms)
            if res is None: continue
            ex, lo, hi, p, n, n1, ns = res
            for j, k in enumerate(Y):
                rows.append({"family": v, "contrast": f"{P} vs EUR", "outcome": k, "excess": ex[j], "lo": lo[j], "hi": hi[j], "p": p[j],
                             "n_cases": n, "n_r1": n1, "strata": ns})
    for r1, r2 in SETTING:
        res = contrast(cases, "setting", r1, r2, versions["total"], Y, rng, a.perms)
        if res is None: continue
        ex, lo, hi, p, n, n1, ns = res
        for j, k in enumerate(Y):
            rows.append({"family": "setting", "contrast": f"{r1} vs {r2}", "outcome": k, "excess": ex[j], "lo": lo[j], "hi": hi[j],
                         "p": p[j], "n_cases": n, "n_r1": n1, "strata": ns})
    for fam in {r["family"] for r in rows}:
        sub = [r for r in rows if r["family"] == fam]
        for r, q in zip(sub, bh([r["p"] for r in sub])): r["q"] = q
    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(f"{r['family']:15s} {r['contrast']:42s} {r['outcome']:42s} {r['excess']:+.4f}  p={r['p']:.4f} q={r['q']:.4f}  n={r['n_cases']} ({r['n_r1']})")


if __name__ == "__main__":
    main()
