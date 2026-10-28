#!/usr/bin/env python3
"""Ancestry versus place of care, pooled over all non-European ancestries.

Setting labels (setting_labels.py): <A>anc_<A>care (ancestry A cared for in its own region), <A>anc_EURcare (ancestry A
cared for in Europe or North America), EURorlocal_EURcare (reference). Pooled contrasts over A in MENA, SASIA, LATAM,
AFR, EASIA, each with the ancestry in the stratum key so that every comparison is within one ancestry where it involves
two non-reference groups:
  own-region care vs reference        strata gene x corpus x ancestry
  care in Europe/North America vs reference
  own-region care vs care in Europe/North America (same ancestry)
and the same three contrasts per ancestry. Outcomes as in population_model.py terms: share of test-based findings,
log(1+ test-based), log(1+ bedside), log(1+ signs), log(1+ symptoms). Mantel-Haenszel weights, permutation between
reporting groups within strata, BH q within the pooled family and within the per-ancestry family.

  setting_pooled.py --cases fine2/A.tsv [...] --setting setting2/A.tsv [...] --termclass T --groups G --perms 2000 --out P
"""
import argparse, csv, json
import numpy as np

from population_model import load_groups, group_of
from symptom_sign import contrast, bh

ANC = ["MENA", "SASIA", "LATAM", "AFR", "EASIA"]
TEST = {"laboratory", "imaging", "physiological_test", "pathology"}
BED = {"clinical_exam", "symptom_history"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", action="append", required=True); ap.add_argument("--setting", action="append", required=True)
    ap.add_argument("--termclass", required=True); ap.add_argument("--groups", required=True)
    ap.add_argument("--perms", type=int, default=2000); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(43)
    G = load_groups(a.groups)
    tc = {x["id"]: x["class"] for x in map(json.loads, open(a.termclass))}
    setl = {}
    for path in a.setting:
        for r in csv.DictReader(open(path), delimiter="\t"): setl[r["id"]] = r["region"]
    cases = []
    for path in a.cases:
        for r in csv.DictReader(open(path), delimiter="\t"):
            st = setl.get(r["id"], "NA")
            if st == "EURorlocal_EURcare": anc, kind = "EUR", "ref"
            elif "anc_" in st:
                anc, care = st.split("anc_"); care = care.replace("care", "")
                if anc not in ANC: continue
                kind = "eur" if care == "EUR" else ("own" if care == anc else None)
                if kind is None: continue
            else: continue
            ts = set(r["terms"].split(";"))
            nt = sum(tc.get(t) in TEST for t in ts); nb = sum(tc.get(t) in BED for t in ts)
            if nt + nb == 0: continue
            cases.append({"gene": r["gene"], "corpus": r["corpus"], "anc": anc, "kind": kind, "group": group_of(G, r["pmid"]),
                          "nt": nt, "nb": nb, "ns": sum(tc.get(t) == "clinical_exam" for t in ts),
                          "ny": sum(tc.get(t) == "symptom_history" for t in ts)})
    Y = {"share of test-based findings": lambda c: c["nt"] / (c["nt"] + c["nb"]), "log(1+ test-based)": lambda c: np.log1p(c["nt"]),
         "log(1+ bedside)": lambda c: np.log1p(c["nb"]), "log(1+ signs)": lambda c: np.log1p(c["ns"]),
         "log(1+ symptoms)": lambda c: np.log1p(c["ny"])}
    rows = []
    specs = [("own-region care vs reference", "own", "ref"), ("Europe/North America care vs reference", "eur", "ref"),
             ("own-region care vs Europe/North America care", "own", "eur")]
    for fam, ancs in (("pooled", ANC), ("per ancestry", None)):
        for A in (["all"] if ancs else ANC):
            use = [c for c in cases if c["anc"] in (ancs if ancs else [A]) or c["kind"] == "ref"]
            for name, k1, k2 in specs:
                pool = [dict(c, lab=c["kind"]) for c in use if c["kind"] in (k1, k2)]
                # the reference has no ancestry; strata of non-reference groups carry their ancestry, and the reference
                # is copied into each ancestry's strata only through gene x corpus (the reference enters every stratum
                # of its gene x corpus), so key on gene x corpus x ancestry for own vs eur and on gene x corpus otherwise
                if k2 == "ref":
                    key = lambda c: (c["gene"], c["corpus"])
                else:
                    key = lambda c: (c["gene"], c["corpus"], c["anc"])
                res = contrast(pool, "lab", k1, k2, key, Y, rng, a.perms)
                if res is None: continue
                ex, lo, hi, p, n, n1, ns = res
                for j, k in enumerate(Y):
                    rows.append({"family": fam, "ancestry": A, "contrast": name, "outcome": k, "excess": ex[j], "lo": lo[j],
                                 "hi": hi[j], "p": p[j], "n_cases": n, "n_r1": n1, "strata": ns})
    for fam in ("pooled", "per ancestry"):
        sub = [r for r in rows if r["family"] == fam]
        for r, q in zip(sub, bh([r["p"] for r in sub])): r["q"] = q
    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for r in rows:
        sc = 100 if r["outcome"].startswith("share") else 1
        print(f"{r['family']:12s} {r['ancestry']:5s} {r['contrast']:46s} {r['outcome']:30s} {sc * r['excess']:+8.3f} "
              f"[{sc * r['lo']:+.3f}, {sc * r['hi']:+.3f}] p={r['p']:.4f} q={r['q']:.4f} n={r['n_cases']} ({r['n_r1']})")


if __name__ == "__main__":
    main()
