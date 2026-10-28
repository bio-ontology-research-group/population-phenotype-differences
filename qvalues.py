#!/usr/bin/env python3
"""Benjamini-Hochberg q-values for the v4 test families; writes results/v4/qvalues.tsv.
Families (BH within each):
  outcomes:<version>   core case-level outcomes (test share, test/bedside counts, finding classes), 5 populations x 9
  organ:<version>      organ-system outcomes (findings, test-based, bedside per system), 5 populations
  strata               core outcomes in age-band and sex subsets
  similarity           all-pairs contrasts (A-A vs A-B; B-B vs A-B), total and patient-carried
  within_vs_outside    P-P vs P-rest contrasts (if present)
  recovery             gene recovery, population vs EUR, total and patient-carried
version in {total, pc, adjage, adjsex}."""
import csv, glob, os, re, collections
import numpy as np


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


rows = []
for f in sorted(glob.glob("results/v4/terms_all/*_outcomes.csv")):
    run = os.path.basename(f)[:-13]
    m = re.match(r"(main|adjage|adjsex|age_\w+?|sex_\w+?)_(MENA|EASIA|SASIA|LATAM|AFR)(?:_(total|pc))?$", run)
    kind, pop, ver = m.group(1), m.group(2), m.group(3)
    version = ver if kind == "main" else kind if kind in ("adjage", "adjsex") else kind
    for r in csv.DictReader(open(f)):
        fam = ("organ:" + version if r["outcome"].startswith("organ[") else
               "outcomes:" + version if kind in ("main", "adjage", "adjsex") else "strata")
        if fam == "strata" and r["outcome"].startswith("organ["): continue
        rows.append({"family": fam, "population": pop, "comparison": f"{pop} vs EUR", "analysis": run, "outcome": r["outcome"],
                     "estimate": float(r["diff"]) - float(r["null_mean"]), "p": float(r["p"]), "n_r1": r["n_r1"], "groups": r["groups"]})
rows = [r for r in rows if not (r["family"] == "strata" and r["outcome"].startswith("organ["))]
for f, fam in (("results/v4/allpairs/contrasts.tsv", "similarity"), ("results/v4/within/contrasts.tsv", "within_vs_outside")):
    if os.path.exists(f):
        for r in csv.DictReader(open(f), delimiter="\t"):
            rows.append({"family": fam, "population": r["pair"], "comparison": r["pair"], "analysis": r["version"],
                         "outcome": r["contrast"], "estimate": float(r["excess"]), "p": float(r["p"]), "n_r1": r["pairs"], "groups": r["units"]})
for f in glob.glob("results/v4/recovery/recovery_*_recovery.csv"):
    for r in csv.DictReader(open(f)):
        rows.append({"family": "recovery", "population": r["r1"], "comparison": f"{r['r1']} vs EUR", "analysis": r["null"],
                     "outcome": "mean log2 rank difference", "estimate": float(r["stat"]) - float(r["null_mean"]),
                     "p": float(r["p_two_sided"]), "n_r1": r["queries"], "groups": r["strata"]})
fams = collections.defaultdict(list)
for i, r in enumerate(rows): fams[r["family"]].append(i)
for fam, idx in fams.items():
    for i, q in zip(idx, bh([rows[i]["p"] for i in idx])): rows[i]["q"] = q; rows[i]["family_size"] = len(idx)
with open("results/v4/qvalues.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["family", "family_size", "comparison", "population", "analysis", "outcome", "estimate", "p", "q", "n_r1", "groups"], delimiter="\t")
    w.writeheader(); w.writerows(sorted(rows, key=lambda r: (r["family"], r["q"])))
for fam, idx in sorted(fams.items()): print(fam, len(idx), "q<0.05:", sum(rows[i]["q"] < 0.05 for i in idx))
