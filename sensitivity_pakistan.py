#!/usr/bin/env python3
"""Sensitivity of the MENA and South Asian test-share differences to the placement of Pakistan and Afghanistan.

The main analysis groups Pakistan and Afghanistan with MENA (fine_regions.py, as in origin_analyses.region_of).
Here every case whose stated origin, or (when the population was assigned from the place of care) whose place of
care, is Pakistan or Afghanistan is moved from MENA to SASIA, and the total-version Mantel-Haenszel difference in the
share of test-based findings (s_i = 0 for cases without a classified finding, as in population_model.py terms) against
EUR is recomputed for MENA and SASIA, with 2,000 permutations between reporting groups within gene x corpus strata.
Both labelings are reported.

  sensitivity_pakistan.py --index results/phenopackets/index.tsv --cases A.tsv [...] --termclass T --groups G --out P
"""
import argparse, csv, json
import numpy as np

from population_model import load_groups, group_of
from symptom_sign import contrast

TEST = {"laboratory", "imaging", "physiological_test", "pathology"}
BED = {"clinical_exam", "symptom_history"}
MOVE = ("pakistan", "afghan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True); ap.add_argument("--cases", action="append", required=True)
    ap.add_argument("--termclass", required=True); ap.add_argument("--groups", required=True)
    ap.add_argument("--perms", type=int, default=2000); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(11)
    G = load_groups(a.groups)
    tc = {x["id"]: x["class"] for x in map(json.loads, open(a.termclass))}
    moved = set()
    for r in csv.DictReader(open(a.index), delimiter="\t"):
        if r["region"] != "MENA": continue
        o = (r.get("origin") or "").lower(); pc = (r.get("place_of_care") or "").lower()
        if any(m in o for m in MOVE) or (not o.strip() and any(m in pc for m in MOVE)): moved.add(r["id"])
    cases = []
    for path in a.cases:
        for r in csv.DictReader(open(path), delimiter="\t"):
            ts = set(r["terms"].split(";")); nt = sum(tc.get(t) in TEST for t in ts); nb = sum(tc.get(t) in BED for t in ts)
            cases.append({"gene": r["gene"], "corpus": r["corpus"], "region": r["region"], "group": group_of(G, r["pmid"]),
                          "alt": "SASIA" if r["id"] in moved and r["region"] == "MENA" else r["region"], "s": nt / max(nt + nb, 1)})
    n_moved = sum(c["alt"] != c["region"] for c in cases)
    print(f"cases moved MENA -> SASIA: {n_moved} (index matches {len(moved)})")
    rows = []
    for lab, name in (("region", "main (Pakistan, Afghanistan in MENA)"), ("alt", "Pakistan, Afghanistan in SASIA")):
        for P in ("MENA", "SASIA"):
            ex, lo, hi, p, n, n1, ns = contrast(cases, lab, P, "EUR", lambda c: (c["gene"], c["corpus"]), {"share": lambda c: c["s"]}, rng, a.perms)
            rows.append({"labeling": name, "population": P, "diff_points": 100 * ex[0], "null_lo": 100 * lo[0], "null_hi": 100 * hi[0],
                         "p": p[0], "n_cases": n, "n_pop": n1, "strata": ns, "cases_moved": n_moved})
            print(f"{name:34s} {P:6s} {100 * ex[0]:+6.2f} [{100 * lo[0]:+.2f}, {100 * hi[0]:+.2f}] p={p[0]:.4f} n={n} ({n1}) strata={ns}")
    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
