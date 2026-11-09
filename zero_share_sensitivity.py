#!/usr/bin/env python3
"""Sensitivity of the main test-share comparison to cases with no classified finding.

The main comparison assigns a test share of 0 to cases whose observed classes are all mixed or
unclassified (2.5% to 4.8% of cases per population).  `population_model.py terms --require-classified`
(jobs/zeroshare_um.sbatch) repeats the total-version comparison without these cases.  This script
puts the two side by side: excess (Mantel-Haenszel difference minus permutation mean) in percentage
points, permutation p and Benjamini-Hochberg q over the five populations.

  python3 zero_share_sensitivity.py --main results/v4/terms_all --sens results/v4/terms_noclass \
      --out results/v4/terms_noclass/zero_share_sensitivity
"""
import argparse, csv, os
import numpy as np

POPS = ["MENA", "SASIA", "LATAM", "AFR", "EASIA"]
LABEL = {"MENA": "MENA", "EASIA": "East Asia", "SASIA": "South Asia", "LATAM": "Latin America", "AFR": "Sub-Saharan Africa"}
OUTCOME = "share of findings needing a test"


def read(path):
    for r in csv.DictReader(open(path)):
        if r["outcome"] == OUTCOME:
            return {"excess": 100 * (float(r["diff"]) - float(r["null_mean"])), "p": float(r["p"]), "n": int(r["n_cases"]), "n_r1": int(r["n_r1"])}
    raise SystemExit(f"{OUTCOME} not in {path}")


def bh(ps):
    ps = np.array(ps); m = len(ps); order = np.argsort(ps); q = np.empty(m)
    prev = 1.0
    for rank, i in list(enumerate(order, 1))[::-1]:
        prev = min(prev, m * ps[i] / rank); q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--main", required=True); ap.add_argument("--sens", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    for P in POPS:
        m = read(os.path.join(a.main, f"main_{P}_total_outcomes.csv")); s = read(os.path.join(a.sens, f"noclass_{P}_total_outcomes.csv"))
        rows.append({"population": P, "label": LABEL[P], "n_main": m["n"], "n_r1_main": m["n_r1"], "excess_main": m["excess"], "p_main": m["p"],
                     "n_sens": s["n"], "n_r1_sens": s["n_r1"], "excess_sens": s["excess"], "p_sens": s["p"], "change": s["excess"] - m["excess"]})
    qm, qs = bh([r["p_main"] for r in rows]), bh([r["p_sens"] for r in rows])
    for r, a_, b_ in zip(rows, qm, qs): r["q_main"], r["q_sens"] = a_, b_
    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(a.out + ".tex", "w") as f:
        f.write("\\begin{tabular}{lrrrrrr}\n\\toprule\n& \\multicolumn{3}{c}{All cases ($s_i=0$ if no classified finding)} & \\multicolumn{3}{c}{Cases with a classified finding} \\\\\n")
        f.write("\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\nPopulation & Cases & Difference & $q$ & Cases & Difference & $q$ \\\\\n\\midrule\n")
        for r in rows:
            f.write(f"{r['label']} & {r['n_r1_main']:,} & {r['excess_main']:+.1f} & {r['p_main']:.4f} & {r['n_r1_sens']:,} & {r['excess_sens']:+.1f} & {r['p_sens']:.4f} \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    for r in rows:
        print(f"{r['label']:18s} main {r['excess_main']:+.2f} (n={r['n_r1_main']}, q={r['q_main']:.4f})  excluding {r['excess_sens']:+.2f} (n={r['n_r1_sens']}, q={r['q_sens']:.4f})  change {r['change']:+.2f}")


if __name__ == "__main__":
    main()
