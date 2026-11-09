#!/usr/bin/env python3
"""Control for the bedside restriction: remove the same number of findings at random.

For every paired case (both bedside and test-based findings), `gda_impact.py rank --random-control`
ranked the causative gene with a random subset of the case's findings of the same size as its
bedside profile (seeded per case).  This script joins those ranks with the main ranking output
(rank_all, rank_bed) and reports, per population,

  bedside ratio   2^mean(log2 rank_bed  - log2 rank_all)   (as in Figure 3a)
  random ratio    2^mean(log2 rank_rand - log2 rank_all)   (same number of findings removed, type unmatched)
  type effect     2^mean(log2 rank_bed  - log2 rank_rand)  (bedside restriction beyond random removal)

with 95% intervals from 1,000 bootstrap resamples of cases, and top-10 recovery under each profile.

  python3 random_control.py --gda 'results/v4/gda/gda_chunk*.csv' --rand 'results/v4/gda_rand/gdarand_chunk*.csv' \
      --out results/v4/gda_rand/random_control
"""
import argparse, csv, glob
import numpy as np

POPS = ["EUR", "MENA", "EASIA", "SASIA", "LATAM", "AFR"]
LABEL = {"EUR": "EUR", "MENA": "MENA", "EASIA": "East Asia", "SASIA": "South Asia", "LATAM": "Latin America", "AFR": "Sub-Saharan Africa"}


def load(pattern):
    return {r["id"]: r for f in sorted(glob.glob(pattern)) for r in csv.DictReader(open(f))}


def boot(x, rng, B=1000):
    m = x.mean(); idx = rng.integers(0, len(x), (B, len(x)))
    bs = x[idx].mean(axis=1)
    return m, np.percentile(bs, 2.5), np.percentile(bs, 97.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gda", required=True); ap.add_argument("--rand", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(7)
    gda, rnd = load(a.gda), load(a.rand)
    rows = []; mism = 0
    for i, r in rnd.items():
        g = gda.get(i)
        if not g or not g["rank_bed"] or g["rank_bed"] == "None": continue
        mism += g["rank_all"] != r["rank_all"]  # floating-point tie order; the main output's rank_all is used throughout
        rows.append((g["region"], float(g["rank_all"]), float(g["rank_bed"]), float(r["rank_rand"]), int(r["n_rand_test"]), int(r["n_test"])))
    print(f"{len(rows)} paired cases joined ({len(rnd)} random-control rows, {sum(1 for g in gda.values() if g['rank_bed'] not in ('', 'None') and int(g['n_test']) > 0)} paired cases in the main output); rank_all differs between the two runs in {mism} cases")
    out = []
    for P in POPS:
        sub = [x for x in rows if x[0] == P]
        if not sub: continue
        ra, rb, rr = (np.array([x[k] for x in sub]) for k in (1, 2, 3))
        la, lb, lr = np.log2(ra), np.log2(rb), np.log2(rr)
        b = boot(lb - la, rng); c = boot(lr - la, rng); d = boot(lb - lr, rng)
        frac_test_kept = np.mean([x[4] / x[5] for x in sub])
        out.append({"population": P, "label": LABEL[P], "n": len(sub),
                    "bedside_ratio": 2 ** b[0], "bedside_lo": 2 ** b[1], "bedside_hi": 2 ** b[2],
                    "random_ratio": 2 ** c[0], "random_lo": 2 ** c[1], "random_hi": 2 ** c[2],
                    "type_ratio": 2 ** d[0], "type_lo": 2 ** d[1], "type_hi": 2 ** d[2],
                    "top10_all": 100 * np.mean(ra <= 10), "top10_bed": 100 * np.mean(rb <= 10), "top10_rand": 100 * np.mean(rr <= 10),
                    "mean_share_of_test_findings_kept_by_random": frac_test_kept})
    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    with open(a.out + ".tex", "w") as f:
        f.write("\\begin{tabular}{lrrrrrrr}\n\\toprule\n")
        f.write("Population & Cases & Bedside only & Random subset & Bedside vs random & Top ten, full & Top ten, bedside & Top ten, random \\\\\n\\midrule\n")
        for o in out:
            f.write(f"{o['label']} & {o['n']:,} & {o['bedside_ratio']:.2f} ({o['bedside_lo']:.2f}--{o['bedside_hi']:.2f}) & "
                    f"{o['random_ratio']:.2f} ({o['random_lo']:.2f}--{o['random_hi']:.2f}) & {o['type_ratio']:.2f} ({o['type_lo']:.2f}--{o['type_hi']:.2f}) & "
                    f"{o['top10_all']:.1f} & {o['top10_bed']:.1f} & {o['top10_rand']:.1f} \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    for o in out:
        print(f"{o['label']:18s} n={o['n']:6,} bedside {o['bedside_ratio']:.2f} [{o['bedside_lo']:.2f},{o['bedside_hi']:.2f}]  "
              f"random {o['random_ratio']:.2f} [{o['random_lo']:.2f},{o['random_hi']:.2f}]  bed/rand {o['type_ratio']:.2f} [{o['type_lo']:.2f},{o['type_hi']:.2f}]  "
              f"top10 {o['top10_all']:.1f} / {o['top10_bed']:.1f} / {o['top10_rand']:.1f}  test kept by random {100*o['mean_share_of_test_findings_kept_by_random']:.0f}%")


if __name__ == "__main__":
    main()
