#!/usr/bin/env python3
"""Is the bedside-restriction loss specific to our scorer?  Summary of the Exomiser / LIRICAL re-ranking.

tool_restriction_run.py ranked the causative gene of the paired cases (both bedside and test-based findings) with
Exomiser (phenotype-only; hiPhive default and hiPhive human-only) and LIRICAL (phenotype-only) under three profiles:
all findings, bedside-only findings, and a random subset of the bedside size (export_profiles.py, same draw as
`gda_impact.py rank --random-control`).  As in random_control.py, per tool configuration and population:

  bedside ratio   2^mean(log2 rank_bed  - log2 rank_all)
  random ratio    2^mean(log2 rank_rand - log2 rank_all)
  type effect     2^mean(log2 rank_bed  - log2 rank_rand)

with 95% intervals from 1,000 bootstrap resamples of cases (rng seed 7), top-10 recovery under each profile, and a
pooled row.  Cases whose causative gene is not in the tool's gene universe (--exclude) are left out.  Also the Spearman
correlation, per configuration, of the tool's log2 rank_all with our scorer's (results/v4/gda) and of the per-case
log2(bed/all) of the tool with ours.  Sensitivity (<out>_scored.csv): the pooled ratios restricted to cases whose
causative gene received a non-zero score under all three profiles, so that no rank depends on the midrank convention.

  python3 tool_restriction.py --percase 'results/v5/tool_restriction/percase_*.csv.gz' --gda 'results/v4/gda/gda_chunk*.csv' \
      --exclude results/v5/tool_restriction/not_in_universe.csv --out results/v5/tool_restriction/tool_restriction
"""
import argparse, csv, glob, gzip, collections
import numpy as np
from scipy.stats import spearmanr

POPS = ["EUR", "MENA", "EASIA", "SASIA", "LATAM", "AFR"]
LABEL = {"EUR": "EUR", "MENA": "MENA", "EASIA": "East Asia", "SASIA": "South Asia", "LATAM": "Latin America",
         "AFR": "Sub-Saharan Africa", "ALL": "Pooled"}
CONFIGS = [("exomiser", "hiphive_default", "Exomiser hiPhive (all species, PPI)"), ("exomiser", "hiphive_human", "Exomiser hiPhive (human only)"),
           ("lirical", "lirical_pheno", "LIRICAL (phenotype-only)")]


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def boot(x, rng, B=1000):
    m = x.mean(); idx = rng.integers(0, len(x), (B, len(x)))
    bs = x[idx].mean(axis=1)
    return m, np.percentile(bs, 2.5), np.percentile(bs, 97.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--percase", required=True); ap.add_argument("--gda", required=True)
    ap.add_argument("--exclude"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(7)
    excl = set()
    if a.exclude:
        excl = {(r["tool"], r["id"]) for r in csv.DictReader(open(a.exclude))}
    R = collections.defaultdict(dict); meta = {}; A = collections.Counter()
    for f in sorted(glob.glob(a.percase)):
        for r in csv.DictReader(opener(f)):
            R[(r["tool"], r["config"], r["id"])][r["profile"]] = float(r["rank"])
            meta[r["id"]] = r["region"]; A[(r["tool"], r["config"], r["id"])] += int(r["absent"])
    ours = {}
    for f in sorted(glob.glob(a.gda)):
        for r in csv.DictReader(open(f)):
            if r["rank_bed"] not in ("", "None"): ours[r["id"]] = (float(r["rank_all"]), float(r["rank_bed"]))
    out = []; cors = []; sens = []; rng2 = np.random.default_rng(7)
    for tool, cfg, lab in CONFIGS:
        rows = [(meta[i], v["all"], v["bed"], v["rand"], i) for (t, c, i), v in R.items()
                if t == tool and c == cfg and (t, i) not in excl and all(p in v for p in ("all", "bed", "rand"))]
        if not rows: continue
        incomplete = sum(1 for (t, c, i), v in R.items() if t == tool and c == cfg and (t, i) not in excl and len(v) < 3)
        print(f"{lab}: {len(rows)} cases with all three profiles ({incomplete} incomplete, {sum(1 for (t, i) in excl if t == tool)} excluded, not in universe)")
        for P in POPS + ["ALL"]:
            sub = [x for x in rows if P == "ALL" or x[0] == P]
            if not sub: continue
            ra, rb, rr = (np.array([x[k] for x in sub]) for k in (1, 2, 3))
            la, lb, lr = np.log2(ra), np.log2(rb), np.log2(rr)
            b = boot(lb - la, rng); c = boot(lr - la, rng); d = boot(lb - lr, rng)
            out.append({"tool": tool, "config": cfg, "population": P, "label": LABEL[P], "n": len(sub),
                        "bedside_ratio": 2 ** b[0], "bedside_lo": 2 ** b[1], "bedside_hi": 2 ** b[2],
                        "random_ratio": 2 ** c[0], "random_lo": 2 ** c[1], "random_hi": 2 ** c[2],
                        "type_ratio": 2 ** d[0], "type_lo": 2 ** d[1], "type_hi": 2 ** d[2],
                        "top10_all": 100 * np.mean(ra <= 10), "top10_bed": 100 * np.mean(rb <= 10), "top10_rand": 100 * np.mean(rr <= 10)})
        sc = [x for x in rows if A[(tool, cfg, x[4])] == 0]
        la, lb, lr = (np.log2(np.array([x[k] for x in sc])) for k in (1, 2, 3))
        b = boot(lb - la, rng2); c = boot(lr - la, rng2); d = boot(lb - lr, rng2)
        sens.append({"tool": tool, "config": cfg, "n": len(sc), "n_unscored_any_profile": len(rows) - len(sc),
                     "bedside_ratio": 2 ** b[0], "bedside_lo": 2 ** b[1], "bedside_hi": 2 ** b[2],
                     "random_ratio": 2 ** c[0], "random_lo": 2 ** c[1], "random_hi": 2 ** c[2],
                     "type_ratio": 2 ** d[0], "type_lo": 2 ** d[1], "type_hi": 2 ** d[2]})
        j = [x for x in rows if x[4] in ours]
        ta = np.log2([x[1] for x in j]); tl = np.log2([x[2] for x in j]) - ta
        oa = np.log2([ours[x[4]][0] for x in j]); ol = np.log2([ours[x[4]][1] for x in j]) - oa
        r1 = spearmanr(ta, oa); r2 = spearmanr(tl, ol)
        cors.append({"tool": tool, "config": cfg, "n": len(j), "spearman_log2rank_all": r1.statistic, "p_rank_all": r1.pvalue,
                     "spearman_log2_bed_over_all": r2.statistic, "p_bed_over_all": r2.pvalue})
    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    with open(a.out + "_spearman.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cors[0])); w.writeheader(); w.writerows(cors)
    with open(a.out + "_scored.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(sens[0])); w.writeheader(); w.writerows(sens)
    with open(a.out + ".tex", "w") as f:
        f.write("\\begin{tabular}{llrrrrrrr}\n\\toprule\n")
        f.write("Tool & Population & Cases & Bedside only & Random subset & Bedside vs random & Top ten, full & Top ten, bedside & Top ten, random \\\\\n\\midrule\n")
        for tool, cfg, lab in CONFIGS:
            sub = [o for o in out if o["config"] == cfg]
            for k, o in enumerate(sub):
                f.write(f"{lab if k == 0 else ''} & {o['label']} & {o['n']:,} & {o['bedside_ratio']:.2f} ({o['bedside_lo']:.2f}--{o['bedside_hi']:.2f}) & "
                        f"{o['random_ratio']:.2f} ({o['random_lo']:.2f}--{o['random_hi']:.2f}) & {o['type_ratio']:.2f} ({o['type_lo']:.2f}--{o['type_hi']:.2f}) & "
                        f"{o['top10_all']:.1f} & {o['top10_bed']:.1f} & {o['top10_rand']:.1f} \\\\\n")
            if sub and cfg != CONFIGS[-1][1]: f.write("\\midrule\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    for o in out:
        print(f"{o['config']:16s} {o['label']:18s} n={o['n']:6,} bedside {o['bedside_ratio']:.2f} [{o['bedside_lo']:.2f},{o['bedside_hi']:.2f}]  "
              f"random {o['random_ratio']:.2f} [{o['random_lo']:.2f},{o['random_hi']:.2f}]  bed/rand {o['type_ratio']:.2f} [{o['type_lo']:.2f},{o['type_hi']:.2f}]  "
              f"top10 {o['top10_all']:.1f} / {o['top10_bed']:.1f} / {o['top10_rand']:.1f}")
    for c in cors:
        print(f"{c['config']:16s} n={c['n']:,} Spearman log2 rank_all vs ours {c['spearman_log2rank_all']:.3f}; log2(bed/all) vs ours {c['spearman_log2_bed_over_all']:.3f}")


if __name__ == "__main__":
    main()
