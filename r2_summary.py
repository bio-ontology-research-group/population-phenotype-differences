#!/usr/bin/env python3
"""Summarise the variance decomposition of within-gene case-to-case similarity (population_model.py pairs --r2-only).

For every population pair: the between-gene share of the variance of Z, the within-gene R2 of the full pair model
(patient-carried version: depth, curation process, reporting group, author region, population), and for each layer
the R2 when it is the only term (alone) and the drop in R2 when it is removed from the full model (unique), with
bootstrap intervals over reporting groups. Writes results/v4/r2/summary.csv (per pair, unique shares of the explained
within-gene variance) and prints medians and ranges over the pairs.

  r2_summary.py [--res results/v4]
"""
import argparse, csv, glob, os, statistics


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--res", default="results/v4"); a = ap.parse_args()
    layers = ["phenotyping depth", "curation process", "reporting group", "author region", "population"]
    rows = []
    for f in sorted(glob.glob(os.path.join(a.res, "r2", "popmodel_r2_*_r2.csv"))):
        d = {r["quantity"]: r for r in csv.DictReader(open(f))}
        if not d: continue
        pair = next(iter(d.values()))["pair"]; full = float(d["within-gene R2, full model"]["value"])
        r = {"pair": pair, "n_pairs": next(iter(d.values()))["n_pairs"], "n_genes": next(iter(d.values()))["n_genes"],
             "between_gene_share": float(d["between-gene share of variance"]["value"]), "within_gene_R2_full": full}
        for L in layers:
            r[L + " unique"] = float(d[L + ": unique"]["value"]); r[L + " alone"] = float(d[L + ": alone"]["value"])
            r[L + " unique lo"] = float(d[L + ": unique"]["lo"]); r[L + " unique hi"] = float(d[L + ": unique"]["hi"])
        pheno = sum(r[L + " unique"] for L in layers[:4])
        r["phenotyping unique share of R2"] = pheno / full; r["population unique share of R2"] = r["population unique"] / full
        rows.append(r)
    with open(os.path.join(a.res, "r2", "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    def med(k): v = [r[k] for r in rows]; return f"median {statistics.median(v):.4f} range {min(v):.4f} to {max(v):.4f}"
    print(f"{len(rows)} pairs")
    for k in ["between_gene_share", "within_gene_R2_full"] + [L + " unique" for L in layers] + [L + " alone" for L in layers] + \
             ["phenotyping unique share of R2", "population unique share of R2"]:
        print(f"{k:38s} {med(k)}")
    for r in rows:
        print(f"{r['pair']:12s} R2 {r['within_gene_R2_full']:.3f}  depth {r['phenotyping depth unique']:.4f}  curation {r['curation process unique']:.4f}  "
              f"group {r['reporting group unique']:.4f}  author {r['author region unique']:.4f}  population {r['population unique']:.4f} "
              f"[{r['population unique lo']:.4f}, {r['population unique hi']:.4f}]  pheno share {r['phenotyping unique share of R2']:.2f} pop share {r['population unique share of R2']:.2f}")


if __name__ == "__main__":
    main()
