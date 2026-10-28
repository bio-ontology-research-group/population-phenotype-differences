#!/usr/bin/env python3
"""Collect the per-run outputs of the cluster jobs into the tables the paper and qvalues.py use.

  allpairs/contrasts.tsv     from allpairs/popmodel_fine_allpairs.json: for every population pair and version (total,
                             patient-carried) the two contrasts A-A vs A-B and B-B vs A-B, excess = coefficient minus
                             permutation mean, permutation p, pairs, genes, permutation units (parsed from the job
                             log when present) and Benjamini-Hochberg q over the 60 contrasts
  within/contrasts.tsv       the same for within/popmodel_vsrest_allpairs.json (P-P vs P-rest; q is added by qvalues.py)
  terms_all/all_outcomes.tsv concatenation of terms_all/*_outcomes.csv with the run name and excess = diff - null_mean
  setting_all/test_share_q.tsv  the share-of-test-based-findings row of every setting_*_outcomes.csv, in percentage
                             points, with q over the contrasts
  gda/gda_stats_q.csv        gda/gda_stats.csv with q within each analysis version (total, patient-carried)

  collect_tables.py --res results/v4 [--out results/v4]
"""
import argparse, csv, glob, json, os, re

from symptom_sign import bh


def allpairs_table(json_path, log_path, out_path, with_q):
    runs = json.load(open(json_path))
    units = {}
    if log_path and os.path.exists(log_path):
        pair = None
        for line in open(log_path):
            m = re.match(r"#+ (\w+) vs (\w+)", line)
            if m: pair = (m.group(1), m.group(2))
            m = re.match(r"\s*\[(total|patient-carried)\s*\] population_\S+\s.*units (\d+)", line)
            if m and pair: units[(pair, m.group(1))] = m.group(2)   # one permutation plan per pair and version
    rows = []
    for d in runs:
        r1, r2 = d["r1"], d["r2"]
        for ver in ("total", "patient-carried"):
            for name in (f"{r1}{r1}_vs_{r1}{r2}", f"{r2}{r2}_vs_{r1}{r2}"):
                if f"{ver}:{name}_coef" not in d: continue
                ex = d[f"{ver}:{name}_coef"] - d[f"{ver}:{name}_nullmean"]
                if r2 == "REST" and not name.startswith(r1 + r1): continue   # within-vs-rest: P-P vs P-rest only
                rows.append({"pair": f"{r1}-{r2}", "version": ver, "contrast": name, "excess": round(ex, 3),
                             "p": round(d[f"{ver}:{name}_p"], 4), "pairs": d["n_pairs"], "genes": d["n_genes"],
                             "units": units.get(((r1, r2), ver), "")})
    fields = ["pair", "version", "contrast", "excess", "p", "pairs", "genes", "units"]
    if with_q:
        for r, q in zip(rows, bh([float(r["p"]) for r in rows])): r["q_BH"] = round(float(q), 4)
        fields.append("q_BH")
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t"); w.writeheader(); w.writerows(rows)
    print(out_path, len(rows))


def outcomes_table(res, out_path):
    rows = []
    for f in sorted(glob.glob(os.path.join(res, "terms_all", "*_outcomes.csv"))):
        run = os.path.basename(f)[:-len("_outcomes.csv")]
        for r in csv.DictReader(open(f)):
            r = {"run": run, **r}; r["excess"] = float(r["diff"]) - float(r["null_mean"]); rows.append(r)
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    print(out_path, len(rows))


def setting_table(res, out_path):
    rows = []
    for f in sorted(glob.glob(os.path.join(res, "setting_all", "setting_*_outcomes.csv"))):
        name = os.path.basename(f)[len("setting_"):-len("_outcomes.csv")]
        for r in csv.DictReader(open(f)):
            if r["outcome"].startswith("share"):
                rows.append({"contrast": name, "diff_points": f"{100 * (float(r['diff']) - float(r['null_mean'])):+.1f}",
                             "p": f"{float(r['p']):.4f}", "n_cases": r["n_cases"], "n_r1": r["n_r1"]})
    for r, q in zip(rows, bh([float(r["p"]) for r in rows])): r["q_BH"] = f"{q:.4f}"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["contrast", "diff_points", "p", "q_BH", "n_cases", "n_r1"], delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(out_path, len(rows))


def gda_table(res, out_path):
    rows = list(csv.DictReader(open(os.path.join(res, "gda", "gda_stats.csv"))))
    for ver in sorted({r["analysis"] for r in rows}):
        idx = [i for i, r in enumerate(rows) if r["analysis"] == ver and r["p"] not in ("", "nan")]
        for i, q in zip(idx, bh([float(rows[i]["p"]) for i in idx])): rows[i]["q"] = round(float(q), 4)
    for r in rows: r.setdefault("q", "")
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) + (["q"] if "q" not in rows[0] else [])); w.writeheader(); w.writerows(rows)
    print(out_path, len(rows))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", default="results/v4"); ap.add_argument("--out", help="output root (default: --res)")
    a = ap.parse_args(); out = a.out or a.res
    for d in ("allpairs", "within", "terms_all", "setting_all", "gda"): os.makedirs(os.path.join(out, d), exist_ok=True)
    logs = glob.glob(os.path.join(a.res, "allpairs", "allpairs-*.log"))
    allpairs_table(os.path.join(a.res, "allpairs", "popmodel_fine_allpairs.json"), logs[0] if logs else None,
                   os.path.join(out, "allpairs", "contrasts.tsv"), with_q=True)
    allpairs_table(os.path.join(a.res, "within", "popmodel_vsrest_allpairs.json"), None,
                   os.path.join(out, "within", "contrasts.tsv"), with_q=False)
    outcomes_table(a.res, os.path.join(out, "terms_all", "all_outcomes.tsv"))
    setting_table(a.res, os.path.join(out, "setting_all", "test_share_q.tsv"))
    gda_table(a.res, os.path.join(out, "gda", "gda_stats_q.csv"))


if __name__ == "__main__":
    main()
