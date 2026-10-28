#!/usr/bin/env python3
"""Test-share gap by corpus, and the abstract-versus-full-text bias by population.

(1) Population minus EUR difference in the share of test-based findings, Mantel-Haenszel within gene x corpus strata,
    permutation of labels between reporting groups (symptom_sign.contrast), separately for the Store, PMC-Patients,
    Store + PMC-Patients (full text) and abstracts.
(2) In Store publications, the per-unit difference between the test share of the abstract extraction and of the
    full-text curation (validate_against_store.py units), compared between each population and EUR: difference of
    means, label permutation p (5,000) and bootstrap 95% interval (2,000).

  corpus_replication.py --fine2 results/v4/fine2 --units results/v4/validation/abs_vs_store_units.csv --termclass T
                        --groups G --perms 1000 --out P   -> P_corpus.csv, P_abstract_bias.csv
"""
import argparse, csv, json, os
import numpy as np

from population_model import load_groups, group_of
from symptom_sign import contrast

POPS = ["MENA", "SASIA", "LATAM", "AFR", "EASIA"]
TEST = {"laboratory", "imaging", "physiological_test", "pathology"}
BED = {"clinical_exam", "symptom_history"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fine2", required=True); ap.add_argument("--units", required=True); ap.add_argument("--termclass", required=True)
    ap.add_argument("--groups", required=True); ap.add_argument("--perms", type=int, default=1000); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(5)
    G = load_groups(a.groups)
    tc = {x["id"]: x["class"] for x in map(json.loads, open(a.termclass))}
    F = lambda n: os.path.join(a.fine2, n)
    corpora = [("Store", [F("cases_store.tsv")]), ("PMC-Patients", [F("cases_pmcp_map2.tsv")]),
               ("Store + PMC-Patients", [F("cases_store.tsv"), F("cases_pmcp_map2.tsv")]), ("abstracts", [F("cases_abs.tsv")])]
    rows = []
    for name, files in corpora:
        cases = []
        for f in files:
            for r in csv.DictReader(open(f), delimiter="\t"):
                ts = set(r["terms"].split(";")); nt = sum(tc.get(t) in TEST for t in ts); nb = sum(tc.get(t) in BED for t in ts)
                if nt + nb: cases.append({"gene": r["gene"], "corpus": r["corpus"], "region": r["region"], "group": group_of(G, r["pmid"]),
                                          "s": nt / (nt + nb)})
        for P in POPS:
            res = contrast(cases, "region", P, "EUR", lambda c: (c["gene"], c["corpus"]), {"share": lambda c: c["s"]}, rng, a.perms)
            if res is None: continue
            ex, lo, hi, p, n, n1, ns = res
            rows.append({"corpus": name, "population": P, "diff_points": 100 * ex[0], "null_lo": 100 * lo[0], "null_hi": 100 * hi[0],
                         "p": p[0], "n_cases": n, "n_pop": n1, "strata": ns})
            print(f"{name:22s} {P:6s} {100 * ex[0]:+6.1f}  p={p[0]:.4f}  n={n} ({n1})")
    with open(a.out + "_corpus.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    U = [r for r in csv.DictReader(open(a.units)) if r["level"] == "publication" and r["share_store"] != "nan" and r["share_pred"] != "nan"]
    d = {}
    for r in U: d.setdefault(r["region"], []).append(float(r["share_pred"]) - float(r["share_store"]))
    e = np.array(d["EUR"]); out = [{"population": "EUR", "n": len(e), "bias_points": 100 * e.mean(), "diff_vs_EUR": "", "lo": "", "hi": "", "p": ""}]
    for P in POPS:
        x = np.array(d.get(P, []))
        if len(x) < 3: continue
        obs = x.mean() - e.mean(); allv = np.concatenate([x, e]); null = []
        for _ in range(5000):
            rng.shuffle(allv); null.append(allv[:len(x)].mean() - allv[len(x):].mean())
        p = (1 + np.sum(np.abs(null) >= abs(obs))) / 5001
        bs = [rng.choice(x, len(x)).mean() - rng.choice(e, len(e)).mean() for _ in range(2000)]
        out.append({"population": P, "n": len(x), "bias_points": 100 * x.mean(), "diff_vs_EUR": 100 * obs,
                    "lo": 100 * np.percentile(bs, 2.5), "hi": 100 * np.percentile(bs, 97.5), "p": p})
        print(f"bias {P}: {100 * x.mean():+.1f} vs EUR {100 * e.mean():+.1f}; diff {100 * obs:+.1f} p={p:.3f}")
    with open(a.out + "_abstract_bias.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)


if __name__ == "__main__":
    main()
