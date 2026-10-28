#!/usr/bin/env python3
"""Follow-up analyses on the per-case outputs of run_analyses.py.

  ci      : bootstrap 95% CIs for Hits@1 / MRR per cohort (baseline_percase.csv, extrinsic)
  matched : circularity-free, depth-matched comparison (depth_percase.csv x lpo_percase.csv)
  source  : source-distinguishability (identifiability) test: can a classifier tell two
            sources apart from the HPO profile alone? AUC under 5-fold CV, all cases and
            shared-gene cases, plus a random-split placebo within one source.
"""
import argparse, csv, os, json, random, collections
import numpy as np
from rank_engine import Engine

random.seed(1); np.random.seed(1)


def boot(vals, fn, B=2000):
    v = np.array(vals, dtype=float); n = len(v)
    stats = [fn(v[np.random.randint(0, n, n)]) for _ in range(B)]
    return fn(v), np.percentile(stats, 2.5), np.percentile(stats, 97.5)


def cmd_ci(a):
    rows = [r for r in csv.DictReader(open(os.path.join(a.out, "baseline_percase.csv"))) if r["ic_type"] == "extrinsic"]
    groups = collections.defaultdict(list)
    for r in rows:
        groups[r["cohort"]].append(int(r["rank"]))
        groups[r["cohort"] + "/" + r["source"]].append(int(r["rank"]))
    out = []
    for g, ranks in sorted(groups.items()):
        h, lo, hi = boot(ranks, lambda v: 100 * np.mean(v == 1))
        m, mlo, mhi = boot(ranks, lambda v: np.mean(1 / v))
        h10, lo10, hi10 = boot(ranks, lambda v: 100 * np.mean(v <= 10))
        out.append({"group": g, "n": len(ranks), "hits1": round(h, 2), "hits1_lo": round(lo, 2), "hits1_hi": round(hi, 2),
                    "hits10": round(h10, 2), "hits10_lo": round(lo10, 2), "hits10_hi": round(hi10, 2),
                    "mrr": round(m, 4), "mrr_lo": round(mlo, 4), "mrr_hi": round(mhi, 4)})
        print(f"  {g:28s} n={len(ranks):5d} Hits@1 {h:5.2f} [{lo:.2f}, {hi:.2f}]  Hits@10 {h10:5.2f} [{lo10:.2f}, {hi10:.2f}]  MRR {m:.4f} [{mlo:.4f}, {mhi:.4f}]")
    with open(os.path.join(a.out, "baseline_ci.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)


def cmd_matched(a):
    lpo = {r["id"]: r for r in csv.DictReader(open(os.path.join(a.out, "lpo_percase.csv")))}
    depth = list(csv.DictReader(open(os.path.join(a.out, "depth_percase.csv"))))
    base = [r for r in csv.DictReader(open(os.path.join(a.out, "baseline_percase.csv"))) if r["ic_type"] == "extrinsic"]
    def summ(name, ranks):
        h, lo, hi = boot(ranks, lambda v: 100 * np.mean(v == 1)); m, mlo, mhi = boot(ranks, lambda v: np.mean(1 / v))
        h10 = 100 * np.mean(np.array(ranks) <= 10)
        print(f"  {name:62s} n={len(ranks):5d} Hits@1 {h:5.2f} [{lo:.2f}, {hi:.2f}] Hits@10 {h10:5.2f} MRR {m:.4f} [{mlo:.4f}, {mhi:.4f}]")
        return {"subset": name, "n": len(ranks), "hits1": round(h, 2), "hits1_lo": round(lo, 2), "hits1_hi": round(hi, 2), "hits10": round(h10, 2), "mrr": round(m, 4), "mrr_lo": round(mlo, 4), "mrr_hi": round(mhi, 4)}
    res = []
    nonc = [r for r in depth if lpo[r["id"]]["pmid_in_hpoa"] == "0"]
    circ = [r for r in depth if lpo[r["id"]]["pmid_in_hpoa"] == "1"]
    res.append(summ("Store, PMID not in HPOA, full depth", [int(r["rank_full"]) for r in nonc if r["rep"] == "0"]))
    res.append(summ("Store, PMID not in HPOA, random subset at clinical depth", [int(r["rank_random_subset"]) for r in nonc]))
    res.append(summ("Store, PMID not in HPOA, top-IC subset at clinical depth", [int(r["rank_topic_subset"]) for r in nonc]))
    res.append(summ("Store, PMID in HPOA, full depth", [int(r["rank_full"]) for r in circ if r["rep"] == "0"]))
    res.append(summ("Store, PMID in HPOA, random subset at clinical depth", [int(r["rank_random_subset"]) for r in circ]))
    res.append(summ("Store, PMID in HPOA, LPO reference, full depth", [int(lpo[r["id"]]["rank_lpo"]) for r in circ if r["rep"] == "0"]))
    for coh in ("Saudi-literature", "Saudi-clinical", "US-clinical"):
        res.append(summ(f"{coh} (as recorded)", [int(r["rank"]) for r in base if r["cohort"] == coh]))
    # mean IC of the subsets vs clinical cohorts
    print("  mean IC: Store random subset %.3f, top-IC subset %.3f, Saudi-literature %.3f, Saudi-clinical %.3f" % (
        np.mean([float(r["mean_ic_random"]) for r in depth]), np.mean([float(r["mean_ic_top"]) for r in depth]),
        np.mean([float(r["mean_ic"]) for r in base if r["cohort"] == "Saudi-literature"]),
        np.mean([float(r["mean_ic"]) for r in base if r["cohort"] == "Saudi-clinical"])))
    with open(os.path.join(a.out, "matched_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0].keys())); w.writeheader(); w.writerows(res)


def cmd_source(a):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score
    import scipy.sparse as sp
    from run_analyses import load_pavs, load_store, evaluable
    eng = Engine(os.path.join(a.hpo_dir, "hp.obo"), os.path.join(a.hpo_dir, "genes_to_phenotype.txt"))
    cases = evaluable(load_pavs(a.pavs) + load_store(a.store), eng)
    def label(c):
        return c["source"] if c["source"] != "Store" else "Store"
    def X_of(cs):
        rows, cols = [], []
        for i, c in enumerate(cs):
            anc = set()
            for t in c["terms"]:
                anc |= eng.anc_sets[t]
            for t in anc:
                rows.append(i); cols.append(eng.tid[t])
        return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(cs), len(eng.tid)))
    def cv_auc(cs, y, seed=0):
        X = X_of(cs); y = np.array(y)
        skf = StratifiedKFold(5, shuffle=True, random_state=seed)
        pred = np.zeros(len(y))
        for tr, te in skf.split(X, y):
            clf = LogisticRegression(C=0.3, max_iter=2000, class_weight="balanced")
            clf.fit(X[tr], y[tr]); pred[te] = clf.predict_proba(X[te])[:, 1]
        return roc_auc_score(y, pred)
    groups = {"A+B": [c for c in cases if c["source"] in "AB"], "F": [c for c in cases if c["source"] == "F"],
              "P": [c for c in cases if c["source"] == "P"], "M": [c for c in cases if c["source"] == "M"],
              "Q": [c for c in cases if c["source"] == "Q"], "Store": [c for c in cases if c["source"] == "Store"]}
    pairs = [("F", "P"), ("A+B", "P"), ("A+B", "F"), ("P", "Q"), ("F", "Q"), ("A+B", "Q"), ("M", "P"), ("M", "Store"), ("P", "Store"), ("Q", "Store")]
    out = []
    for g1, g2 in pairs:
        c1, c2 = groups[g1], groups[g2]
        n = min(len(c1), len(c2), 1500)
        s1 = random.sample(c1, n); s2 = random.sample(c2, n)
        auc_all = cv_auc(s1 + s2, [0] * n + [1] * n)
        # shared-gene subset: only cases whose gene occurs in both groups (controls disease mix)
        shared = {c["gene"] for c in c1} & {c["gene"] for c in c2}
        sc1 = [c for c in c1 if c["gene"] in shared]; sc2 = [c for c in c2 if c["gene"] in shared]
        m = min(len(sc1), len(sc2), 1500)
        auc_shared = cv_auc(random.sample(sc1, m) + random.sample(sc2, m), [0] * m + [1] * m) if m >= 40 else float("nan")
        # depth-only model: n_terms and mean IC
        print(f"  {g1:5s} vs {g2:5s} n={n:4d} AUC(all)={auc_all:.3f}  shared-gene n={m:4d} AUC={auc_shared:.3f}", flush=True)
        out.append({"g1": g1, "g2": g2, "n_per_group": n, "auc_all": round(auc_all, 3), "n_shared": m, "auc_shared_gene": round(auc_shared, 3) if m >= 40 else ""})
    # placebo: random halves of P, and of Store
    for g in ("P", "Store"):
        cs = random.sample(groups[g], min(len(groups[g]), 1500)); y = [random.randint(0, 1) for _ in cs]
        auc = cv_auc(cs, y)
        print(f"  placebo random split within {g}: AUC={auc:.3f}")
        out.append({"g1": g, "g2": g + " (random split)", "n_per_group": len(cs) // 2, "auc_all": round(auc, 3), "n_shared": "", "auc_shared_gene": ""})
    with open(os.path.join(a.out, "source_distinguishability.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["ci", "matched", "source"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--hpo-dir"); ap.add_argument("--pavs"); ap.add_argument("--store")
    a = ap.parse_args()
    {"ci": cmd_ci, "matched": cmd_matched, "source": cmd_source}[a.cmd](a)


if __name__ == "__main__":
    main()
