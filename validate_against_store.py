#!/usr/bin/env python3
"""Extraction performance against Phenopacket Store curation.

The Store curators assigned HPO classes to each patient from the full text of the publication. Our pipeline
(language-model extraction and HPO mapping) was run on the same publications, either on their PubMed abstracts
(--pred-corpus abstracts) or on their PMC-Patients summaries. Units:
  publication x gene   union of the observed HPO classes of all patients with that gene in the publication;
  patient              publications in which both sides have exactly one patient with that gene.
Metrics per unit, then averaged (macro) and pooled (micro):
  exact     precision = |P & G| / |P|, recall = |P & G| / |G|;
  lenient   a predicted class counts as correct when it equals, or is an ancestor or a descendant of, a curated class
            (and symmetrically for recall), ignoring matches through classes at or above the organ-system level;
  F1 from precision and recall.
Breakdowns: by the Store population label of the patients, and by finding class (test-based, sign, symptom;
termclass.jsonl) of the classes. Also reported: the share of test-based findings on both sides and its difference,
by population, because a population-dependent bias in this share would bias the main result.

  validate_against_store.py --store-cases fine2/cases_store.tsv --pred-cases P.tsv --obo hp.obo --termclass T --out O
"""
import argparse, collections, csv, json
import numpy as np

from rank_engine import parse_obo

TEST = {"laboratory", "imaging", "physiological_test", "pathology"}
POPS = ["EUR", "MENA", "SASIA", "LATAM", "AFR", "EASIA"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store-cases", required=True); ap.add_argument("--pred-cases", required=True)
    ap.add_argument("--obo", required=True); ap.add_argument("--termclass", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    terms, alt = parse_obo(a.obo)
    anc = {}

    def A(t):
        if t in anc: return anc[t]
        s = {t}
        for p in terms[t]["is_a"]:
            if p in terms: s |= A(p)
        anc[t] = s; return s
    import sys; sys.setrecursionlimit(100000)
    top = {"HP:0000001", "HP:0000118"} | {t for t in terms if "HP:0000118" in terms[t]["is_a"]}
    tc = {x["id"]: x["class"] for x in map(json.loads, open(a.termclass))}

    def kind(t):
        c = tc.get(t)
        return "test" if c in TEST else ("sign" if c == "clinical_exam" else ("symptom" if c == "symptom_history" else "other"))

    def load(path):
        by = collections.defaultdict(list)
        for r in csv.DictReader(open(path), delimiter="\t"):
            ts = {alt.get(t, t) for t in r["terms"].split(";") if t}
            ts = {t for t in ts if t in terms and not terms[t]["obsolete"]}
            if ts: by[(r["pmid"], r["gene"])].append((ts, r["region"]))
        return by
    G = load(a.store_cases); P = load(a.pred_cases)

    def match(x, Y):
        if x in Y: return True
        ax = A(x) - top
        return any((y in ax) or (x in (A(y) - top)) for y in Y)
    units = []
    for k in sorted(set(G) & set(P)):
        reg = collections.Counter(r for _, r in G[k]).most_common(1)[0][0]
        for level, g, p in (("publication", set().union(*(t for t, _ in G[k])), set().union(*(t for t, _ in P[k]))),) + \
                ((("patient", G[k][0][0], P[k][0][0]),) if len(G[k]) == 1 and len(P[k]) == 1 else ()):
            u = {"pmid": k[0], "gene": k[1], "level": level, "region": reg, "n_store": len(g), "n_pred": len(p)}
            u["tp_exact"] = len(g & p)
            u["tp_len_p"] = sum(match(x, g) for x in p); u["tp_len_g"] = sum(match(x, p) for x in g)
            for kd in ("test", "sign", "symptom"):
                gk = {t for t in g if kind(t) == kd}; pk = {t for t in p if kind(t) == kd}
                u[f"g_{kd}"] = len(gk); u[f"p_{kd}"] = len(pk)
                u[f"rec_len_{kd}"] = sum(match(x, p) for x in gk); u[f"prec_len_{kd}"] = sum(match(x, g) for x in pk)
            nt_g = u["g_test"]; nb_g = u["g_sign"] + u["g_symptom"]; nt_p = u["p_test"]; nb_p = u["p_sign"] + u["p_symptom"]
            u["share_store"] = nt_g / (nt_g + nb_g) if nt_g + nb_g else float("nan")
            u["share_pred"] = nt_p / (nt_p + nb_p) if nt_p + nb_p else float("nan")
            units.append(u)
    with open(a.out + "_units.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(units[0])); w.writeheader(); w.writerows(units)
    rows = []

    def summ(level, name, U):
        if not U: return
        def micro(num, den): s = sum(u[den] for u in U); return sum(u[num] for u in U) / s if s else float("nan")
        r = {"level": level, "group": name, "units": len(U), "publications": len({u["pmid"] for u in U}),
             "store_classes": np.mean([u["n_store"] for u in U]), "pred_classes": np.mean([u["n_pred"] for u in U]),
             "precision_exact": micro("tp_exact", "n_pred"), "recall_exact": micro("tp_exact", "n_store"),
             "precision_lenient": micro("tp_len_p", "n_pred"), "recall_lenient": micro("tp_len_g", "n_store")}
        for m in ("exact", "lenient"):
            p_, r_ = r[f"precision_{m}"], r[f"recall_{m}"]; r[f"f1_{m}"] = 2 * p_ * r_ / (p_ + r_) if p_ + r_ else float("nan")
        for kd in ("test", "sign", "symptom"):
            r[f"recall_lenient_{kd}"] = micro(f"rec_len_{kd}", f"g_{kd}"); r[f"precision_lenient_{kd}"] = micro(f"prec_len_{kd}", f"p_{kd}")
        d = [u["share_pred"] - u["share_store"] for u in U if u["share_pred"] == u["share_pred"] and u["share_store"] == u["share_store"]]
        r["test_share_store"] = np.nanmean([u["share_store"] for u in U]); r["test_share_pred"] = np.nanmean([u["share_pred"] for u in U])
        r["test_share_diff"] = np.mean(d) if d else float("nan"); r["test_share_n"] = len(d)
        rows.append(r)
        print(f"{level:11s} {name:6s} units {len(U):4d}  P/R/F1 exact {r['precision_exact']:.2f}/{r['recall_exact']:.2f}/{r['f1_exact']:.2f}  "
              f"lenient {r['precision_lenient']:.2f}/{r['recall_lenient']:.2f}/{r['f1_lenient']:.2f}  recall test/sign/symptom "
              f"{r['recall_lenient_test']:.2f}/{r['recall_lenient_sign']:.2f}/{r['recall_lenient_symptom']:.2f}  "
              f"classes store {r['store_classes']:.1f} pred {r['pred_classes']:.1f}  test share store {r['test_share_store']:.3f} "
              f"pred {r['test_share_pred']:.3f} diff {r['test_share_diff']:+.3f}")
    for level in ("publication", "patient"):
        U = [u for u in units if u["level"] == level]
        summ(level, "all", U)
        for g in POPS: summ(level, g, [u for u in U if u["region"] == g])
    with open(a.out + "_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
