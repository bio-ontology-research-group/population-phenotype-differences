#!/usr/bin/env python3
"""Score the pipeline against an independent clinician reference (prepare_human_review.py packets).

A reviewer's findings.csv (document_id, patient_label, source_phrase, status, hpo_id, hpo_label, observation_basis)
is compared, per source document, with the common-pipeline output for the same document
(data/v5_private/paired_mapping_final: mentions.jsonl, selected.jsonl) and with the per-class finding classes
(results/termclass.jsonl plus termclass_extra.jsonl).  Reported:

  ids        validity of the reviewer's HPO identifiers against the ontology release
  class      reviewer's observation basis per finding against the pipeline class of the same HPO class
             (bedside = history, examination; test-based = laboratory, imaging, physiological_test, pathology)
  extraction present HPO classes of the pipeline against the reviewer's, per document (exact and lenient:
             equal, ancestor or descendant below organ-system level), micro-averaged
  share      test share per document under the reviewer's basis, under the pipeline classes applied to the
             reviewer's HPO classes, and in the pipeline output
  paired     full text minus abstract, per publication, in the reviewer's reference

  python3 clinician_validation.py --findings data/v5_private/human_review/returned/Paul/findings_PNS.csv \
      --key data/v5_private/human_review/PRIVATE_alignment_key.csv --mapping data/v5_private/paired_mapping_final \
      --termclass results/termclass.jsonl --obo data/hpo/hp.obo --out results/v5/clinician_validation/PNS
"""
import argparse, collections, csv, json, sys
from pathlib import Path
import numpy as np
from rank_engine import parse_obo

TEST = {"laboratory", "imaging", "physiological_test", "pathology"}
BED = {"clinical_exam", "symptom_history"}
BASIS = {"history": "symptom_history", "examination": "clinical_exam", "laboratory": "laboratory", "imaging": "imaging",
         "physiological_test": "physiological_test", "pathology": "pathology"}
POPS = ["EUR", "MENA", "SASIA", "LATAM", "AFR", "EASIA"]


def grp(c):
    return "test" if c in TEST else ("bed" if c in BED else None)


def kappa(pairs):
    n = len(pairs); cats = sorted({a for a, _ in pairs} | {b for _, b in pairs})
    po = sum(a == b for a, b in pairs) / n
    pe = sum((sum(a == c for a, _ in pairs) / n) * (sum(b == c for _, b in pairs) / n) for c in cats)
    return po, (po - pe) / (1 - pe) if pe < 1 else float("nan")


def boot_kappa(pairs, docs, rng, B=2000):
    by = collections.defaultdict(list)
    for p, d in zip(pairs, docs): by[d].append(p)
    keys = sorted(by); out = []
    for _ in range(B):
        s = [p for k in rng.choice(len(keys), len(keys)) for p in by[keys[k]]]
        out.append(kappa(s))
    a = np.array(out); return np.nanpercentile(a[:, 0], [2.5, 97.5]), np.nanpercentile(a[:, 1], [2.5, 97.5])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("findings", "key", "mapping", "termclass", "obo", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(); sys.setrecursionlimit(100000); rng = np.random.default_rng(7)
    terms, alt = parse_obo(a.obo); anc = {}

    def A(t):
        if t not in anc:
            s = {t}
            for p in terms[t]["is_a"]:
                if p in terms: s |= A(p)
            anc[t] = s
        return anc[t]
    top = {"HP:0000001", "HP:0000118"} | {t for t in terms if "HP:0000118" in terms[t]["is_a"]}

    def lenient(x, Y):
        return x in Y or any(y in A(x) - top or x in A(y) - top for y in Y)
    tc = {}
    for f in (Path(a.mapping) / "termclass_extra.jsonl", Path(a.termclass)):
        for x in map(json.loads, open(f)):
            if x.get("class"): tc[x["id"]] = x["class"]
    key = {r["document_id"]: r for r in csv.DictReader(open(a.key))}
    rows = list(csv.DictReader(open(a.findings, encoding="utf-8-sig")))
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True); rep = {}

    # ids
    ids = collections.Counter()
    for r in rows:
        h = r["hpo_id"].strip()
        ids["empty" if not h else "current" if h in terms and not terms[h]["obsolete"] else "alt_id" if h in alt else
            "obsolete" if h in terms else "unknown"] += 1
        r["hpo"] = alt.get(h, h) if (h in terms and not terms[h]["obsolete"]) or h in alt else None
        r["region"] = key[r["document_id"]]["region"]; r["source"] = key[r["document_id"]]["source"]
    rep["rows"] = len(rows); rep["documents"] = len({r["document_id"] for r in rows}); rep["hpo_ids"] = dict(ids)
    rep["status"] = dict(collections.Counter(r["status"] for r in rows))
    rep["basis"] = dict(collections.Counter(r["observation_basis"] for r in rows))

    # class agreement (all findings with a valid HPO class, a reviewer basis and a pipeline class)
    conf = collections.Counter(); pairs = []; pdocs = []; pregion = []; exact = []; detail = []
    for r in rows:
        b = BASIS.get(r["observation_basis"]); c = tc.get(r["hpo"]) if r["hpo"] else None
        conf[(r["observation_basis"], c or ("no_hpo" if not r["hpo"] else "unclassified"))] += 1
        if b and grp(c):
            pairs.append((grp(b), grp(c))); pdocs.append(key[r["document_id"]]["pmid"]); pregion.append(r["region"]); exact.append((b, c))
            if grp(b) != grp(c): detail.append({k: r[k] for k in ("document_id", "source_phrase", "hpo_id", "hpo_label", "observation_basis")} | {"pipeline_class": c})
    po, k = kappa(pairs); ci_po, ci_k = boot_kappa(pairs, pdocs, rng)
    rep["class"] = {"n": len(pairs), "agreement": po, "agreement_ci": list(ci_po), "kappa": k, "kappa_ci": list(ci_k),
                    "exact_class_agreement": sum(x == y for x, y in exact) / len(exact),
                    "reviewer_test_pipeline_bed": sum(p == ("test", "bed") for p in pairs),
                    "reviewer_bed_pipeline_test": sum(p == ("bed", "test") for p in pairs),
                    "by_region": {P: {"n": sum(g == P for g in pregion),
                                      "agreement": float(np.mean([x == y for (x, y), g in zip(pairs, pregion) if g == P])),
                                      "reviewer_test_pipeline_bed": sum(p == ("test", "bed") and g == P for p, g in zip(pairs, pregion)),
                                      "reviewer_bed_pipeline_test": sum(p == ("bed", "test") and g == P for p, g in zip(pairs, pregion))}
                                  for P in POPS}}
    rep["class_confusion"] = {f"{x}|{y}": n for (x, y), n in sorted(conf.items())}
    with open(str(out) + "_class_disagreements.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(detail[0])); w.writeheader(); w.writerows(detail)

    # pipeline output per document
    sel = {x["key"]: alt.get(x["hpo"], x["hpo"]) for x in map(json.loads, open(Path(a.mapping) / "selected.jsonl")) if x.get("ok") and x.get("hpo")}
    pred = collections.defaultdict(set); npred_phr = collections.Counter(); nmapped = collections.Counter()
    for m in map(json.loads, open(Path(a.mapping) / "mentions.jsonl")):
        if m["status"] != "present": continue
        npred_phr[m["document_key"]] += 1
        h = sel.get(m["phrase_key"])
        if h in terms and not terms[h]["obsolete"]: pred[m["document_key"]].add(h); nmapped[m["document_key"]] += 1
    ref = collections.defaultdict(set); refb = collections.defaultdict(list)
    for r in rows:
        if r["status"] != "present": continue
        if r["hpo"]: ref[r["document_id"]].add(r["hpo"])
        if BASIS.get(r["observation_basis"]): refb[r["document_id"]].append(grp(BASIS[r["observation_basis"]]))

    def share(ts):
        t = sum(tc.get(x) in TEST for x in ts); b = sum(tc.get(x) in BED for x in ts)
        return (100 * t / (t + b) if t + b else None), t, b
    docs = []
    for d, k_ in sorted(key.items()):
        dk = f"{k_['pmid']}:{k_['gene']}:{k_['source']}"; R = ref[d]; P = pred[dk]
        sr, tr, br = share(R); sp, tp, bp = share(P)
        nb = len(refb[d]); sb = 100 * sum(x == "test" for x in refb[d]) / nb if nb else None
        docs.append({"document_id": d, "pmid": k_["pmid"], "region": k_["region"], "source": k_["source"],
                     "ref_terms": len(R), "pred_terms": len(P), "exact_tp": len(R & P),
                     "lenient_ref_found": sum(lenient(x, P) for x in R), "lenient_pred_supported": sum(lenient(x, R) for x in P),
                     "ref_test": tr, "ref_bed": br, "pred_test": tp, "pred_bed": bp,
                     "share_ref_basis": sb, "share_ref_pipeline_class": sr, "share_pred": sp})
    with open(str(out) + "_documents.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(docs[0])); w.writeheader(); w.writerows(docs)

    def prf(D):
        R = sum(d["ref_terms"] for d in D); P = sum(d["pred_terms"] for d in D)
        return {"documents": len(D), "ref_terms": R, "pred_terms": P,
                "exact_precision": sum(d["exact_tp"] for d in D) / P if P else None, "exact_recall": sum(d["exact_tp"] for d in D) / R if R else None,
                "lenient_precision": sum(d["lenient_pred_supported"] for d in D) / P if P else None,
                "lenient_recall": sum(d["lenient_ref_found"] for d in D) / R if R else None}
    rep["extraction"] = {"all": prf(docs)} | {s: prf([d for d in docs if d["source"] == s]) for s in ("abstract", "fulltext")} | \
        {P: prf([d for d in docs if d["region"] == P]) for P in POPS}

    def typed_recall(D, T):
        num = den = 0
        for d in D:
            dk = f"{key[d['document_id']]['pmid']}:{key[d['document_id']]['gene']}:{d['source']}"
            for x in ref[d["document_id"]]:
                if tc.get(x) in T: den += 1; num += lenient(x, pred[dk])
        return {"n": den, "lenient_recall": num / den if den else None}
    rep["recall_by_type"] = {s: {"test": typed_recall([d for d in docs if d["source"] == s], TEST), "bed": typed_recall([d for d in docs if d["source"] == s], BED)}
                             for s in ("abstract", "fulltext")}

    def diffs(D, x, y):
        v = np.array([d[x] - d[y] for d in D if d[x] is not None and d[y] is not None])
        if len(v) < 2: return {"n": len(v)}
        b = [rng.choice(v, len(v)).mean() for _ in range(2000)]
        return {"n": len(v), "mean_diff": float(v.mean()), "ci": [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]}
    rep["share"] = {"pipeline_minus_reviewer_basis": {s: diffs([d for d in docs if s in ("all", d["source"])], "share_pred", "share_ref_basis") for s in ("all", "abstract", "fulltext")},
                    "pipeline_class_minus_reviewer_basis_on_reviewer_terms": diffs(docs, "share_ref_pipeline_class", "share_ref_basis"),
                    "pipeline_minus_reviewer_basis_by_region": {P: diffs([d for d in docs if d["region"] == P], "share_pred", "share_ref_basis") for P in POPS},
                    "mean_share_reviewer_basis": {s: float(np.mean([d["share_ref_basis"] for d in docs if d["source"] == s and d["share_ref_basis"] is not None])) for s in ("abstract", "fulltext")},
                    "mean_share_pipeline": {s: float(np.mean([d["share_pred"] for d in docs if d["source"] == s and d["share_pred"] is not None])) for s in ("abstract", "fulltext")}}
    x = [(d["share_pred"], d["share_ref_basis"]) for d in docs if d["share_pred"] is not None and d["share_ref_basis"] is not None]
    rep["share"]["pearson_r"] = float(np.corrcoef(*zip(*x))[0, 1]); rep["share"]["n_r"] = len(x)

    # paired full text minus abstract in the reviewer's reference
    by = collections.defaultdict(dict)
    for d in docs: by[d["pmid"]][d["source"]] = d
    pr = [{"ft_terms": v["fulltext"]["ref_terms"], "ab_terms": v["abstract"]["ref_terms"], "ft_share": v["fulltext"]["share_ref_basis"], "ab_share": v["abstract"]["share_ref_basis"]} for v in by.values() if len(v) == 2]
    rep["paired_reviewer"] = {"publications": len(pr), "terms_fulltext_minus_abstract": diffs(pr, "ft_terms", "ab_terms"), "share_fulltext_minus_abstract": diffs(pr, "ft_share", "ab_share")}
    json.dump(rep, open(str(out) + ".json", "w"), indent=1); print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
