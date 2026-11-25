#!/usr/bin/env python3
"""Agreement between two reviewers of the same packets (prepare_human_review.py), and of the pipeline with both.

Inputs are the normalised findings sheets (clinician_review_normalize.py).  Per source document:

  terms      present HPO classes of reviewer A against reviewer B: exact and lenient (equal, ancestor or
             descendant below organ-system level) overlap, micro-averaged, as positive agreement
             2|A and B| / (|A| + |B|); the same for the pipeline against each reviewer, against the classes
             both reviewers recorded (consensus) and against the classes either recorded (union)
  status     present against explicitly absent, for classes both reviewers recorded in a document
  basis      bedside (history, examination) against test-based (laboratory, imaging, physiological test,
             pathology) for classes both recorded as present with a single basis: agreement and Cohen's kappa,
             publication bootstrap; and the pipeline class against the basis on which both reviewers agree
  share      test share per document under each reviewer's basis: correlation and mean difference, and
             pipeline minus the mean of the two reviewers, overall, by source and by regional group

  python3 clinician_agreement.py --a A.norm.csv --b B.norm.csv --key KEY.csv --mapping DIR --termclass F --obo F --out PREFIX
"""
import argparse, collections, csv, json, sys
from pathlib import Path
import numpy as np
from rank_engine import parse_obo
from clinician_validation import TEST, BED, BASIS, POPS, grp, kappa, boot_kappa


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("a", "b", "key", "mapping", "termclass", "obo", "out"): ap.add_argument("--" + k, required=True)
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

    def len_hits(X, Y):
        return sum(x in Y or any(y in A(x) - top or x in A(y) - top for y in Y) for x in X)
    tc = {}
    for f in (Path(a.mapping) / "termclass_extra.jsonl", Path(a.termclass)):
        for x in map(json.loads, open(f)):
            if x.get("class"): tc[x["id"]] = x["class"]
    key = {r["document_id"]: r for r in csv.DictReader(open(a.key))}

    def load(f):
        st = collections.defaultdict(dict); ba = collections.defaultdict(lambda: collections.defaultdict(set))
        for r in csv.DictReader(open(f, encoding="utf-8-sig")):
            h = r["hpo_id"].strip(); h = alt.get(h, h)
            if r["document_id"] not in key or h not in terms or terms[h]["obsolete"]: continue
            d = r["document_id"]
            if r["status"] == "present" or h not in st[d]: st[d][h] = r["status"]
            if r["status"] == "present" and BASIS.get(r["observation_basis"]): ba[d][h].add(grp(BASIS[r["observation_basis"]]))
        return st, ba
    (sa, ba), (sb, bb) = load(a.a), load(a.b)
    # pipeline output
    sel = {x["key"]: alt.get(x["hpo"], x["hpo"]) for x in map(json.loads, open(Path(a.mapping) / "selected.jsonl")) if x.get("ok") and x.get("hpo")}
    pipe = collections.defaultdict(set); doc_of = {f'{k["pmid"]}:{k["gene"]}:{k["source"]}': d for d, k in key.items()}
    for m in map(json.loads, open(Path(a.mapping) / "mentions.jsonl")):
        d = doc_of.get(m["document_key"]); h = sel.get(m["phrase_key"])
        if d and m["status"] == "present" and h in terms and not terms[h]["obsolete"]: pipe[d].add(h)
    docs = sorted(key); P = lambda s, d: {h for h, v in s[d].items() if v == "present"}
    rep = {"documents_a": len([d for d in docs if sa[d]]), "documents_b": len([d for d in docs if sb[d]])}

    def overlap(X, Y, sub=docs):
        nx = sum(len(X(d)) for d in sub); ny = sum(len(Y(d)) for d in sub)
        ex = sum(len(X(d) & Y(d)) for d in sub); lx = sum(len_hits(X(d), Y(d)) for d in sub); ly = sum(len_hits(Y(d), X(d)) for d in sub)
        return {"n_x": nx, "n_y": ny, "exact_x_in_y": ex / nx if nx else None, "exact_y_in_x": ex / ny if ny else None, "exact_positive_agreement": 2 * ex / (nx + ny),
                "lenient_x_in_y": lx / nx if nx else None, "lenient_y_in_x": ly / ny if ny else None, "lenient_positive_agreement": (lx + ly) / (nx + ny)}
    both = [d for d in docs if P(sa, d) and P(sb, d)]
    pa, pb = (lambda d: P(sa, d)), (lambda d: P(sb, d)); pp = lambda d: pipe[d]
    cons = lambda d: {x for x in P(sa, d) if len_hits({x}, P(sb, d))} | {x for x in P(sb, d) if len_hits({x}, P(sa, d))}
    uni = lambda d: P(sa, d) | P(sb, d)
    rep["terms"] = {"documents_both": len(both)}
    for name, sub in (("all", both), ("abstract", [d for d in both if key[d]["source"] == "abstract"]), ("fulltext", [d for d in both if key[d]["source"] == "fulltext"])):
        rep["terms"][name] = {"a_vs_b": overlap(pa, pb, sub), "pipeline_vs_a": overlap(pp, pa, sub), "pipeline_vs_b": overlap(pp, pb, sub),
                              "pipeline_vs_consensus": overlap(pp, cons, sub), "pipeline_vs_union": overlap(pp, uni, sub)}
    # status and basis on shared classes
    st = [(sa[d][h], sb[d][h]) for d in docs for h in set(sa[d]) & set(sb[d]) if {sa[d][h], sb[d][h]} <= {"present", "explicitly_absent"}]
    rep["status"] = {"n": len(st), "agreement": kappa(st)[0], "kappa": kappa(st)[1]}
    pairs, pd, trip = [], [], []
    for d in docs:
        for h in set(ba[d]) & set(bb[d]):
            if len(ba[d][h]) == 1 and len(bb[d][h]) == 1:
                x, y = next(iter(ba[d][h])), next(iter(bb[d][h])); pairs.append((x, y)); pd.append(key[d]["pmid"])
                if grp(tc.get(h)): trip.append((x, y, grp(tc[h]), key[d]["pmid"], key[d]["region"]))
    po, k = kappa(pairs); ci_po, ci_k = boot_kappa(pairs, pd, rng)
    rep["basis_a_vs_b"] = {"n": len(pairs), "agreement": po, "agreement_ci": list(ci_po), "kappa": k, "kappa_ci": list(ci_k)}
    ag = [(x, c, pm, rg) for x, y, c, pm, rg in trip if x == y]
    po, k = kappa([(x, c) for x, c, _, _ in ag]); ci_po, ci_k = boot_kappa([(x, c) for x, c, _, _ in ag], [pm for _, _, pm, _ in ag], rng)
    rep["pipeline_class_vs_agreed_basis"] = {"n": len(ag), "agreement": po, "agreement_ci": list(ci_po), "kappa": k, "kappa_ci": list(ci_k),
        "reviewers_bed_pipeline_test": sum(x == "bed" and c == "test" for x, c, _, _ in ag), "reviewers_test_pipeline_bed": sum(x == "test" and c == "bed" for x, c, _, _ in ag),
        "by_region": {R: {"n": sum(r == R for *_, r in ag), "agreement": float(np.mean([x == c for x, c, _, r in ag if r == R])) if any(r == R for *_, r in ag) else None,
                          "reviewers_bed_pipeline_test": sum(x == "bed" and c == "test" and r == R for x, c, _, r in ag)} for R in POPS}}
    dis = [(c) for x, y, c, _, _ in trip if x != y]
    rep["pipeline_class_where_reviewers_disagree"] = {"n": len(dis), "pipeline_test": sum(c == "test" for c in dis)}

    # test share per document
    def share(b, d):
        v = [next(iter(s)) for s in b[d].values() if len(s) == 1]
        return 100 * sum(x == "test" for x in v) / len(v) if len(v) >= 3 else None
    def pshare(d):
        v = [grp(tc.get(h)) for h in pipe[d]]; v = [x for x in v if x]
        return 100 * sum(x == "test" for x in v) / len(v) if len(v) >= 3 else None
    rows = [{"document_id": d, "pmid": key[d]["pmid"], "source": key[d]["source"], "region": key[d]["region"], "a": share(ba, d), "b": share(bb, d), "pipeline": pshare(d)} for d in docs]

    def diff(rs, f):
        by = collections.defaultdict(list)
        for r in rs:
            v = f(r)
            if v is not None: by[r["pmid"]].append(v)
        ks = sorted(by)
        if not ks: return {"n": 0}
        m = float(np.mean([v for k in ks for v in by[k]])); n = sum(len(by[k]) for k in ks)
        bs = [np.mean([v for i in rng.choice(len(ks), len(ks)) for v in by[ks[i]]]) for _ in range(2000)]
        return {"n": n, "mean_diff": m, "ci": [float(x) for x in np.percentile(bs, [2.5, 97.5])]}
    sub = lambda f: {"all": diff(rows, f), "abstract": diff([r for r in rows if r["source"] == "abstract"], f), "fulltext": diff([r for r in rows if r["source"] == "fulltext"], f),
                     "by_region": {R: diff([r for r in rows if r["region"] == R], f) for R in POPS},
                     "by_region_fulltext": {R: diff([r for r in rows if r["region"] == R and r["source"] == "fulltext"], f) for R in POPS}}
    ok = [r for r in rows if r["a"] is not None and r["b"] is not None]
    rep["share"] = {"n_both": len(ok), "pearson_a_b": float(np.corrcoef([r["a"] for r in ok], [r["b"] for r in ok])[0, 1]),
                    "a_minus_b": sub(lambda r: r["a"] - r["b"] if None not in (r["a"], r["b"]) else None),
                    "pipeline_minus_mean_reviewer": sub(lambda r: r["pipeline"] - (r["a"] + r["b"]) / 2 if None not in (r["a"], r["b"], r["pipeline"]) else None)}
    okp = [r for r in ok if r["pipeline"] is not None]
    rep["share"]["pearson_pipeline_mean_reviewer"] = float(np.corrcoef([r["pipeline"] for r in okp], [(r["a"] + r["b"]) / 2 for r in okp])[0, 1])
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    with open(str(out) + "_documents.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    json.dump(rep, open(str(out) + ".json", "w"), indent=1); print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
