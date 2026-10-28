#!/usr/bin/env python3
"""Analyses for the PAVS prioritization benchmark paper.

  validate  : compare the vectorised engine with the per-case ranks in the PAVS repo
  baseline  : extrinsic/intrinsic Resnik ranks for all evaluable cases (negated terms dropped)
  lpo       : leave-publication-out circularity control for Phenopacket Store cases
  depth     : depth-matched re-evaluation of Store cases (subsample to clinical term counts)

All outputs go to --out (CSV). Reference files come from one HPO release (--hpo-dir).
"""
import argparse, json, os, re, sys, csv, random, collections
import numpy as np
from rank_engine import Engine, metrics

random.seed(20261024)
np.random.seed(20261024)


def load_pavs(path, keep_excluded=False):
    pp = json.load(open(path))
    out = []
    for p in pp:
        pid = p["id"]; src = pid[5]
        if src not in "ABFPMQ":
            continue
        hpos = [f["type"]["id"] for f in p.get("phenotypicFeatures", []) if keep_excluded or not f.get("excluded")]
        gene = None
        try:
            gene = p["interpretations"][0]["diagnosis"]["genomicInterpretations"][0]["variantInterpretation"]["variationDescriptor"]["geneContext"]["symbol"]
        except (KeyError, IndexError):
            pass
        pmid = None
        for er in p.get("metaData", {}).get("externalReferences", []):
            m = re.match(r"PMID:(\d+)", er.get("id", ""))
            if m:
                pmid = m.group(1)
        cohort = {"A": "Saudi-clinical", "B": "Saudi-clinical", "F": "Saudi-clinical", "P": "Saudi-clinical",
                  "M": "Saudi-literature", "Q": "US-clinical"}[src]
        out.append({"id": pid, "cohort": cohort, "source": src, "hpos": hpos, "gene": gene, "pmid": pmid})
    return out


def load_store(root, keep_excluded=False):
    out = []
    for dp, _, fs in os.walk(root):
        for f in fs:
            if not f.endswith(".json"):
                continue
            try:
                p = json.load(open(os.path.join(dp, f)))
            except Exception:
                continue
            hpos = [x["type"]["id"] for x in p.get("phenotypicFeatures", []) if keep_excluded or not x.get("excluded")]
            gene = None
            for i in p.get("interpretations", []):
                for gi in i.get("diagnosis", {}).get("genomicInterpretations", []):
                    gene = gi.get("variantInterpretation", {}).get("variationDescriptor", {}).get("geneContext", {}).get("symbol")
                    if gene:
                        break
                if gene:
                    break
            pmid = None
            for er in p.get("metaData", {}).get("externalReferences", []):
                m = re.match(r"PMID:(\d+)", er.get("id", ""))
                if m:
                    pmid = m.group(1)
            if pmid is None:
                m = re.match(r"PMID_(\d+)_", p.get("id", ""))
                if m:
                    pmid = m.group(1)
            out.append({"id": p.get("id", f), "cohort": "Literature", "source": "Store", "hpos": hpos, "gene": gene, "pmid": pmid,
                        "cohort_dir": os.path.basename(dp)})
    return out


def evaluable(cases, eng):
    out = []
    for c in cases:
        hp = eng.map_terms(c["hpos"])
        if hp and c["gene"] in eng.gid:
            c = dict(c); c["terms"] = hp
            out.append(c)
    return out


def mean_ic(eng, terms):
    return float(np.mean([eng.ic[eng.tid[t]] for t in terms]))


def write_rows(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})


def cmd_validate(a, eng, pavs, store):
    ref = {}
    for r in csv.DictReader(open(a.ref_results)):
        ref[r["pavs_id"]] = (r["gene"], int(r["rank"]), r["group"])
    pool = [c for c in pavs + store if c["id"] in ref]
    random.shuffle(pool)
    sample = pool[: a.n]
    exact = 0; rows = []
    for c in sample:
        hp = eng.map_terms(c["hpos"])
        if not hp or c["gene"] not in eng.gid:
            continue
        rank, score, auc = eng.rank_gene(hp, c["gene"])
        rows.append({"id": c["id"], "gene": c["gene"], "ref_gene": ref[c["id"]][0], "ref_rank": ref[c["id"]][1], "rank": rank, "group": ref[c["id"]][2]})
        exact += rank == ref[c["id"]][1]
    rr = np.array([r["ref_rank"] for r in rows]); mr = np.array([r["rank"] for r in rows])
    print(f"validate: {len(rows)} cases, exact rank match {exact} ({100*exact/len(rows):.1f}%), "
          f"same-gene {sum(r['gene']==r['ref_gene'] for r in rows)}, hits1 ref {np.mean(rr==1):.3f} new {np.mean(mr==1):.3f}, "
          f"MRR ref {np.mean(1/rr):.4f} new {np.mean(1/mr):.4f}, spearman {np.corrcoef(np.argsort(np.argsort(rr)), np.argsort(np.argsort(mr)))[0,1]:.4f}")
    write_rows(os.path.join(a.out, "validate.csv"), rows, ["id", "group", "gene", "ref_gene", "ref_rank", "rank"])


def cmd_baseline(a, eng, pavs, store):
    rows = []
    for ic_type in ("extrinsic", "intrinsic"):
        eng.set_ic(ic_type)
        for c in evaluable(pavs + store, eng):
            rank, score, auc = eng.rank_gene(c["terms"], c["gene"])
            rows.append({"id": c["id"], "cohort": c["cohort"], "source": c["source"], "pmid": c["pmid"], "gene": c["gene"],
                         "n_terms": len(c["terms"]), "mean_ic": round(mean_ic(eng, c["terms"]), 4), "ic_type": ic_type,
                         "rank": rank, "score": round(score, 5), "auc": round(auc, 5)})
        print(f"baseline {ic_type}: {sum(1 for r in rows if r['ic_type']==ic_type)} cases", flush=True)
    write_rows(os.path.join(a.out, "baseline_percase.csv"), rows,
               ["id", "cohort", "source", "pmid", "gene", "n_terms", "mean_ic", "ic_type", "rank", "score", "auc"])
    summ = []
    for ic_type in ("extrinsic", "intrinsic"):
        for cohort in sorted({r["cohort"] for r in rows}):
            sub = [r for r in rows if r["ic_type"] == ic_type and r["cohort"] == cohort]
            m = metrics([r["rank"] for r in sub], [r["auc"] for r in sub])
            m.update({"ic_type": ic_type, "cohort": cohort, "mean_terms": round(np.mean([r["n_terms"] for r in sub]), 2),
                      "mean_ic": round(np.mean([r["mean_ic"] for r in sub]), 3)})
            summ.append(m)
            print(f"  {ic_type:9s} {cohort:16s} n={m['n']:5d} Hits@1={m['hits1']:.2f} Hits@10={m['hits10']:.2f} MRR={m['mrr']:.4f} AUC={m['auc']:.4f} terms={m['mean_terms']} IC={m['mean_ic']}")
    write_rows(os.path.join(a.out, "baseline_summary.csv"), summ, ["ic_type", "cohort", "n", "hits1", "hits10", "mrr", "median_rank", "auc", "mean_terms", "mean_ic"])


def build_support(eng, hpoa, g2p):
    """(gene, hpo) -> set of references supporting it, via g2p disease_id and phenotype.hpoa."""
    dref = collections.defaultdict(set)  # (disease, hpo) -> refs
    for line in open(hpoa):
        if line.startswith("#") or line.startswith("database_id"):
            continue
        c = line.rstrip("\n").split("\t")
        for r in c[4].split(";"):
            r = r.strip()
            if r:
                dref[(c[0], c[3])].add(r)
    support = collections.defaultdict(set)
    for line in open(g2p):
        if line.startswith("ncbi_gene_id") or line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        gene, hpo, dis = c[1], c[2], c[5] if len(c) > 5 else ""
        support[(gene, hpo)] |= dref.get((dis, hpo), set())
    return support


def gene_score_sets(eng, P_terms, G_terms):
    """BMA Resnik for one gene profile given as term list (python-set version of the engine)."""
    if not G_terms:
        return 0.0
    ancP = set();
    for p in P_terms:
        ancP |= eng.anc_sets[p]
    ancG = set()
    for g in G_terms:
        ancG |= eng.anc_sets[g]
    s1 = np.mean([max((eng.ic[eng.tid[a]] for a in eng.anc_sets[p] & ancG), default=0.0) for p in P_terms])
    s2 = np.mean([max((eng.ic[eng.tid[a]] for a in eng.anc_sets[g] & ancP), default=0.0) for g in G_terms])
    return (s1 + s2) / 2


def cmd_lpo(a, eng, pavs, store):
    eng.set_ic("extrinsic")
    support = build_support(eng, os.path.join(a.hpo_dir, "phenotype.hpoa"), os.path.join(a.hpo_dir, "genes_to_phenotype.txt"))
    # pairs solely supported by one PMID
    sole = collections.defaultdict(list)  # pmid -> [(gene,hpo)]
    for (g, h), refs in support.items():
        pm = {r[5:] for r in refs if r.startswith("PMID:")}
        if len(pm) == 1 and len(refs) == 1:
            sole[next(iter(pm))].append((g, h))
    pm_any = collections.defaultdict(set)  # pmid -> genes it supports at all
    for (g, h), refs in support.items():
        for r in refs:
            if r.startswith("PMID:"):
                pm_any[r[5:]].add(g)
    print(f"lpo: {len(support)} gene-term pairs, {sum(len(v) for v in sole.values())} solely supported by one PMID across {len(sole)} PMIDs", flush=True)
    rows = []
    cases = evaluable(store, eng)
    for i, c in enumerate(cases):
        full_scores = eng.bma_all_genes(c["terms"])
        order = np.lexsort((np.array(eng.genes), -full_scores))
        rank_full = int(np.where(order == eng.gid[c["gene"]])[0][0] + 1)
        pm = c["pmid"]
        in_hpoa = pm in pm_any
        supports_gene = c["gene"] in pm_any.get(pm, set())
        removed = sole.get(pm, [])
        affected = collections.defaultdict(set)
        for g, h in removed:
            affected[g].add(h)
        gene_terms = set(eng.g2t[c["gene"]])
        case_terms_in_profile = sum(1 for t in c["terms"] if t in gene_terms)
        case_terms_sole = sum(1 for t in c["terms"] if t in affected.get(c["gene"], set()))
        rank_lpo = rank_full; score_lpo = float(full_scores[eng.gid[c["gene"]]])
        if affected:
            scores = full_scores.copy()
            for g, hs in affected.items():
                if g not in eng.gid:
                    continue
                newG = [t for t in eng.g2t[g] if t not in hs]
                scores[eng.gid[g]] = gene_score_sets(eng, c["terms"], newG)
            order = np.lexsort((np.array(eng.genes), -scores))
            rank_lpo = int(np.where(order == eng.gid[c["gene"]])[0][0] + 1)
            score_lpo = float(scores[eng.gid[c["gene"]]])
        rows.append({"id": c["id"], "pmid": pm, "gene": c["gene"], "n_terms": len(c["terms"]), "pmid_in_hpoa": int(in_hpoa),
                     "pmid_supports_gene": int(supports_gene), "n_gene_terms": len(gene_terms),
                     "n_gene_terms_removed": len(affected.get(c["gene"], set())), "case_terms_in_profile": case_terms_in_profile,
                     "case_terms_sole_support": case_terms_sole, "rank_full": rank_full, "rank_lpo": rank_lpo,
                     "score_full": round(float(full_scores[eng.gid[c["gene"]]]), 5), "score_lpo": round(score_lpo, 5)})
        if i % 1000 == 0:
            print(f"lpo: {i}/{len(cases)}", flush=True)
    write_rows(os.path.join(a.out, "lpo_percase.csv"), rows, list(rows[0].keys()))
    def summ(name, sub):
        if not sub:
            return
        mf = metrics([r["rank_full"] for r in sub]); ml = metrics([r["rank_lpo"] for r in sub])
        print(f"  {name:38s} n={len(sub):5d} Hits@1 full {mf['hits1']:.2f} -> LPO {ml['hits1']:.2f}; MRR {mf['mrr']:.4f} -> {ml['mrr']:.4f}; "
              f"Hits@10 {mf['hits10']:.2f} -> {ml['hits10']:.2f}")
    summ("all Store cases", rows)
    summ("PMID in HPOA (any gene)", [r for r in rows if r["pmid_in_hpoa"]])
    summ("PMID supports the causative gene", [r for r in rows if r["pmid_supports_gene"]])
    summ("PMID solely supports >=1 term of gene", [r for r in rows if r["n_gene_terms_removed"] > 0])
    summ("PMID not in HPOA", [r for r in rows if not r["pmid_in_hpoa"]])


def cmd_depth(a, eng, pavs, store):
    eng.set_ic("extrinsic")
    clin = [c for c in evaluable(pavs, eng) if c["cohort"] == "Saudi-clinical"]
    target = [len(c["terms"]) for c in clin]
    cases = evaluable(store, eng)
    rows = []
    for i, c in enumerate(cases):
        T = c["terms"]
        rank_full, _, _ = eng.rank_gene(T, c["gene"])
        for rep in range(a.reps):
            n = random.choice(target)
            if n >= len(T):
                sub_rand = T; sub_top = T
            else:
                sub_rand = random.sample(T, n)
                sub_top = sorted(T, key=lambda t: -eng.ic[eng.tid[t]])[:n]
            r_rand, _, _ = eng.rank_gene(sub_rand, c["gene"])
            r_top, _, _ = eng.rank_gene(sub_top, c["gene"])
            rows.append({"id": c["id"], "pmid": c["pmid"], "gene": c["gene"], "n_full": len(T), "n_sub": min(n, len(T)), "rep": rep,
                         "rank_full": rank_full, "rank_random_subset": r_rand, "rank_topic_subset": r_top,
                         "mean_ic_full": round(mean_ic(eng, T), 3), "mean_ic_random": round(mean_ic(eng, sub_rand), 3), "mean_ic_top": round(mean_ic(eng, sub_top), 3)})
        if i % 1000 == 0:
            print(f"depth: {i}/{len(cases)}", flush=True)
    write_rows(os.path.join(a.out, "depth_percase.csv"), rows, list(rows[0].keys()))
    for k in ("rank_full", "rank_random_subset", "rank_topic_subset"):
        m = metrics([r[k] for r in rows])
        print(f"  Store {k:20s} Hits@1 {m['hits1']:.2f} Hits@10 {m['hits10']:.2f} MRR {m['mrr']:.4f}")
    m = metrics([eng.rank_gene(c['terms'], c['gene'])[0] for c in clin])
    print(f"  Saudi-clinical reference      Hits@1 {m['hits1']:.2f} Hits@10 {m['hits10']:.2f} MRR {m['mrr']:.4f} (mean terms {np.mean(target):.2f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["validate", "baseline", "lpo", "depth"])
    ap.add_argument("--hpo-dir", required=True)
    ap.add_argument("--pavs", required=True)
    ap.add_argument("--store", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref-results")
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--keep-excluded", action="store_true", help="count negated features as present (replicates the original script)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    eng = Engine(os.path.join(a.hpo_dir, "hp.obo"), os.path.join(a.hpo_dir, "genes_to_phenotype.txt"))
    print(f"engine: {len(eng.terms)} terms, {len(eng.genes)} genes", flush=True)
    pavs = load_pavs(a.pavs, a.keep_excluded)
    store = load_store(a.store, a.keep_excluded)
    print(f"cases: PAVS {len(pavs)}, Store {len(store)}", flush=True)
    {"validate": cmd_validate, "baseline": cmd_baseline, "lpo": cmd_lpo, "depth": cmd_depth}[a.cmd](a, eng, pavs, store)


if __name__ == "__main__":
    main()
