#!/usr/bin/env python3
"""Within-curator population contrast on the Phenopacket Store.

Labels: --labels is a TSV with columns pmid, region (MENA, EUR, EASIA, OTHER, UNKNOWN)
and optionally n_patients. Two sources:
  make-labels-proxy : from phenopacket_store_0.1.26_pmid_country.tsv (first-author country)
  make-labels-llm   : from origin.jsonl produced by store_origin_extract.py (patient origin)

Analyses (all on Store cases, extrinsic Resnik, HPO 2026-01-08):
  rank    : ranking metrics by region, full and LPO reference, with depth
  profile : per-gene term-profile divergence MENA vs EUR with a paper-level permutation null
  popref  : population-specific reference test: for each gene with enough cases on both sides,
            rank held-out MENA cases against (a) a profile built from other MENA papers,
            (b) a size-matched profile from EUR papers, (c) HPOA; leave-one-paper-out
"""
import argparse, csv, json, os, re, random, collections
import numpy as np
from rank_engine import Engine, metrics
from run_analyses import load_store, evaluable

random.seed(7); np.random.seed(7)

MENA = {"Saudi Arabia", "United Arab Emirates", "Qatar", "Kuwait", "Oman", "Bahrain", "Yemen", "Jordan", "Lebanon", "Syria", "Iraq",
        "Egypt", "Morocco", "Tunisia", "Algeria", "Libya", "Sudan", "Iran", "Turkey", "Türkiye", "Pakistan", "Palestine", "Afghanistan",
        "Arabian Peninsula", "Middle East", "Arab", "Dubai"}
EASIA = {"China", "Japan", "South Korea", "Korea", "Taiwan", "Republic of Korea", "Singapore", "Thailand", "Vietnam", "Malaysia", "Indonesia", "Philippines", "Hong Kong"}
EURNA = {"USA", "United States", "United States of America", "Canada", "United Kingdom", "UK", "England", "Scotland", "Wales", "Ireland", "France", "Germany",
         "Netherlands", "Belgium", "Switzerland", "Austria", "Italy", "Spain", "Portugal", "Sweden", "Norway", "Denmark", "Finland", "Iceland",
         "Poland", "Czechia", "Czech Republic", "Hungary", "Greece", "Australia", "New Zealand", "Russia", "Ukraine", "Romania", "Bulgaria", "Serbia",
         "Croatia", "Slovenia", "Slovakia", "Estonia", "Latvia", "Lithuania", "Israel", "Luxembourg", "Cyprus", "Malta",
         "US", "The Netherlands", "the Netherlands", "Europe", "Scandinavia", "Faroe Islands", "Greenland", "Republic of Ireland",
         "United Kingdom and Ireland", "United Kingdom and Republic of Ireland", "Bosnia and Herzegovina", "Macedonia", "Kosovo"}


def region_of(country):
    if not country:
        return "UNKNOWN"
    c = country.strip()
    if c in MENA: return "MENA"
    if c in EASIA: return "EASIA"
    if c in EURNA: return "EUR"
    return "OTHER"


def cmd_make_labels_proxy(a):
    with open(a.out_labels, "w") as f:
        f.write("pmid\tregion\tdetail\n")
        for r in csv.DictReader(open(a.proxy), delimiter="\t"):
            f.write(f"{r['pmid']}\t{region_of(r['first_country'])}\t{r['first_country']}\n")


def cmd_make_labels_llm(a):
    """Majority region over patient_origins (weighted by n_patients when given); recruitment country as fallback."""
    with open(a.out_labels, "w") as f:
        f.write("pmid\tregion\tdetail\n")
        for line in open(a.llm):
            d = json.loads(line)
            if not d.get("ok"):
                f.write(f"{d['pmid']}\tUNKNOWN\tllm_failed\n"); continue
            r = d["result"] or {}
            w = collections.Counter()
            for o in r.get("patient_origins", []) or []:
                reg = region_of(o.get("country") or "")
                if reg == "UNKNOWN":
                    continue
                n = o.get("n_patients") or 1
                w[reg] += n
            src = "origin"
            if not w:
                for o in r.get("recruitment_countries", []) or []:
                    reg = region_of(o.get("country") or "")
                    if reg != "UNKNOWN":
                        w[reg] += 1
                src = "recruitment"
            if not w:
                f.write(f"{d['pmid']}\tUNKNOWN\tnone\n"); continue
            top, n = w.most_common(1)[0]
            mixed = len(w) > 1 and n / sum(w.values()) < 0.7
            f.write(f"{d['pmid']}\t{'MIXED' if mixed else top}\t{src}:{dict(w)}\n")


def load_labels(path):
    return {r["pmid"]: r["region"] for r in csv.DictReader(open(path), delimiter="\t")}


def setup(a):
    eng = Engine(os.path.join(a.hpo_dir, "hp.obo"), os.path.join(a.hpo_dir, "genes_to_phenotype.txt"))
    cases = evaluable(load_store(a.store), eng)
    lab = load_labels(a.labels)
    for c in cases:
        c["region"] = lab.get(c["pmid"], "UNKNOWN")
    return eng, cases


def cmd_rank(a):
    eng, cases = setup(a)
    lpo = {r["id"]: int(r["rank_lpo"]) for r in csv.DictReader(open(os.path.join(a.results, "lpo_percase.csv")))}
    base = {r["id"]: r for r in csv.DictReader(open(os.path.join(a.results, "baseline_percase.csv"))) if r["ic_type"] == "extrinsic"}
    rows = []
    for reg in ("EUR", "EASIA", "MENA", "OTHER", "MIXED", "UNKNOWN"):
        sub = [c for c in cases if c["region"] == reg and c["id"] in base]
        if not sub:
            continue
        rf = [int(base[c["id"]]["rank"]) for c in sub]; rl = [lpo.get(c["id"], int(base[c["id"]]["rank"])) for c in sub]
        mf, ml = metrics(rf), metrics(rl)
        papers = len({c["pmid"] for c in sub})
        rows.append({"region": reg, "n": len(sub), "papers": papers, "terms": round(np.mean([len(c["terms"]) for c in sub]), 2),
                     "mean_ic": round(np.mean([float(base[c["id"]]["mean_ic"]) for c in sub]), 3),
                     "hits1_full": round(mf["hits1"], 2), "hits1_lpo": round(ml["hits1"], 2), "mrr_full": round(mf["mrr"], 4), "mrr_lpo": round(ml["mrr"], 4)})
        print(f"  {reg:8s} n={len(sub):5d} papers={papers:4d} terms={rows[-1]['terms']:5.2f} IC={rows[-1]['mean_ic']:.3f} Hits@1 full {mf['hits1']:.2f} LPO {ml['hits1']:.2f}  MRR {mf['mrr']:.4f} / {ml['mrr']:.4f}")
    with open(os.path.join(a.results, f"origin_rank_{a.tag}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


def profile(eng, cs):
    v = np.zeros(len(eng.tid))
    for c in cs:
        anc = set()
        for t in c["terms"]:
            anc |= eng.anc_sets[t]
        for t in anc:
            v[eng.tid[t]] += 1
    return v / max(1, len(cs))


def jsd(p, q):
    p = p + 1e-9; q = q + 1e-9
    p = p / p.sum(); q = q / q.sum(); m = (p + q) / 2
    return 0.5 * np.sum(p * np.log(p / m)) + 0.5 * np.sum(q * np.log(q / m))


def cmd_profile(a):
    eng, cases = setup(a)
    bygene = collections.defaultdict(list)
    for c in cases:
        if c["region"] in ("MENA", "EUR"):
            bygene[c["gene"]].append(c)
    rows = []
    for g, cs in sorted(bygene.items()):
        m = [c for c in cs if c["region"] == "MENA"]; e = [c for c in cs if c["region"] == "EUR"]
        pm = {c["pmid"] for c in m}; pe = {c["pmid"] for c in e}
        if len(m) < a.min_cases or len(e) < a.min_cases or len(pm) < 2 or len(pe) < 2:
            continue
        obs = jsd(profile(eng, m), profile(eng, e))
        # paper-level permutation: shuffle region labels across papers, keeping paper sizes
        papers = collections.defaultdict(list)
        for c in cs:
            papers[c["pmid"]].append(c)
        plist = list(papers)
        null = []
        for _ in range(a.perms):
            random.shuffle(plist)
            k = len(pm); pm_perm = set(plist[:k])
            mm = [c for p in pm_perm for c in papers[p]]; ee = [c for p in plist[k:] for c in papers[p]]
            if not mm or not ee:
                continue
            null.append(jsd(profile(eng, mm), profile(eng, ee)))
        pval = (1 + sum(1 for x in null if x >= obs)) / (1 + len(null))
        # also the depth difference and the most divergent terms (direct terms only)
        def freq(cs2):
            f = collections.Counter()
            for c in cs2:
                for t in set(c["terms"]):
                    f[t] += 1
            return {t: n / len(cs2) for t, n in f.items()}
        fm, fe = freq(m), freq(e)
        diff = sorted(set(fm) | set(fe), key=lambda t: -abs(fm.get(t, 0) - fe.get(t, 0)))[:3]
        top = "; ".join(f"{t} {fm.get(t,0):.2f} vs {fe.get(t,0):.2f}" for t in diff)
        rows.append({"gene": g, "n_mena": len(m), "papers_mena": len(pm), "n_eur": len(e), "papers_eur": len(pe),
                     "terms_mena": round(np.mean([len(c["terms"]) for c in m]), 2), "terms_eur": round(np.mean([len(c["terms"]) for c in e]), 2),
                     "jsd": round(obs, 4), "jsd_null_median": round(float(np.median(null)), 4) if null else "", "p_paper_perm": round(pval, 4), "top_diff_terms": top})
        print(f"  {g:10s} MENA {len(m):3d}/{len(pm):2d}p EUR {len(e):3d}/{len(pe):2d}p terms {rows[-1]['terms_mena']:.1f}/{rows[-1]['terms_eur']:.1f} JSD {obs:.3f} null {rows[-1]['jsd_null_median']} p={pval:.3f}")
    if rows:
        with open(os.path.join(a.results, f"origin_profile_{a.tag}.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"  genes tested: {len(rows)}; p<0.05 at paper level: {sum(1 for r in rows if r['p_paper_perm'] < 0.05)}")


def rank_with_profile(eng, case_terms, gene, gene_terms_override):
    """Rank gene among all genes when its profile is replaced by gene_terms_override (direct terms)."""
    scores = eng.bma_all_genes(case_terms)
    from run_analyses import gene_score_sets
    scores[eng.gid[gene]] = gene_score_sets(eng, case_terms, gene_terms_override)
    order = np.lexsort((np.array(eng.genes), -scores))
    return int(np.where(order == eng.gid[gene])[0][0] + 1)


def cmd_popref(a):
    """For each gene: held-out MENA case (leave-one-paper-out) ranked with the gene's profile replaced by
    the union of direct terms of (a) other MENA papers' cases, (b) a size-matched sample of EUR cases, (c) HPOA (unchanged)."""
    eng, cases = setup(a)
    bygene = collections.defaultdict(list)
    for c in cases:
        if c["region"] in ("MENA", "EUR"):
            bygene[c["gene"]].append(c)
    rows = []
    for g, cs in sorted(bygene.items()):
        m = [c for c in cs if c["region"] == "MENA"]; e = [c for c in cs if c["region"] == "EUR"]
        pm = {c["pmid"] for c in m}
        if len(m) < a.min_cases or len(e) < a.min_cases or len(pm) < 2:
            continue
        for c in m:
            train_m = [x for x in m if x["pmid"] != c["pmid"]]
            if not train_m:
                continue
            train_e = random.sample(e, min(len(e), len(train_m)))
            prof_m = sorted({t for x in train_m for t in x["terms"]})
            prof_e = sorted({t for x in train_e for t in x["terms"]})
            r_m = rank_with_profile(eng, c["terms"], g, prof_m)
            r_e = rank_with_profile(eng, c["terms"], g, prof_e)
            r_h, _, _ = eng.rank_gene(c["terms"], g)
            rows.append({"gene": g, "id": c["id"], "pmid": c["pmid"], "n_train_mena": len(train_m), "n_train_eur": len(train_e),
                         "rank_mena_ref": r_m, "rank_eur_ref": r_e, "rank_hpoa": r_h})
        # symmetric: held-out EUR cases with EUR-ref vs MENA-ref
        pe = {c["pmid"] for c in e}
        for c in random.sample(e, min(len(e), len(m))):
            train_e = [x for x in e if x["pmid"] != c["pmid"]]
            if not train_e:
                continue
            train_m = random.sample(m, min(len(m), len(train_e)))
            prof_e = sorted({t for x in train_e for t in x["terms"]}); prof_m = sorted({t for x in train_m for t in x["terms"]})
            rows.append({"gene": g, "id": c["id"], "pmid": c["pmid"], "n_train_mena": len(train_m), "n_train_eur": len(train_e),
                         "rank_mena_ref": rank_with_profile(eng, c["terms"], g, prof_m), "rank_eur_ref": rank_with_profile(eng, c["terms"], g, prof_e),
                         "rank_hpoa": eng.rank_gene(c["terms"], g)[0], "heldout_region": "EUR"})
    for r in rows:
        r.setdefault("heldout_region", "MENA")
    if not rows:
        print("  no gene meets the thresholds"); return
    with open(os.path.join(a.results, f"origin_popref_{a.tag}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "id", "pmid", "heldout_region", "n_train_mena", "n_train_eur", "rank_mena_ref", "rank_eur_ref", "rank_hpoa"]); w.writeheader(); w.writerows(rows)
    for reg in ("MENA", "EUR"):
        sub = [r for r in rows if r["heldout_region"] == reg]
        if not sub:
            continue
        print(f"  held-out {reg} cases n={len(sub)} genes={len({r['gene'] for r in sub})}: Hits@1 MENA-ref {100*np.mean([r['rank_mena_ref']==1 for r in sub]):.1f}  EUR-ref {100*np.mean([r['rank_eur_ref']==1 for r in sub]):.1f}  HPOA {100*np.mean([r['rank_hpoa']==1 for r in sub]):.1f}; "
              f"MRR {np.mean([1/r['rank_mena_ref'] for r in sub]):.3f} / {np.mean([1/r['rank_eur_ref'] for r in sub]):.3f} / {np.mean([1/r['rank_hpoa'] for r in sub]):.3f}; "
              f"same-region ref better in {100*np.mean([(r['rank_mena_ref']<r['rank_eur_ref']) if reg=='MENA' else (r['rank_eur_ref']<r['rank_mena_ref']) for r in sub]):.1f}% of cases, worse in "
              f"{100*np.mean([(r['rank_mena_ref']>r['rank_eur_ref']) if reg=='MENA' else (r['rank_eur_ref']>r['rank_mena_ref']) for r in sub]):.1f}%")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["make-labels-proxy", "make-labels-llm", "rank", "profile", "popref"])
    ap.add_argument("--proxy"); ap.add_argument("--llm"); ap.add_argument("--out-labels")
    ap.add_argument("--labels"); ap.add_argument("--hpo-dir"); ap.add_argument("--store"); ap.add_argument("--results", default="results")
    ap.add_argument("--tag", default="proxy"); ap.add_argument("--min-cases", type=int, default=5); ap.add_argument("--perms", type=int, default=500)
    a = ap.parse_args()
    {"make-labels-proxy": cmd_make_labels_proxy, "make-labels-llm": cmd_make_labels_llm, "rank": cmd_rank, "profile": cmd_profile, "popref": cmd_popref}[a.cmd](a)


if __name__ == "__main__":
    main()
