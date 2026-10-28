#!/usr/bin/env python3
"""Consequences of regional documentation for gene-disease annotations and phenotype-based variant interpretation.

rank   For every case of the six populations (Store, PMC-Patients, abstracts), rank the causative gene among all HPOA
       genes by symmetric BMA Resnik similarity between the case and each gene's HPO annotation profile
       (genes_to_phenotype.txt, extrinsic IC), as phenotype-driven variant prioritisation does. Annotations supported
       only by the case's own publication are removed for that case (leave-publication-out). Rankings use all terms,
       only bedside/history terms and only test-based terms (termclass.jsonl). The case's coverage by the annotation
       profile of its gene is the fraction of its terms that equal or lie below an annotated term of the gene that
       is more specific than an organ system. Parallel over --chunks; writes OUT_chunk<k>.csv.
stats  Population vs EUR differences in the per-case outcomes (log2 rank with all, bedside-only and test-only terms,
       top-10 indicator, coverage) with Mantel-Haenszel weights within gene x corpus strata (total) or gene x corpus x
       author region strata (patient-carried) and permutation of labels between reporting groups (as in
       population_model.py terms); and within each population the paired loss from removing the test-based
       findings (log2 rank bedside-only minus all terms), for cases with both kinds. Writes OUT_stats.csv.

  gda_impact.py rank  --cases A.tsv [...] --hpo-dir hpo --termclass T --chunks 40 --out P
  gda_impact.py stats --percase 'P_chunk*.csv' --cases A.tsv [...] --groups G --perms 2000 --out P
"""
import argparse, collections, csv, glob, json, math, os, sys
import numpy as np

from rank_engine import Engine

POPS = ["EUR", "MENA", "SASIA", "LATAM", "AFR", "EASIA"]
TEST = {"laboratory", "imaging", "physiological_test", "pathology"}
BED = {"clinical_exam", "symptom_history"}
ROOT_PA = "HP:0000118"


def build_support(hpoa, g2p):
    """(gene, hpo) -> set of references, through the disease identifiers of genes_to_phenotype.txt."""
    dref = collections.defaultdict(set)
    for line in open(hpoa):
        if line.startswith("#") or line.startswith("database_id"): continue
        c = line.rstrip("\n").split("\t")
        for r in c[4].split(";"):
            if r.strip(): dref[(c[0], c[3])].add(r.strip())
    support = collections.defaultdict(set)
    for line in open(g2p):
        if line.startswith("ncbi_gene_id") or line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        support[(c[1], c[2])] |= dref.get((c[5], c[2]), set())
    return support


def gene_score_sets(eng, P, G):
    if not G: return 0.0
    ancP = set().union(*(eng.anc_sets[p] for p in P)); ancG = set().union(*(eng.anc_sets[g] for g in G))
    s1 = np.mean([max((eng.ic[eng.tid[a]] for a in eng.anc_sets[p] & ancG), default=0.0) for p in P])
    s2 = np.mean([max((eng.ic[eng.tid[a]] for a in eng.anc_sets[g] & ancP), default=0.0) for g in G])
    return (s1 + s2) / 2


def read_cases(paths, pops=POPS):
    out = []
    for p in paths:
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r["region"] in pops: out.append(r)
    return out


def cmd_rank(a):
    eng = Engine(os.path.join(a.hpo_dir, "hp.obo"), os.path.join(a.hpo_dir, "genes_to_phenotype.txt"))
    eng.set_ic("extrinsic")
    support = build_support(os.path.join(a.hpo_dir, "phenotype.hpoa"), os.path.join(a.hpo_dir, "genes_to_phenotype.txt"))
    sole = collections.defaultdict(list)
    for (g, h), refs in support.items():
        pm = {r[5:] for r in refs if r.startswith("PMID:")}
        if len(pm) == 1 and len(refs) == 1: sole[next(iter(pm))].append((g, h))
    tc = {x["id"]: x["class"] for x in map(json.loads, open(a.termclass))}
    organ = {t for t in eng.terms if ROOT_PA in eng.terms[t]["is_a"]} | {ROOT_PA, "HP:0000001"}
    cases = [r for r in read_cases(a.cases) if r["gene"] in eng.gid]
    cases = cases[a.chunk::a.chunks]
    gnames = np.array(eng.genes)
    rows = []

    def rank(terms, affected):
        if not terms: return None
        sc = eng.bma_all_genes(terms)
        for g, hs in affected.items():
            if g in eng.gid: sc[eng.gid[g]] = gene_score_sets(eng, terms, [t for t in eng.g2t[g] if t not in hs])
        order = np.lexsort((gnames, -sc))
        return int(np.where(order == eng.gid[g0])[0][0] + 1)
    for i, r in enumerate(cases):
        terms = eng.map_terms(r["terms"].split(";"))
        if not terms: continue
        g0 = r["gene"]
        affected = collections.defaultdict(set)
        for g, h in sole.get(r["pmid"], []): affected[g].add(h)
        bed = [t for t in terms if tc.get(t) in BED]; tst = [t for t in terms if tc.get(t) in TEST]
        prof = [t for t in eng.g2t[g0] if t not in affected.get(g0, set()) and t not in organ]
        cov = np.mean([any(p in eng.anc_sets[t] for p in prof) for t in terms])
        rows.append({"id": r["id"], "gene": g0, "region": r["region"], "corpus": r["corpus"], "n_terms": len(terms),
                     "n_bed": len(bed), "n_test": len(tst), "n_removed": len(affected.get(g0, ())),
                     "rank_all": rank(terms, affected), "rank_bed": rank(bed, affected), "rank_test": rank(tst, affected),
                     "coverage": round(float(cov), 4), "n_genes": len(eng.genes)})
        if i % 500 == 0: print(f"chunk {a.chunk}: {i}/{len(cases)}", flush=True)
    with open(f"{a.out}_chunk{a.chunk}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def cmd_stats(a):
    from population_model import load_groups, group_of, make_plan, draw
    rng = np.random.default_rng(31)
    G = load_groups(a.groups)
    meta = {r["id"]: r for r in read_cases(a.cases)}
    pc = [r for f in sorted(glob.glob(a.percase)) for r in csv.DictReader(open(f))]
    for r in pc:
        m = meta[r["id"]]; r["author"] = m.get("author_region") or "UNKNOWN"; r["group"] = group_of(G, m["pmid"])
        for k in ("rank_all", "rank_bed", "rank_test"):
            r[k] = float(r[k]) if r[k] not in ("", "None") else float("nan")
        r["coverage"] = float(r["coverage"])
        for k in ("n_terms", "n_bed", "n_test", "n_removed"): r[k] = int(r[k])
    out = []
    lg = lambda x: np.log2(x)
    outcomes = {"log2 rank, all terms": lambda r: lg(r["rank_all"]), "top-10, all terms": lambda r: float(r["rank_all"] <= 10),
                "log2 rank, bedside terms only": lambda r: lg(r["rank_bed"]),
                "log2 rank, test-based terms only": lambda r: lg(r["rank_test"]),
                "coverage by gene annotation": lambda r: r["coverage"],
                "loss from removing test-based findings (log2 rank)": lambda r: lg(r["rank_bed"]) - lg(r["rank_all"])}
    # within-population paired loss from removing the test-based findings
    for P in POPS:
        sub = [r for r in pc if r["region"] == P and not math.isnan(r["rank_bed"]) and r["n_test"] > 0]
        d = np.array([lg(r["rank_bed"]) - lg(r["rank_all"]) for r in sub])
        t10a = np.mean([r["rank_all"] <= 10 for r in sub]); t10b = np.mean([r["rank_bed"] <= 10 for r in sub])
        bs = [np.mean(rng.choice(d, len(d))) for _ in range(1000)]
        out.append({"analysis": "paired", "population": P, "outcome": "log2 rank bedside-only minus all terms", "estimate": d.mean(),
                    "lo": np.percentile(bs, 2.5), "hi": np.percentile(bs, 97.5), "p": "", "n_cases": len(sub),
                    "extra": f"top-10 all {100 * t10a:.1f}% -> bedside-only {100 * t10b:.1f}%"})
        print(f"paired {P}: n={len(sub)} loss {d.mean():+.3f} log2; top-10 {100 * t10a:.1f}% -> {100 * t10b:.1f}%; "
              f"mean coverage {np.mean([r['coverage'] for r in pc if r['region'] == P]):.3f}")
    for P in POPS[1:]:
        for version in ("total", "patient-carried"):
            key = (lambda r: (r["gene"], r["corpus"])) if version == "total" else (lambda r: (r["gene"], r["corpus"], r["author"]))
            cs = [r for r in pc if r["region"] in (P, "EUR")]
            sreg = collections.defaultdict(set)
            for r in cs: sreg[key(r)].add(r["region"])
            cs = [r for r in cs if len(sreg[key(r)]) == 2]
            sid = {s: k for k, s in enumerate(sorted({key(r) for r in cs}))}; S = np.array([sid[key(r)] for r in cs]); ns = len(sid)
            lab0 = np.array([r["region"] for r in cs])
            plan = make_plan(cs, range(len(cs)), lambda i: key(cs[i]), lambda i: cs[i]["group"], lab0)
            Y = np.array([[f(r) for f in outcomes.values()] for r in cs], float)
            ok = ~np.isnan(Y)

            def stat(L):
                res = []
                for j in range(Y.shape[1]):
                    m = ok[:, j]; L1 = ((L == P) & m).astype(float); L2 = ((L != P) & m).astype(float)
                    n1 = np.bincount(S, weights=L1, minlength=ns); n2 = np.bincount(S, weights=L2, minlength=ns)
                    w = np.where((n1 > 0) & (n2 > 0), n1 * n2 / np.maximum(n1 + n2, 1), 0)
                    y = np.where(m, Y[:, j], 0)
                    c = (L1 * w[S] / np.maximum(n1[S], 1) - L2 * w[S] / np.maximum(n2[S], 1)) / max(w.sum(), 1e-12)
                    res.append(c @ y)
                return np.array(res)
            t0 = stat(lab0); null = np.array([stat(draw(plan, lab0, rng)) for _ in range(a.perms)]); mu = null.mean(0)
            p = (1 + (np.abs(null - mu) >= np.abs(t0 - mu)).sum(0)) / (1 + a.perms)
            lo, hi = np.percentile(null - mu, [2.5, 97.5], axis=0)
            for j, k in enumerate(outcomes):
                out.append({"analysis": version, "population": P, "outcome": k, "estimate": t0[j] - mu[j], "lo": lo[j], "hi": hi[j],
                            "p": p[j], "n_cases": len(cs), "extra": f"{int((lab0 == P).sum())} {P}; {ns} strata"})
                print(f"{version:15s} {P:6s} {k:52s} {t0[j] - mu[j]:+.4f}  p={p[j]:.4f}")
    with open(a.out + "_stats.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["rank", "stats"])
    ap.add_argument("--cases", action="append", required=True); ap.add_argument("--hpo-dir"); ap.add_argument("--termclass")
    ap.add_argument("--chunks", type=int, default=1); ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--percase"); ap.add_argument("--groups"); ap.add_argument("--perms", type=int, default=2000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    {"rank": cmd_rank, "stats": cmd_stats}[a.cmd](a)


if __name__ == "__main__":
    main()
