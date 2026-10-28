#!/usr/bin/env python3
"""Joint case-to-case model over all corpora (PAVS clinical and curated literature, Phenopacket
Store, PMC-Patients), adjusting for provenance (who recorded the phenotypes, and how).

Commands
  build-pavs --pavs PAVS_phenopackets.json --country pavs_pmid_country.tsv --hgnc H --out cases_pavs.tsv
      PAVS sources A, B, F, P (Saudi clinical; corpus pavs_clin), M (Saudi literature curated by
      PAVS; pavs_lit), Q (US clinical, Ziats 2020; pavs_us). D (DDG2P panel) is excluded.
      Patient region MENA for A/B/F/P/M, EUR for Q. Author region: clinical sources by site;
      M by first-author country.
  authors --cases T.tsv [T.tsv ...] --have A.tsv --out A_all.tsv
      PubMed author lists for every PMID in the case tables (esummary), merged with --have.
  pairs --cases ... --hpo-dir D --cache S.npy --authors A.tsv --r1 MENA --r2 EUR --perms N --out PREFIX
      Within-gene pairs of r1/r2 cases from different, author-disjoint publications. OLS of
      hubness-normalised similarity on: gene fixed effects; depth (log n_i + log n_j, product,
      absolute difference); provenance-pair fixed effects (unordered corpus pair); author-region
      terms (number of r1-author cases in the pair, same author region, any unknown author region);
      patient-region terms: number of r1 cases in the pair (r1 cases more/less similar to everyone)
      and same patient region, the latter split by provenance class (literature-literature,
      clinical-literature, clinical-clinical). Inference by permuting patient-region labels between
      publications within (gene, corpus, author region), so that the same-region coefficients are
      tested beyond author region and provenance; and within (gene, corpus) for the combined effect.
"""
import argparse, collections, csv, json, os, re, sys, time, urllib.parse, urllib.request
import numpy as np

CLIN = {"pavs_clin", "pavs_us"}


def cmd_build_pavs(a):
    from run_analyses import load_pavs
    from pmcp_to_cases import load_hgnc
    from origin_analyses import region_of
    hg = load_hgnc(a.hgnc)
    country = {r["pmid"]: r["first_country"] for r in csv.DictReader(open(a.country), delimiter="\t")}
    n = collections.Counter()
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "pmid", "gene", "region", "region_source", "author_region", "corpus", "source", "terms"])
        for c in load_pavs(a.pavs):
            g = hg.get((c["gene"] or "").upper().strip())
            if not g or not c["hpos"]: n["no_gene_or_terms"] += 1; continue
            s = c["source"]
            corpus = "pavs_us" if s == "Q" else "pavs_lit" if s == "M" else "pavs_clin"
            region = "EUR" if s == "Q" else "MENA"
            if s == "M":
                author = region_of(country.get(c["pmid"] or "", "")) if c["pmid"] else "UNKNOWN"
            else:
                author = region
            pmid = c["pmid"] or f"PAVSM_{c['id']}"
            w.writerow([c["id"], pmid, g, region, "pavs", author, corpus, s, ";".join(dict.fromkeys(c["hpos"]))])
            n[corpus] += 1
    print(dict(n))


def cmd_authors(a):
    have = {}
    if a.have:
        for l in open(a.have):
            p, _, au = l.rstrip("\n").partition("\t"); have[p] = au
    pm = set()
    for path in a.cases:
        pm |= {r["pmid"] for r in csv.DictReader(open(path), delimiter="\t") if r["pmid"].isdigit()}
    todo = sorted(pm - set(have))
    print(f"{len(pm)} pmids, {len(todo)} to fetch", flush=True)
    for i in range(0, len(todo), 200):
        b = todo[i:i + 200]
        data = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(b), "retmode": "json"}).encode()
        j = None
        for t in range(4):
            try:
                j = json.load(urllib.request.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", data, timeout=60)); break
            except Exception:
                time.sleep(3 * (t + 1))
        for p in b:
            r = (j or {}).get("result", {}).get(p) or {}
            have[p] = ";".join(x["name"] for x in r.get("authors", []) if x.get("authtype") == "Author")
        time.sleep(0.4)
    with open(a.out, "w") as f:
        for p in sorted(have): f.write(f"{p}\t{have[p]}\n")
    print(f"wrote {len(have)} ({sum(bool(v) for v in have.values())} with authors)")


def cmd_pairs(a):
    import case_pairs as cp
    rng = np.random.default_rng(13)
    sim, cases = cp.setup(a)
    N = len(cases)
    reg = np.array([c["region"] for c in cases]); corp = np.array([c["corpus"] or "?" for c in cases])
    aut = np.array([c["author"] for c in cases]); pm = np.array([c["pmid"] for c in cases])
    gene = np.array([c["gene"] for c in cases])
    use = np.where(np.isin(reg, [a.r1, a.r2]))[0]
    byg = collections.defaultdict(list)
    for i in use: byg[gene[i]].append(i)
    I, J = [], []
    for g, idx in byg.items():
        if len({reg[i] for i in idx}) < 2: continue
        for k, i in enumerate(idx):
            for j in idx[k + 1:]:
                if not cp.related(pm[i], pm[j]): I.append(i); J.append(j)
    I = np.array(I); J = np.array(J)
    y = sim.Z[I, J].astype(float) if a.sim == "z" else sim.S[I, J].astype(float)
    G = np.unique(gene[I], return_inverse=True)[1]
    print(f"pairs {len(I)} in {G.max() + 1} genes", flush=True)
    lit = lambda c: c not in CLIN
    pcls = np.array(["lit-lit" if lit(corp[i]) and lit(corp[j]) else "clin-clin" if not lit(corp[i]) and not lit(corp[j]) else "clin-lit"
                     for i, j in zip(I, J)])
    ctab = collections.Counter(zip(pcls, [("same" if reg[i] == reg[j] else "cross") for i, j in zip(I, J)]))
    print("pairs by provenance class x patient region:", dict(ctab))
    # fixed covariates
    L = np.log(np.maximum(sim.lens, 1))
    li, lj = L[I], L[J]
    cols = {"depth_sum": li + lj, "depth_prod": li * lj, "depth_absdiff": np.abs(li - lj)}
    ppair = np.array(["|".join(sorted((corp[i], corp[j]))) for i, j in zip(I, J)])
    for v in sorted(set(ppair))[1:]: cols[f"prov[{v}]"] = (ppair == v).astype(float)
    ar1 = (aut == a.r1).astype(float); unk = np.isin(aut, ["UNKNOWN", "MIXED", "OTHER"])
    cols["author_n_r1"] = ar1[I] + ar1[J]
    cols["author_same"] = ((aut[I] == aut[J]) & ~unk[I] & ~unk[J]).astype(float)
    cols["author_unknown"] = (unk[I] | unk[J]).astype(float)
    classes = [c for c in ("lit-lit", "clin-lit", "clin-clin") if (pcls == c).any()]

    def demean(M):
        M = np.atleast_2d(M.T).T
        s = np.zeros((G.max() + 1, M.shape[1])); np.add.at(s, G, M)
        n = np.bincount(G).astype(float)
        return M - (s / n[:, None])[G]
    F = demean(np.column_stack(list(cols.values())))
    yd = demean(y)[:, 0]
    # FWL: residualise on F once
    Fq, _ = np.linalg.qr(F)
    def resid(M): return M - Fq @ (Fq.T @ M)
    yr = resid(yd[:, None])[:, 0]
    lnames = [f"patient_n_{a.r1}"] + [f"patient_same[{c}]" for c in classes] + ["patient_same[all]"]

    def label_coefs(R):
        r1 = (R == a.r1).astype(float); same = (R[I] == R[J]).astype(float)
        X = np.column_stack([r1[I] + r1[J]] + [same * (pcls == c) for c in classes])
        Xr = resid(demean(X)); b = np.linalg.lstsq(Xr, yr, rcond=None)[0]
        X2 = np.column_stack([r1[I] + r1[J], same]); X2r = resid(demean(X2)); b2 = np.linalg.lstsq(X2r, yr, rcond=None)[0]
        return np.concatenate([b, b2[1:]])
    b0 = label_coefs(reg)
    out = []
    for null in ("gene+corpus", "gene+corpus+author"):
        strata = collections.defaultdict(lambda: collections.defaultdict(list))
        for i in use:
            key = (gene[i], corp[i]) if null == "gene+corpus" else (gene[i], corp[i], aut[i])
            strata[key][pm[i]].append(i)
        plan = [(list(p.values()), [reg[v[0]] for v in p.values()]) for p in strata.values() if len(p) > 1]
        plan = [(ix, lb) for ix, lb in plan if len(set(lb)) > 1]
        B = []
        for b in range(a.perms):
            R = reg.copy()
            for ix, lb in plan:
                lb = list(lb); rng.shuffle(lb)
                for v, l in zip(ix, lb): R[v] = l
            B.append(label_coefs(R))
        B = np.array(B); mu = B.mean(0)
        for k, nm in enumerate(lnames):
            p = (1 + np.sum(np.abs(B[:, k] - mu[k]) >= abs(b0[k] - mu[k]))) / (1 + len(B))
            out.append({"null": null, "term": nm, "coef": b0[k], "null_mean": mu[k], "excess": b0[k] - mu[k],
                        "lo": np.percentile(B[:, k], 2.5), "hi": np.percentile(B[:, k], 97.5), "p": p, "units_permuted": len(plan)})
    print(f"{'null':20s} {'term':28s} {'coef':>8s} {'nullmean':>9s} {'excess':>8s}  95% null            p")
    for r in out:
        print(f"{r['null']:20s} {r['term']:28s} {r['coef']:+8.3f} {r['null_mean']:+9.3f} {r['excess']:+8.3f}  [{r['lo']:+.3f}, {r['hi']:+.3f}]  {r['p']:.4f}")
    # fixed-effect estimates for context (author and provenance terms), from the full OLS
    X = np.column_stack([F, resid(demean(np.column_stack([(reg == a.r1).astype(float)[I] + (reg == a.r1).astype(float)[J],
                                                         (reg[I] == reg[J]).astype(float)])))])
    bf = np.linalg.lstsq(np.column_stack([F, demean(np.column_stack([(reg == a.r1).astype(float)[I] + (reg == a.r1).astype(float)[J], (reg[I] == reg[J]).astype(float)]))]), yd, rcond=None)[0]
    print("full OLS (gene FE) coefficients:")
    for nm, v in zip(list(cols) + [f"patient_n_{a.r1}", "patient_same"], bf): print(f"   {nm:40s} {v:+.3f}")
    with open(a.out + "_merged_pairs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    # provenance calibration: the same publication recorded by two processes
    if a.calibrate:
        byp = collections.defaultdict(list)
        for i in range(N): byp[pm[i], gene[i]].append(i)
        same_pub = collections.defaultdict(list)
        for (p, g), idx in byp.items():
            for k, i in enumerate(idx):
                for j in idx[k + 1:]:
                    if corp[i] != corp[j]: same_pub["|".join(sorted((corp[i], corp[j])))].append(sim.Z[i, j])
        within = collections.defaultdict(list)
        for k in range(len(I)):
            if corp[I[k]] == corp[J[k]]: within[corp[I[k]]].append(y[k])
        print("provenance calibration (z-similarity):")
        for k, v in sorted(same_pub.items()): print(f"   same publication+gene, {k:22s} pairs {len(v):6d} mean {np.mean(v):+.3f}")
        for k, v in sorted(within.items()): print(f"   different publication, same gene, both {k:12s} pairs {len(v):6d} mean {np.mean(v):+.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build-pavs", "authors", "pairs"])
    ap.add_argument("--pavs"); ap.add_argument("--country"); ap.add_argument("--hgnc"); ap.add_argument("--out", required=True)
    ap.add_argument("--cases", action="append"); ap.add_argument("--have"); ap.add_argument("--hpo-dir")
    ap.add_argument("--cache"); ap.add_argument("--authors"); ap.add_argument("--sim", choices=["z", "raw"], default="z")
    ap.add_argument("--r1", default="MENA"); ap.add_argument("--r2", default="EUR"); ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--calibrate", action="store_true")
    a = ap.parse_args()
    {"build-pavs": cmd_build_pavs, "authors": cmd_authors, "pairs": cmd_pairs}[a.cmd](a)


if __name__ == "__main__":
    main()
