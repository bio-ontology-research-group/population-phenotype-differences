#!/usr/bin/env python3
"""Which phenotypes and genes carry population differences between case reports?

Commands
  termdb  --obo hp.obo --out DIR
      FUNC termdb tables (term.txt, term2term.txt, graph_path.txt) from an OBO file,
      using the HPO accession as the internal id (the OntoFUNC convention).
  terms   --cases T.tsv --obo hp.obo --r1 MENA --r2 EUR [--strata gene|gene+author]
          [--min-cases 10] [--perms 2000] [--func-seeds 5] --out PREFIX
      Per HPO term (ancestor-propagated, phenotypic abnormality subtree): Mantel-Haenszel
      weighted within-gene difference in the proportion of r1 vs r2 cases showing the
      term. Null: region labels permuted between publications within stratum; FWER from
      the permutation distribution of max |z|. Also writes FUNC binomial input files
      (entity = HPO term, counts = most-specific annotations in r1 and r2 cases after
      balancing r1 and r2 case numbers within each gene).
  genes   --cases T.tsv --pairs P.csv --g2p genes_to_phenotype.txt --r1 MENA --r2 EUR
          [--strata gene|gene+author] [--perms 2000] --out PREFIX
      Per gene: depth-adjusted same-region minus cross-region pair similarity (from
      case_pairs.py pairs output), standardised against the publication-level permutation
      null; heterogeneity test (sum of squared per-gene z vs its permutation
      distribution); FUNC Wilcoxon input (gene -> HPOA phenotype and inheritance terms,
      value = per-gene z); AR-only vs AD-only comparison.
"""
import argparse, collections, csv, os
import numpy as np
from rank_engine import parse_obo

ROOT_PA = "HP:0000118"


def names(obo):
    out, cur = {}, None
    for l in open(obo):
        if l.startswith("id: "): cur = l[4:].strip()
        elif l.startswith("name: ") and cur: out[cur] = l[6:].strip()
        elif l.startswith("[Typedef]"): cur = None
    return out


def ancestors_fn(terms, graph_path=None):
    if graph_path:  # inferred closure from a FUNC/OntoFUNC termdb (e.g. built with ELK)
        anc = collections.defaultdict(set)
        for l in open(graph_path):
            f = l.split("\t")
            if f[1].startswith("HP:") and f[2].startswith("HP:"): anc[f[2]].add(f[1])
        for t in terms: anc[t].add(t)
        fz = {t: frozenset(v) for t, v in anc.items()}
        return lambda t: fz.get(t, frozenset([t]))
    memo = {}
    def anc(t):
        if t not in memo:
            s = {t}
            for p in terms.get(t, {}).get("is_a", []): s |= anc(p)
            memo[t] = frozenset(s)
        return memo[t]
    return anc


def read_cases(path, terms, alt, r1, r2):
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["region"] not in (r1, r2): continue
        ts = [alt.get(t, t) for t in r["terms"].split(";")]
        ts = [t for t in ts if t in terms and not terms[t]["obsolete"]]
        if ts: out.append({"id": r["id"], "pmid": r["pmid"], "gene": r["gene"], "region": r["region"],
                           "author": r.get("author_region") or "UNKNOWN", "terms": ts})
    return out


def strata_of(cases, how):
    return [(c["gene"],) if how == "gene" else (c["gene"], c["author"]) for c in cases]


def permuter(cases, strata, r1, rng):
    """Returns observed label vector and a function drawing a publication-level permutation within strata."""
    L0 = np.array([c["region"] == r1 for c in cases], int)
    groups = collections.defaultdict(lambda: collections.defaultdict(list))
    for i, (c, s) in enumerate(zip(cases, strata)): groups[s][c["pmid"]].append(i)
    plan = []
    for s, pubs in groups.items():
        idx = list(pubs.values())
        if len(idx) > 1: plan.append((idx, [L0[ii[0]] for ii in idx]))
    def draw():
        L = L0.copy()
        for idx, labs in plan:
            labs = list(labs); rng.shuffle(labs)
            for ii, l in zip(idx, labs): L[ii] = l
        return L
    return L0, draw


# ---------------------------------------------------------------- termdb
def cmd_termdb(a):
    terms, _ = parse_obo(a.obo); nm = names(a.obo)
    live = {t for t, v in terms.items() if not v["obsolete"]}
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "term.txt"), "w") as f:
        for t in sorted(live): f.write(f"{t}\t{nm.get(t, t)}\tHP:0000001\t{t}\tend\n")
    n = 0
    with open(os.path.join(a.out, "term2term.txt"), "w") as f:
        for t in sorted(live):
            for p in terms[t]["is_a"]:
                if p in live: f.write(f"{n}\tis-a\t{p}\t{t}\n"); n += 1
    n = 0
    with open(os.path.join(a.out, "graph_path.txt"), "w") as f:
        for t in sorted(live):  # BFS upwards for minimal distance
            dist = {t: 0}; frontier = [t]
            while frontier:
                nxt = []
                for x in frontier:
                    for p in terms[x]["is_a"]:
                        if p in live and p not in dist: dist[p] = dist[x] + 1; nxt.append(p)
                frontier = nxt
            for anc, d in dist.items(): f.write(f"{n}\t{anc}\t{t}\tis-a\t{d}\t0\n"); n += 1
    print(f"termdb: {len(live)} terms, {n} paths -> {a.out}")


# ---------------------------------------------------------------- terms
def cmd_terms(a):
    rng = np.random.default_rng(3)
    terms, alt = parse_obo(a.obo); nm = names(a.obo)
    anc = ancestors_fn(terms, os.path.join(a.termdb, "graph_path.txt") if a.termdb else None)
    cases = read_cases(a.cases, terms, alt, a.r1, a.r2)
    regs = collections.defaultdict(set)
    for c in cases: regs[c["gene"]].add(c["region"])
    cases = [c for c in cases if len(regs[c["gene"]]) == 2]
    strata = strata_of(cases, a.strata)
    sreg = collections.defaultdict(set)
    for c, s in zip(cases, strata): sreg[s].add(c["region"])
    keep = [i for i, s in enumerate(strata) if len(sreg[s]) == 2]
    cases = [cases[i] for i in keep]; strata = [strata[i] for i in keep]
    print(f"{len(cases)} cases ({sum(c['region'] == a.r1 for c in cases)} {a.r1}) in {len(set(strata))} strata, "
          f"{len({c['gene'] for c in cases})} genes, {len({c['pmid'] for c in cases})} publications")
    # propagated presence matrix, phenotypic-abnormality subtree, excluding the root
    prop = [{t for x in c["terms"] for t in anc(x)} for c in cases]
    cnt = collections.Counter(t for p in prop for t in p)
    T = sorted(t for t, n in cnt.items() if n >= a.min_cases and t != ROOT_PA and ROOT_PA in anc(t))
    col = {t: j for j, t in enumerate(T)}
    X = np.zeros((len(cases), len(T)), np.float32)
    for i, p in enumerate(prop):
        for t in p:
            if t in col: X[i, col[t]] = 1
    depth = np.log([len(c["terms"]) for c in cases])
    sid = {s: k for k, s in enumerate(sorted(set(strata)))}; S = np.array([sid[s] for s in strata]); ns = len(sid)
    L0, draw = permuter(cases, strata, a.r1, rng)

    def coef(L):
        n1 = np.bincount(S, weights=L, minlength=ns); n = np.bincount(S, minlength=ns); n2 = n - n1
        ok = (n1 > 0) & (n2 > 0); w = np.where(ok, n1 * n2 / np.maximum(n, 1), 0)
        a_ = np.where(L == 1, w[S] / np.maximum(n1[S], 1), -w[S] / np.maximum(n2[S], 1))
        return a_ / w.sum()

    def stat(L):
        c_ = coef(L); return c_ @ X, c_ @ depth
    t0, d0 = stat(L0)
    null = np.empty((a.perms, len(T)), np.float32); dnull = np.empty(a.perms)
    for b in range(a.perms):
        null[b], dnull[b] = stat(draw())
    mu, sd = null.mean(0), null.std(0) + 1e-9
    z = (t0 - mu) / sd
    p = (1 + (np.abs(null - mu) >= np.abs(t0 - mu)).sum(0)) / (1 + a.perms)
    maxz = np.abs((null - mu) / sd).max(1)
    fwer = (1 + (maxz[None, :] >= np.abs(z)[:, None]).sum(1)) / (1 + a.perms)
    # BH FDR
    o = np.argsort(p); q = np.empty_like(p); m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    print(f"depth (log #terms) {a.r1}-{a.r2} within strata: {d0:+.3f}; null 95% [{np.percentile(dnull,2.5):+.3f}, {np.percentile(dnull,97.5):+.3f}]")
    print(f"terms tested {len(T)}; p<0.05 {int((p < 0.05).sum())} (expected {0.05*len(T):.0f}); FDR<0.1 {int((q < 0.1).sum())}; FWER<0.05 {int((fwer < 0.05).sum())}")
    # global test: sum z^2 against its permutation distribution
    zn = (null - mu) / sd; g0 = (z ** 2).sum(); gn = (zn ** 2).sum(1)
    print(f"global (sum z^2) {g0:.0f} vs null mean {gn.mean():.0f}; p={(1 + (gn >= g0).sum()) / (1 + a.perms):.4f}")
    p1 = np.array([X[L0 == 1, j].mean() for j in range(len(T))]); p2 = np.array([X[L0 == 0, j].mean() for j in range(len(T))])
    rows = sorted(zip(T, t0, z, p, q, fwer, p1, p2, X.sum(0)), key=lambda r: r[3])
    with open(a.out + "_terms.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["hpo", "name", "mh_diff", "z", "p", "q_bh", "fwer", f"prev_{a.r1}", f"prev_{a.r2}", "n_cases"])
        for r in rows: w.writerow([r[0], nm.get(r[0], ""), f"{r[1]:.4f}", f"{r[2]:.3f}", f"{r[3]:.4g}", f"{r[4]:.4g}", f"{r[5]:.4g}", f"{r[6]:.3f}", f"{r[7]:.3f}", int(r[8])])
    for r in rows[:15]:
        print(f"  {r[0]} {nm.get(r[0],'')[:45]:45s} diff {r[1]:+.3f} z {r[2]:+.2f} p {r[3]:.4f} q {r[4]:.3f} fwer {r[5]:.3f}  prev {r[6]:.2f}/{r[7]:.2f}")
    # FUNC binomial input: gene-balanced subsamples, most-specific annotations
    def most_specific(ts):
        s = set(ts); return [t for t in s if not any(t != u and t in anc(u) for u in s)]
    byg = collections.defaultdict(lambda: {a.r1: [], a.r2: []})
    for c in cases: byg[c["gene"]][c["region"]].append(c)
    for seed in range(a.func_seeds):
        r = np.random.default_rng(100 + seed); k = collections.Counter()
        for g, d in byg.items():
            m = min(len(d[a.r1]), len(d[a.r2]))
            for reg in (a.r1, a.r2):
                for i in r.choice(len(d[reg]), m, replace=False):
                    for t in most_specific(d[reg][i]["terms"]):
                        if ROOT_PA in anc(t) and t != ROOT_PA: k[t, reg] += 1
        with open(f"{a.out}_func_binom_seed{seed}.tsv", "w") as f:
            for t in sorted({t for t, _ in k}): f.write(f"{t}\t{t}\t{k[t, a.r1]}\t{k[t, a.r2]}\n")


# ---------------------------------------------------------------- genes
def cmd_genes(a):
    rng = np.random.default_rng(5)
    cinfo = {r["id"]: r for r in csv.DictReader(open(a.cases), delimiter="\t")}
    P = list(csv.DictReader(open(a.pairs)))
    sim = np.array([float(r["sim"]) for r in P]); li = np.log([float(r["n_i"]) for r in P]); lj = np.log([float(r["n_j"]) for r in P])
    D = np.column_stack([np.ones_like(sim), li + lj, li * lj]); beta = np.linalg.lstsq(D, sim, rcond=None)[0]
    res = sim - D @ beta
    ids = sorted({r["id_i"] for r in P} | {r["id_j"] for r in P}); ix = {x: k for k, x in enumerate(ids)}
    cases = [{"id": x, "pmid": cinfo[x]["pmid"], "gene": cinfo[x]["gene"], "region": cinfo[x]["region"],
              "author": cinfo[x]["author_region"] or "UNKNOWN"} for x in ids]
    I = np.array([ix[r["id_i"]] for r in P]); J = np.array([ix[r["id_j"]] for r in P])
    genes = sorted({r["gene"] for r in P}); gx = {g: k for k, g in enumerate(genes)}; G = np.array([gx[r["gene"]] for r in P])
    L0, draw = permuter(cases, strata_of(cases, a.strata), a.r1, rng)

    def per_gene(L):
        same = (L[I] == L[J]).astype(float)
        ns = np.bincount(G, same, len(genes)); nc = np.bincount(G, 1 - same, len(genes))
        ss = np.bincount(G, same * res, len(genes)); sc = np.bincount(G, (1 - same) * res, len(genes))
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where((ns > 0) & (nc > 0), ss / ns - sc / nc, np.nan)
    d0 = per_gene(L0); null = np.array([per_gene(draw()) for _ in range(a.perms)])
    mu = np.nanmean(null, 0); sd = np.nanstd(null, 0)
    ok = np.isfinite(d0) & (sd > 1e-9) & (np.isfinite(null).mean(0) > 0.5)
    z = np.where(ok, (d0 - mu) / np.where(ok, sd, 1), np.nan)
    pg = np.array([(1 + np.sum(np.abs(null[:, k] - mu[k]) >= abs(d0[k] - mu[k]))) / (1 + np.isfinite(null[:, k]).sum()) if ok[k] else np.nan for k in range(len(genes))])
    zn = (null - mu) / np.where(sd > 0, sd, 1)
    H0 = np.nansum(z[ok] ** 2); Hn = np.nansum(zn[:, ok] ** 2, 1)
    M0 = np.nanmean(z[ok]); Mn = np.nanmean(zn[:, ok], 1)
    print(f"genes with both pair types: {int(ok.sum())}")
    print(f"mean per-gene z {M0:+.3f} (null 95% [{np.percentile(Mn,2.5):+.3f}, {np.percentile(Mn,97.5):+.3f}])")
    print(f"heterogeneity sum z^2 = {H0:.1f} vs null mean {Hn.mean():.1f}; p={(1 + (Hn >= H0).sum()) / (1 + a.perms):.4f}")
    nsame = np.bincount(G, (L0[I] == L0[J]).astype(float), len(genes)); npair = np.bincount(G, minlength=len(genes))
    # annotations
    g2p = collections.defaultdict(set)
    for r in csv.DictReader(open(a.g2p), delimiter="\t"): g2p[r["gene_symbol"]].add(r["hpo_id"])
    inh = {g: ("AR" if "HP:0000007" in s and "HP:0000006" not in s else "AD" if "HP:0000006" in s and "HP:0000007" not in s else
               "both" if "HP:0000006" in s else "other") for g, s in g2p.items()}
    order = np.argsort(np.where(ok, -np.abs(z), 0))
    with open(a.out + "_genes.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["gene", "same_minus_cross", "z", "p", "pairs", "same_pairs", "inheritance"])
        for k in order:
            if ok[k]: w.writerow([genes[k], f"{d0[k]:.4f}", f"{z[k]:.3f}", f"{pg[k]:.4g}", int(npair[k]), int(nsame[k]), inh.get(genes[k], "NA")])
    print("  top genes by |z|:")
    for k in order[:12]:
        if ok[k]: print(f"   {genes[k]:10s} d {d0[k]:+.3f} z {z[k]:+.2f} p {pg[k]:.4f} pairs {int(npair[k])} inh {inh.get(genes[k],'NA')}")
    zz = {genes[k]: z[k] for k in range(len(genes)) if ok[k]}
    from scipy.stats import mannwhitneyu
    ar = [v for g, v in zz.items() if inh.get(g) == "AR"]; ad = [v for g, v in zz.items() if inh.get(g) == "AD"]
    if ar and ad:
        u = mannwhitneyu(ar, ad)
        print(f"  AR-only genes n={len(ar)} mean z {np.mean(ar):+.3f}; AD-only n={len(ad)} mean z {np.mean(ad):+.3f}; Mann-Whitney p={u.pvalue:.3f}")
    with open(a.out + "_func_wilcoxon.tsv", "w") as f:
        for g, v in sorted(zz.items()):
            for t in sorted(g2p.get(g, ())): f.write(f"{g}\t{t}\t{v:.5f}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["termdb", "terms", "genes"])
    ap.add_argument("--obo"); ap.add_argument("--cases"); ap.add_argument("--pairs"); ap.add_argument("--g2p")
    ap.add_argument("--r1", default="MENA"); ap.add_argument("--r2", default="EUR")
    ap.add_argument("--strata", choices=["gene", "gene+author"], default="gene")
    ap.add_argument("--min-cases", type=int, default=10); ap.add_argument("--perms", type=int, default=2000)
    ap.add_argument("--func-seeds", type=int, default=5); ap.add_argument("--termdb", help="FUNC termdb whose graph_path.txt gives the (inferred) ancestor closure")
    ap.add_argument("--dummy-unused", help=argparse.SUPPRESS); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    {"termdb": cmd_termdb, "terms": cmd_terms, "genes": cmd_genes}[a.cmd](a)


if __name__ == "__main__":
    main()
