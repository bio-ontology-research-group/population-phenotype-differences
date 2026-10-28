#!/usr/bin/env python3
"""Case-to-case tests of population effects (no disease/gene reference profile involved).

Every case is compared with other cases of the same gene, never with an aggregated HPOA/OMIM
profile. Similarity: symmetric best-match-average Resnik between the observed HPO terms of two
cases, with information content estimated from the case corpus itself (ancestor-propagated term
frequency over cases), so HPOA is not used at all. Pairs of cases from the same publication are
never compared.

Tests (regions: patient origin labels; author region = first-author affiliation proxy):
  pairs     : within gene, are same-region pairs more similar than cross-region pairs?
              Statistic: mean over genes of [mean sim(same-region pairs) - mean sim(cross pairs)],
              after regressing out depth (log n_terms of both cases). Null: region labels
              permuted across publications within (gene, author region) strata.
  retrieval : nearest-case gene ranking. A query case ranks all genes by its maximum similarity
              to any library case of each gene (excluding its own publication). For the query's
              own gene the library is restricted to k cases from the same region or k cases from
              the other region (other genes: identical full library), and the two ranks are compared
              per query; inference by bootstrap over publications.

Input: a case table TSV (id, pmid, gene, region, author_region, terms[;-separated HPO IDs]).
Usage:
  case_pairs.py build-store --store DIR --labels labels_llm.tsv --author-labels labels_proxy.tsv --hpo-dir H --out cases_store.tsv
  case_pairs.py pairs|retrieval --cases cases.tsv [--cases more.tsv] --hpo-dir H --out PREFIX --r1 MENA --r2 EUR
"""
import argparse, csv, os, random, collections, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rank_engine import parse_obo

rng = np.random.default_rng(11); random.seed(11)


# ---------------------------------------------------------------- ontology / similarity
class CaseSim:
    def __init__(self, obo, cases):
        self.terms, self.alt = parse_obo(obo)
        anc = {}
        def get(t):
            if t in anc: return anc[t]
            s = {t}
            for p in self.terms[t]["is_a"]:
                if p in self.terms: s |= get(p)
            anc[t] = s; return s
        sys.setrecursionlimit(100000)
        for t in self.terms: get(t)
        self.anc = anc
        # corpus IC over cases (each case counted once per term in its ancestor closure)
        cnt = collections.Counter()
        for c in cases:
            cl = set()
            for t in c["terms"]: cl |= anc[t]
            cnt.update(cl)
        N = len(cases)
        self.used = sorted({t for c in cases for t in c["terms"]})
        self.uid = {t: i for i, t in enumerate(self.used)}
        allanc = sorted({a for t in self.used for a in anc[t]})
        aid = {a: i for i, a in enumerate(allanc)}
        ic = np.array([-np.log2(cnt[a] / N) if cnt[a] else 0.0 for a in allanc])
        B = np.zeros((len(self.used), len(allanc)), dtype=bool)
        for t in self.used:
            B[self.uid[t], [aid[a] for a in anc[t]]] = True
        # Resnik matrix over used terms: R[t1,t2] = max IC of common ancestors
        R = np.zeros((len(self.used), len(self.used)), dtype=np.float32)
        for t in self.used:
            cols = [aid[a] for a in anc[t]]
            R[self.uid[t]] = (B[:, cols] * ic[cols]).max(axis=1)
        self.R = R
        # concatenated term index for all cases
        self.cases = cases
        self.lens = np.array([len(c["terms"]) for c in cases])
        self.off = np.concatenate([[0], np.cumsum(self.lens)[:-1]])
        self.cat = np.array([self.uid[t] for c in cases for t in c["terms"]])

    def full(self, cache=None, workers=1, block=2000):
        """Full case x case similarity (float32) and its hubness-normalised version
        Z_ij = ((S_ij - mu_i)/sd_i + (S_ij - mu_j)/sd_j) / 2, mu/sd over off-diagonal entries.
        With workers > 1 the rows are computed by forked processes into a memmap at `cache`.
        Z is built in row blocks, so peak memory stays near the size of S plus Z."""
        N = len(self.cases)
        if cache and os.path.exists(cache):
            S = np.load(cache, mmap_mode="r")   # read-only map; Z below is the only full matrix held in memory
        elif workers > 1 and cache:
            import multiprocessing as mp
            tmp = cache + ".part"
            np.lib.format.open_memmap(tmp, mode="w+", dtype=np.float32, shape=(N, N)).flush()
            global _FULL_SELF, _FULL_TMP
            _FULL_SELF, _FULL_TMP = self, tmp
            starts = list(range(0, N, block))
            with mp.get_context("fork").Pool(workers) as pool:
                for k, _ in enumerate(pool.imap_unordered(_full_block, [(b, min(b + block, N)) for b in starts])):
                    print(f"  sim blocks {k + 1}/{len(starts)}", flush=True)
            os.replace(tmp, cache)
            S = np.load(cache, mmap_mode=None)
        else:
            S = np.empty((N, N), dtype=np.float32)
            for i in range(N):
                S[i] = self.sim_row(i)
                if i % 2000 == 0: print(f"  sim rows {i}/{N}", flush=True)
            if cache: np.save(cache, S)
        mu = np.empty(N); sd = np.empty(N)
        for b in range(0, N, block):
            x = S[b:b + block].astype(np.float64)
            r = np.arange(b, min(b + block, N)); x[r - b, r] = np.nan
            mu[r] = np.nanmean(x, axis=1); sd[r] = np.nanstd(x, axis=1)
        sd[sd == 0] = 1
        Z = np.empty((N, N), dtype=np.float32)
        for b in range(0, N, block):
            e = min(b + block, N); x = S[b:e].astype(np.float64)
            Z[b:e] = ((x - mu[b:e, None]) / sd[b:e, None] + (x - mu[None, :]) / sd[None, :]) / 2
        self.S, self.Z = S, Z

    def sim_row(self, i):
        """BMA Resnik between case i and every case."""
        Q = [self.uid[t] for t in self.cases[i]["terms"]]
        G = self.R[np.ix_(Q, self.cat)]                      # |Q| x total_terms
        d1 = np.maximum.reduceat(G, self.off, axis=1).mean(axis=0)   # per case: mean over q of best match in case
        d2 = np.add.reduceat(G.max(axis=0), self.off) / self.lens    # per case: mean over case terms of best match in q
        return (d1 + d2) / 2.0


def _full_block(be):
    b, e = be
    S = np.load(_FULL_TMP, mmap_mode="r+")
    for i in range(b, e): S[i] = _FULL_SELF.sim_row(i)
    S.flush(); del S
    return b


def map_terms(sim, hpos):
    out = []
    for h in hpos:
        h = sim.alt.get(h, h) if hasattr(sim, "alt") else h
        out.append(h)
    return out


# ---------------------------------------------------------------- IO
def read_cases(paths, obo_terms, alt):
    cs = []
    for p in paths:
        for r in csv.DictReader(open(p), delimiter="\t"):
            ts = []
            for t in r["terms"].split(";"):
                t = alt.get(t, t)
                if t in obo_terms and t not in ts: ts.append(t)
            if ts and r["gene"]:
                cs.append({"id": r["id"], "pmid": r["pmid"], "gene": r["gene"], "region": r["region"],
                           "author": r.get("author_region", "UNKNOWN") or "UNKNOWN", "terms": ts, "corpus": r.get("corpus", "")})
    return cs


def cmd_build_store(a):
    from run_analyses import load_store
    lab = {r["pmid"]: r["region"] for r in csv.DictReader(open(a.labels), delimiter="\t")}
    aut = {r["pmid"]: r["region"] for r in csv.DictReader(open(a.author_labels), delimiter="\t")} if a.author_labels else {}
    terms, alt = parse_obo(os.path.join(a.hpo_dir, "hp.obo"))
    n = 0
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["id", "pmid", "gene", "region", "author_region", "corpus", "terms"])
        for c in load_store(a.store):
            ts = [alt.get(t, t) for t in c["hpos"]]
            ts = [t for t in dict.fromkeys(ts) if t in terms]
            if not ts or not c["gene"] or not c["pmid"]: continue
            w.writerow([c["id"], c["pmid"], c["gene"], lab.get(c["pmid"], "UNKNOWN"), aut.get(c["pmid"], "UNKNOWN"), "store", ";".join(ts)]); n += 1
    print("wrote", n, "cases")


# ---------------------------------------------------------------- tests
AUTHORS = {}


def related(p, q):
    """Same publication, or (with --authors) publications sharing an author."""
    return p == q or bool(AUTHORS.get(p, set()) & AUTHORS.get(q, set()))


def setup(a):
    if a.authors:
        for l in open(a.authors):
            p, _, au = l.rstrip("\n").partition("\t")
            AUTHORS[p] = {x.strip().lower() for x in au.split(";") if x.strip()}
    terms, alt = parse_obo(os.path.join(a.hpo_dir, "hp.obo"))
    cases = read_cases(a.cases, terms, alt)
    print(f"cases {len(cases)}; regions {collections.Counter(c['region'] for c in cases).most_common()}", flush=True)
    sim = CaseSim(os.path.join(a.hpo_dir, "hp.obo"), cases)
    sim.full(a.cache, workers=getattr(a, "workers", 1))
    M = sim.Z if a.sim == "z" else sim.S
    sim.sim_row = lambda i: M[i]
    return sim, cases


def pair_table(sim, cases, r1, r2):
    """All cross-publication pairs within gene among cases labelled r1 or r2."""
    bygene = collections.defaultdict(list)
    for i, c in enumerate(cases):
        if c["region"] in (r1, r2): bygene[c["gene"]].append(i)
    rows = []
    for g, idx in bygene.items():
        if len(idx) < 3 or len({cases[i]["region"] for i in idx}) < 2: continue
        if len({cases[i]["pmid"] for i in idx}) < 2: continue
        for k, i in enumerate(idx):
            s = sim.sim_row(i)
            for j in idx[k + 1:]:
                if related(cases[i]["pmid"], cases[j]["pmid"]): continue
                rows.append((g, i, j, float(s[j])))
    return rows


def pair_stat(rows, cases, region, depth_coef):
    """mean over genes of (mean same-region residual sim - mean cross-region residual sim)."""
    per = collections.defaultdict(lambda: [[], []])
    for g, i, j, s in rows:
        res = s - depth_coef(i, j)
        same = region[i] == region[j]
        per[g][0 if same else 1].append(res)
    d = [np.mean(v[0]) - np.mean(v[1]) for v in per.values() if v[0] and v[1]]
    return (np.mean(d) if d else np.nan), len(d)


def cmd_pairs(a):
    sim, cases = setup(a)
    rows = pair_table(sim, cases, a.r1, a.r2)
    print(f"pairs {len(rows)} in {len({r[0] for r in rows})} genes", flush=True)
    # depth adjustment: OLS of sim on log n_i + log n_j + log n_i*log n_j (fitted once on all pairs)
    L = np.log([len(c["terms"]) for c in cases])
    X = np.array([[1, L[i] + L[j], L[i] * L[j]] for _, i, j, _ in rows]); y = np.array([r[3] for r in rows])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    dc = lambda i, j: beta[0] + beta[1] * (L[i] + L[j]) + beta[2] * L[i] * L[j]
    region = {i: c["region"] for i, c in enumerate(cases)}
    obs, ng = pair_stat(rows, cases, region, dc)
    # also: same-author-region effect, and patient effect within author strata
    author = {i: c["author"] for i, c in enumerate(cases)}
    obs_auth, _ = pair_stat(rows, cases, author, dc)
    # permutation: shuffle patient-region labels among publications within (gene, author region) strata
    pubs = collections.defaultdict(set)
    for g, i, j, _ in rows:
        pubs[(g, cases[i]["author"])].add(cases[i]["pmid"]); pubs[(g, cases[j]["author"])].add(cases[j]["pmid"])
    pub_region = {c["pmid"]: c["region"] for c in cases}
    null = []
    for _ in range(a.perms):
        newlab = {}
        for (g, au), ps in pubs.items():
            ps = sorted(ps); labs = [pub_region[p] for p in ps]; rng.shuffle(labs)
            for p, l in zip(ps, labs): newlab[(g, p)] = l
        reg = {}
        for g, i, j, _ in rows:
            reg[i] = newlab[(g, cases[i]["pmid"])]; reg[j] = newlab[(g, cases[j]["pmid"])]
        null.append(pair_stat(rows, cases, reg, dc)[0])
    null = np.array(null)
    c0 = null.mean()
    p = (1 + np.sum(np.abs(null - c0) >= abs(obs - c0))) / (1 + len(null))
    # author effect beyond patient region: permute author labels across publications within (gene, patient region)
    pubs2 = collections.defaultdict(set)
    for g, i, j, _ in rows:
        pubs2[(g, cases[i]["region"])].add(cases[i]["pmid"]); pubs2[(g, cases[j]["region"])].add(cases[j]["pmid"])
    pub_author = {c["pmid"]: c["author"] for c in cases}
    nullA = []
    for _ in range(a.perms):
        newlab = {}
        for (g, rr), ps in pubs2.items():
            ps = sorted(ps); labs = [pub_author[q] for q in ps]; rng.shuffle(labs)
            for q, l in zip(ps, labs): newlab[(g, q)] = l
        aa = {}
        for g, i, j, _ in rows:
            aa[i] = newlab[(g, cases[i]["pmid"])]; aa[j] = newlab[(g, cases[j]["pmid"])]
        nullA.append(pair_stat(rows, cases, aa, dc)[0])
    nullA = np.array(nullA); cA = nullA.mean()
    # combined (population + recorder) effect: permute patient-region labels across publications within gene only
    pubs3 = collections.defaultdict(set)
    for g, i, j, _ in rows:
        pubs3[g].add(cases[i]["pmid"]); pubs3[g].add(cases[j]["pmid"])
    nullU = []
    for _ in range(a.perms):
        newlab = {}
        for g, ps in pubs3.items():
            ps = sorted(ps); labs = [pub_region[q] for q in ps]; rng.shuffle(labs)
            for q, l in zip(ps, labs): newlab[(g, q)] = l
        rr = {}
        for g, i, j, _ in rows:
            rr[i] = newlab[(g, cases[i]["pmid"])]; rr[j] = newlab[(g, cases[j]["pmid"])]
        nullU.append(pair_stat(rows, cases, rr, dc)[0])
    nullU = np.array(nullU); cU = nullU.mean()
    pU = (1 + np.sum(np.abs(nullU - cU) >= abs(obs - cU))) / (1 + len(nullU))
    print(f"  COMBINED region effect (null: labels permuted within gene): {obs:+.4f}; null mean {cU:+.4f} "
          f"95% [{np.percentile(nullU,2.5):+.4f}, {np.percentile(nullU,97.5):+.4f}]; excess {obs-cU:+.4f}; p={pU:.4f}")
    pA = (1 + np.sum(np.abs(nullA - cA) >= abs(obs_auth - cA))) / (1 + len(nullA))
    # descriptive: raw mean sims per pair type
    typ = collections.defaultdict(list)
    for g, i, j, s in rows:
        typ["-".join(sorted([cases[i]["region"], cases[j]["region"]]))].append(s)
    print(f"  genes with both pair types: {ng}")
    for k, v in sorted(typ.items()): print(f"  {k:12s} pairs {len(v):7d} mean sim {np.mean(v):.4f}")
    print(f"  PATIENT region same-minus-cross: {obs:+.4f}; permutation null (author strata fixed) mean {c0:+.4f} "
          f"95% [{np.percentile(null,2.5):+.4f}, {np.percentile(null,97.5):+.4f}]; excess {obs-c0:+.4f}; p={p:.4f}")
    print(f"  AUTHOR region same-minus-cross: {obs_auth:+.4f}; null (patient strata fixed) mean {cA:+.4f} "
          f"95% [{np.percentile(nullA,2.5):+.4f}, {np.percentile(nullA,97.5):+.4f}]; excess {obs_auth-cA:+.4f}; p={pA:.4f}")
    with open(a.out + "_pairs.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["gene", "id_i", "id_j", "region_i", "region_j", "author_i", "author_j", "n_i", "n_j", "sim"])
        for g, i, j, s in rows:
            w.writerow([g, cases[i]["id"], cases[j]["id"], cases[i]["region"], cases[j]["region"], cases[i]["author"], cases[j]["author"],
                        len(cases[i]["terms"]), len(cases[j]["terms"]), round(s, 5)])
    with open(a.out + "_pairs_summary.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["r1", "r2", "sim", "pairs", "genes", "patient_stat", "patient_null_mean", "patient_null_lo", "patient_null_hi", "patient_p",
                                       "author_stat", "author_null_mean", "author_null_lo", "author_null_hi", "author_p", "combined_null_mean", "combined_p"])
        w.writerow([a.r1, a.r2, a.sim, len(rows), ng, round(obs, 5), round(c0, 5), round(np.percentile(null, 2.5), 5), round(np.percentile(null, 97.5), 5), round(p, 5),
                    round(obs_auth, 5), round(cA, 5), round(np.percentile(nullA, 2.5), 5), round(np.percentile(nullA, 97.5), 5), round(pA, 5), round(cU, 5), round(pU, 5)])


def cmd_retrieval(a):
    sim, cases = setup(a)
    genes = sorted({c["gene"] for c in cases})
    gidx = {g: k for k, g in enumerate(genes)}
    cg = np.array([gidx[c["gene"]] for c in cases])
    pm = np.array([c["pmid"] for c in cases])
    reg = np.array([c["region"] for c in cases])
    bygene = collections.defaultdict(list)
    for i, c in enumerate(cases): bygene[c["gene"]].append(i)
    rows = []
    for g, idx in bygene.items():
        for qi in idx:
            if reg[qi] not in (a.r1, a.r2): continue
            other = a.r2 if reg[qi] == a.r1 else a.r1
            same_lib = [j for j in idx if reg[j] == reg[qi] and not related(pm[j], pm[qi])]
            oth_lib = [j for j in idx if reg[j] == other and not related(pm[j], pm[qi])]
            k = min(len(same_lib), len(oth_lib), a.k)
            if k < 1: continue
            s = sim.sim_row(qi)
            # best score per gene over library cases (exclude query publication)
            valid = pm != pm[qi]
            if AUTHORS: valid &= np.array([not related(p, pm[qi]) for p in pm])
            best = np.full(len(genes), -1.0)
            np.maximum.at(best, cg[valid], s[valid])
            res = {}
            for name, lib in (("same", same_lib), ("other", oth_lib)):
                ranks = []
                for _ in range(a.reps):
                    sub = random.sample(lib, k)
                    b = best.copy(); b[gidx[g]] = s[sub].max()
                    ranks.append(1 + np.sum(b > b[gidx[g]]))
                res[name] = float(np.mean(ranks))
            rows.append({"id": cases[qi]["id"], "pmid": pm[qi], "gene": g, "region": reg[qi], "author": cases[qi]["author"], "k": k,
                         "rank_same": res["same"], "rank_other": res["other"]})
    if not rows:
        print("no evaluable queries"); return
    with open(a.out + "_retrieval.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    out = []
    for r0 in (a.r1, a.r2):
        R = [r for r in rows if r["region"] == r0]
        if not R: continue
        d = np.array([np.log2(r["rank_other"]) - np.log2(r["rank_same"]) for r in R])
        h_same = np.mean([r["rank_same"] == 1 for r in R]); h_oth = np.mean([r["rank_other"] == 1 for r in R])
        # bootstrap over publications
        pubs = sorted({r["pmid"] for r in R}); byp = collections.defaultdict(list)
        for r, x in zip(R, d): byp[r["pmid"]].append(x)
        bs = []
        for _ in range(2000):
            samp = rng.choice(len(pubs), len(pubs)); v = np.concatenate([byp[pubs[t]] for t in samp]); bs.append(v.mean())
        lo, hi = np.percentile(bs, [2.5, 97.5])
        print(f"  queries {r0:6s} n={len(R):5d} pubs={len(pubs):4d} genes={len({r['gene'] for r in R}):4d}  "
              f"Hits@1 same-region lib {100*h_same:5.1f} other {100*h_oth:5.1f}; mean log2(rank_other/rank_same) {d.mean():+.3f} [{lo:+.3f}, {hi:+.3f}]  "
              f"same better {100*np.mean(d>0):.1f}% worse {100*np.mean(d<0):.1f}%")
        out.append({"query_region": r0, "n": len(R), "pubs": len(pubs), "genes": len({r['gene'] for r in R}), "hits1_same": round(100*h_same, 2),
                    "hits1_other": round(100*h_oth, 2), "mean_log2_ratio": round(d.mean(), 4), "lo": round(lo, 4), "hi": round(hi, 4)})
    with open(a.out + "_retrieval_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build-store", "pairs", "retrieval"])
    ap.add_argument("--store"); ap.add_argument("--labels"); ap.add_argument("--author-labels"); ap.add_argument("--hpo-dir", required=True)
    ap.add_argument("--cases", action="append"); ap.add_argument("--out", required=True)
    ap.add_argument("--r1", default="MENA"); ap.add_argument("--r2", default="EUR")
    ap.add_argument("--authors", help="TSV pmid<TAB>author;author: also exclude pairs/library cases sharing an author")
    ap.add_argument("--cache"); ap.add_argument("--sim", choices=["z", "raw"], default="z")
    ap.add_argument("--perms", type=int, default=1000); ap.add_argument("--k", type=int, default=5); ap.add_argument("--reps", type=int, default=5)
    a = ap.parse_args()
    {"build-store": cmd_build_store, "pairs": cmd_pairs, "retrieval": cmd_retrieval}[a.cmd](a)


if __name__ == "__main__":
    main()
