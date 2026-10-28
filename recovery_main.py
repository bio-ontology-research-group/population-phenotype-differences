#!/usr/bin/env python3
"""Is the causative gene recovered less well, by nearest-case phenotype similarity, for patients of one population?

Case-to-case, no disease profiles: every gene is scored by its best-matching library case
(max similarity), library = cases of the library corpora from any region, excluding cases
from the query's publication or from publications sharing an author. Outcome per query:
log2 rank of the causative gene (also Hits@1/@10). Statistic: within (gene, query corpus)
strata with queries from both regions, mean log2 rank of r1 queries minus r2 queries,
averaged over strata (positive = r1 recovered worse). Nulls permute region labels between
reporting groups within (gene, corpus) ('total') or within (gene, corpus, author region)
('patient-carried'). Labels do not enter the ranks, so permutations are cheap.
Usage: recovery_main.py --cases ... --hpo-dir D --cache S.npy --authors A.tsv --groups G.tsv
       --r1 MENA --r2 EUR --lib-corpus store,pmcp,pavs_lit,pubmed_abs [--query-corpus ...] --out PREFIX
"""
import argparse, collections, csv
import numpy as np
import case_pairs as cp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", action="append", required=True); ap.add_argument("--hpo-dir", required=True)
    ap.add_argument("--cache"); ap.add_argument("--authors"); ap.add_argument("--sim", choices=["z", "raw"], default="z")
    ap.add_argument("--groups", required=True); ap.add_argument("--r1", default="MENA"); ap.add_argument("--r2", default="EUR")
    ap.add_argument("--lib-corpus", required=True); ap.add_argument("--query-corpus")
    ap.add_argument("--perms", type=int, default=2000); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sim, cases = cp.setup(a)
    GRP = {r[0]: r[1] for r in csv.reader(open(a.groups), delimiter="\t")}
    genes = sorted({c["gene"] for c in cases}); gidx = {g: i for i, g in enumerate(genes)}
    cg = np.array([gidx[c["gene"]] for c in cases])
    pms = sorted({c["pmid"] for c in cases}); pid = {p: i for i, p in enumerate(pms)}
    cp_ = np.array([pid[c["pmid"]] for c in cases])
    au2p = collections.defaultdict(set)
    for p in pms:
        for x in cp.AUTHORS.get(p, ()): au2p[x].add(pid[p])
    lc = set(a.lib_corpus.split(",")); qc = set(a.query_corpus.split(",")) if a.query_corpus else None
    lib = np.array([c["corpus"] in lc for c in cases])
    order = np.argsort(cg, kind="stable"); starts = np.flatnonzero(np.r_[True, np.diff(cg[order]) != 0])
    assert len(starts) == len(genes)
    rows = []
    for q, c in enumerate(cases):
        if c["region"] not in (a.r1, a.r2) or (qc and c["corpus"] not in qc): continue
        bad = {pid[c["pmid"]]}
        for x in cp.AUTHORS.get(c["pmid"], ()): bad |= au2p[x]
        badp = np.zeros(len(pms), dtype=bool); badp[list(bad)] = True
        valid = lib & ~badp[cp_]
        if not valid[cg == cg[q]].any(): continue
        s = sim.sim_row(q)
        best = np.maximum.reduceat(np.where(valid, s, -np.inf)[order], starts)
        own = best[cg[q]]; others = np.delete(best, cg[q])
        rank = 1 + int(np.sum(others > own)) + 0.5 * int(np.sum(others == own))
        rows.append({"id": c["id"], "pmid": c["pmid"], "gene": c["gene"], "corpus": c["corpus"], "region": c["region"],
                     "author": c["author"], "group": GRP.get(c["pmid"], "P:" + c["pmid"]), "rank": rank,
                     "n_genes": int(np.isfinite(best).sum())})
        if len(rows) % 2000 == 0: print(f"  queries {len(rows)}", flush=True)
    with open(a.out + "_recovery_percase.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    lr = np.log2([r["rank"] for r in rows]); reg = np.array([r["region"] for r in rows])
    for rg in (a.r1, a.r2):
        for co in sorted({r["corpus"] for r in rows}):
            m = (reg == rg) & np.array([r["corpus"] == co for r in rows])
            if m.sum():
                rk = np.array([r["rank"] for r in rows])[m]
                print(f"{rg:6s} {co:12s} n={m.sum():6d} Hits@1 {100*np.mean(rk<=1):5.1f} Hits@10 {100*np.mean(rk<=10):5.1f} median rank {np.median(rk):.0f}")
    strat = collections.defaultdict(list)
    for i, r in enumerate(rows): strat[(r["gene"], r["corpus"])].append(i)
    strat = {k: np.array(v) for k, v in strat.items() if len({rows[i]["region"] for i in v}) == 2}
    def stat(lab):
        d = []
        for idx in strat.values():
            l = lab[idx]; x1, x2 = lr[idx][l == a.r1], lr[idx][l == a.r2]
            if len(x1) and len(x2): d.append(x1.mean() - x2.mean())
        return np.mean(d), len(d)
    t0, ns = stat(reg)
    nq = sum(len(v) for v in strat.values())
    print(f"{a.r1} minus {a.r2}: mean within-stratum difference in log2 rank {t0:+.3f} over {ns} strata ({nq} queries); "
          f"factor {2**t0:.2f} in rank", flush=True)
    rng = np.random.default_rng(11); out = []
    for null in ("total", "patient-carried"):
        units = collections.defaultdict(lambda: collections.defaultdict(list))
        for k, idx in strat.items():
            for i in idx:
                r = rows[i]; key = k if null == "total" else k + (r["author"],)
                units[key][(r["group"], r["region"])].append(i)
        T = []
        for b in range(a.perms):
            lab = reg.copy()
            for us in units.values():
                if len(us) < 2: continue
                ks = list(us); regs = [k[1] for k in ks]; rng.shuffle(regs)
                for k, rr in zip(ks, regs): lab[us[k]] = rr
            T.append(stat(lab)[0])
        T = np.array(T); p = (1 + np.sum(np.abs(T - T.mean()) >= abs(t0 - T.mean()))) / (1 + len(T))
        print(f"  null {null:15s}: mean {T.mean():+.3f} 95% [{np.percentile(T,2.5):+.3f}, {np.percentile(T,97.5):+.3f}]; "
              f"excess {t0 - T.mean():+.3f}; two-sided p={p:.4f}", flush=True)
        out.append({"r1": a.r1, "r2": a.r2, "strata": ns, "queries": nq, "stat": t0, "null": null, "null_mean": T.mean(),
                    "lo": np.percentile(T, 2.5), "hi": np.percentile(T, 97.5), "p_two_sided": p, "perms": len(T)})
    with open(a.out + "_recovery.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)


if __name__ == "__main__":
    main()
