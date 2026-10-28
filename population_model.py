#!/usr/bin/env python3
"""Separate provenance from population effects in case-report phenotypes.

Layers
  curation   who turned the report into HPO terms (corpus: Store curators, Qwen on PMC-Patients,
             PAVS curators, PAVS clinical records). Held fixed by corpus fixed effects and by
             permuting only within corpus.
  group      the reporting group (senior/last author; publications of one senior author form one
             group). A population effect must replicate across independent groups, so population
             labels are permuted between groups, never between publications of the same group.
  region     region of the reporting group's institution. Mixes recording conventions and local
             clinical practice (screening, vaccination, referral). Reported both ways:
             'total' effect (author region not held fixed) and 'patient-carried' effect (permutation
             within author region, using patients reported from outside their region).
  population patient origin. Any difference that survives the layers above counts, genetic or not.

Commands
  groups  --authors A.tsv --out G.tsv                 pmid -> group (last author)
  pairs   --cases ... --hpo-dir D --cache S.npy --groups G.tsv --r1 MENA --r2 EUR --perms N --out P
      Pair-level decomposition over all within-gene cross-publication pairs (shared-group pairs
      included, so the group effect is estimated): z-similarity ~ gene FE + depth + corpus-pair FE
      + same group + author-region terms + patient terms (n r1 cases in pair; same population).
      Group-cluster bootstrap CIs for all coefficients; group-level permutation tests for the
      population terms, within (gene, corpus) ['total'] and within (gene, corpus, author region)
      ['patient-carried'].
  allpairs --pops MENA,EUR,... (pairs arguments otherwise)
      The pairs decomposition for every population pair, with the direct contrasts r1r1 vs r1r2 and
      r2r2 vs r1r2; writes OUT_allpairs.json. --vs-rest: each population P against the others pooled (P-P vs P-rest).
  terms   --cases ... --obo hp.obo --termdb T --groups G.tsv --r1 MENA --r2 EUR --perms N --out P
      Per HPO term (inferred closure from an OntoFUNC termdb), Mantel-Haenszel within-(gene, corpus)
      difference in the proportion of r1 vs r2 cases with the term; group-level permutations.
"""
import argparse, collections, csv, json, os
import numpy as np


def load_groups(path):
    return {r[0]: r[1] for r in csv.reader(open(path), delimiter="\t")}


def cmd_groups(a):
    with open(a.out, "w") as f:
        for l in open(a.authors):
            p, _, au = l.rstrip("\n").partition("\t")
            au = [x.strip().lower() for x in au.split(";") if x.strip()]
            f.write(f"{p}\t{('LA:' + au[-1]) if au else 'P:' + p}\n")
    print("written", a.out)


def group_of(G, pmid):
    return G.get(pmid, "P:" + pmid)


def perm_plan(units_by_stratum):
    """units_by_stratum: {stratum: {unit: [case idx]}} -> list of (unit idx lists, labels)"""
    return units_by_stratum


def draw(plan, lab0, rng):
    L = lab0.copy()
    for idx, labs in plan:
        labs = list(labs); rng.shuffle(labs)
        for ii, l in zip(idx, labs): L[ii] = l
    return L


def make_plan(cases, use, key_fn, unit_fn, lab0):
    strata = collections.defaultdict(lambda: collections.defaultdict(list))
    for i in use: strata[key_fn(i)][unit_fn(i)].append(i)
    plan = []
    for units in strata.values():
        # a unit may hold both labels within a stratum: split by label so each unit has one label
        blocks = []
        for idx in units.values():
            by = collections.defaultdict(list)
            for i in idx: by[lab0[i]].append(i)
            blocks += list(by.values())
        if len(blocks) > 1 and len({lab0[b[0]] for b in blocks}) > 1:
            plan.append((blocks, [lab0[b[0]] for b in blocks]))
    return plan


# ---------------------------------------------------------------- pairs
def cmd_pairs(a, pre=None):
    import case_pairs as cp
    rng = np.random.default_rng(17)
    sim, cases = pre or cp.setup(a)
    G = load_groups(a.groups)
    reg = np.array([c["region"] for c in cases]); corp = np.array([c["corpus"] or "?" for c in cases])
    aut = np.array([c["author"] for c in cases]); pm = np.array([c["pmid"] for c in cases])
    gene = np.array([c["gene"] for c in cases]); grp = np.array([group_of(G, p) for p in pm])
    use = np.where(np.isin(reg, [a.r1, a.r2]))[0]
    byg = collections.defaultdict(list)
    for i in use: byg[gene[i]].append(i)
    I, J = [], []
    for g, idx in byg.items():
        if len({reg[i] for i in idx}) < 2: continue
        for k, i in enumerate(idx):
            for j in idx[k + 1:]:
                if pm[i] != pm[j]: I.append(i); J.append(j)
    I = np.array(I); J = np.array(J)
    if len(I) == 0: print("no pairs"); return []
    y = (sim.Z if a.sim == "z" else sim.S)[I, J].astype(float)
    gi = np.unique(gene[I], return_inverse=True)[1]
    print(f"pairs {len(I)} in {gi.max() + 1} genes; same-group pairs {(grp[I] == grp[J]).sum()}", flush=True)
    L = np.log(np.maximum(sim.lens, 1)); li, lj = L[I], L[J]
    cols = {"depth_sum": li + lj, "depth_prod": li * lj, "depth_absdiff": np.abs(li - lj)}
    ppair = np.array(["|".join(sorted((corp[i], corp[j]))) for i, j in zip(I, J)])
    lev = sorted(set(ppair)); base = max(lev, key=lambda v: (ppair == v).sum())
    for v in lev:
        if v != base: cols[f"curation[{v}]"] = (ppair == v).astype(float)
    cols["same_group"] = (grp[I] == grp[J]).astype(float)
    unk = np.isin(aut, ["UNKNOWN", "MIXED", "OTHER"]); ar1 = (aut == a.r1).astype(float)
    cols["author_region_n_r1"] = ar1[I] + ar1[J]
    cols["author_region_same"] = ((aut[I] == aut[J]) & ~unk[I] & ~unk[J]).astype(float)
    cols["author_region_unknown"] = (unk[I] | unk[J]).astype(float)

    def demean(M):
        M = np.atleast_2d(M.T).T
        s = np.zeros((gi.max() + 1, M.shape[1])); np.add.at(s, gi, M)
        return M - (s / np.bincount(gi)[:, None])[gi]

    def popcols(R):
        r1 = (R == a.r1).astype(float)
        return np.column_stack([r1[I] + r1[J], (R[I] == R[J]).astype(float)])
    names_fixed = list(cols)
    F_all = demean(np.column_stack(list(cols.values()))); yd = demean(y)[:, 0]
    F_noaut = demean(np.column_stack([v for k, v in cols.items() if not k.startswith("author_region")]))

    def fit(F, R, w=None):
        X = np.column_stack([F, demean(popcols(R))])
        if w is not None: X = X * w[:, None] ** .5; yy = yd * w ** .5
        else: yy = yd
        return np.linalg.lstsq(X, yy, rcond=None)[0]
    if getattr(a, "r2_only", False):
        # variance decomposition of within-gene pair similarity: share of the total variance between genes, and
        # within genes the R^2 of the full model, of each layer alone and the unique R^2 of each layer (drop one)
        P_ = demean(popcols(reg))
        layers = {"phenotyping depth": [k for k in names_fixed if k.startswith("depth")],
                  "curation process": [k for k in names_fixed if k.startswith("curation")],
                  "reporting group": ["same_group"],
                  "author region": [k for k in names_fixed if k.startswith("author_region")]}
        colidx = {k: j for j, k in enumerate(names_fixed)}

        def r2(blocks, w=None):
            X = np.column_stack(blocks); yy = yd
            if w is not None: X = X * w[:, None] ** .5; yy = yd * w ** .5
            b = np.linalg.lstsq(X, yy, rcond=None)[0]; e = yy - X @ b
            return 1 - (e @ e) / (yy @ yy)

        def parts(w=None):
            blk = {L: F_all[:, [colidx[k] for k in ks]] for L, ks in layers.items()}; blk["population"] = P_
            full = r2(list(blk.values()), w)
            res = {"within-gene R2, full model": full}
            for L in blk:
                res[f"{L}: alone"] = r2([blk[L]], w)
                res[f"{L}: unique"] = full - r2([v for k, v in blk.items() if k != L], w)
            return res
        yc = y - y.mean(); gene_share = 1 - (yd @ yd) / (yc @ yc)
        p0 = parts()
        groups_u = np.unique(grp[use]); gidx = {g: k for k, g in enumerate(groups_u)}
        gI = np.array([gidx[grp[i]] for i in I]); gJ = np.array([gidx[grp[j]] for j in J])
        bs = []
        for _ in range(getattr(a, "r2_boot", 0)):
            m = np.bincount(rng.integers(0, len(groups_u), len(groups_u)), minlength=len(groups_u)).astype(float)
            w = m[gI] * m[gJ]
            if w.sum() > 0: bs.append(parts(w))
        rows_ = [{"pair": f"{a.r1}-{a.r2}", "quantity": "between-gene share of variance", "value": gene_share, "lo": "", "hi": "",
                  "n_pairs": len(I), "n_genes": int(gi.max() + 1)}]
        for k, v in p0.items():
            lo, hi = (np.percentile([b[k] for b in bs], [2.5, 97.5]) if bs else ("", ""))
            rows_.append({"pair": f"{a.r1}-{a.r2}", "quantity": k, "value": v, "lo": lo, "hi": hi, "n_pairs": len(I), "n_genes": int(gi.max() + 1)})
            print(f"   {k:40s} {100 * v:7.3f}%")
        print(f"   between-gene share of variance {100 * gene_share:.2f}%")
        with open(a.out + "_r2.csv", "w", newline="") as f:
            w_ = csv.DictWriter(f, fieldnames=list(rows_[0])); w_.writeheader(); w_.writerows(rows_)
        return []
    pop_names = [f"population_n_{a.r1}", "population_same"]
    b_all = fit(F_all, reg); b_noaut = fit(F_noaut, reg)
    # group-cluster bootstrap (resample groups; pair weight = product of case multiplicities)
    groups_u = np.unique(grp[use]); gidx = {g: k for k, g in enumerate(groups_u)}
    gI = np.array([gidx[grp[i]] for i in I]); gJ = np.array([gidx[grp[j]] for j in J])
    BS = []
    for b in range(a.boot):
        m = np.bincount(rng.integers(0, len(groups_u), len(groups_u)), minlength=len(groups_u)).astype(float)
        w = m[gI] * m[gJ]
        if w.sum() == 0: continue
        BS.append(fit(F_all, reg, w))
    BS = np.array(BS)
    print("\nDecomposition (OLS, gene FE; all layers in one model). Coefficients in z-similarity units;")
    print(f"95% CI from {len(BS)} group-cluster bootstraps. Curation baseline = {base}.")
    rows = []
    for k, nm in enumerate(names_fixed + pop_names):
        lo, hi = np.percentile(BS[:, k], [2.5, 97.5])
        print(f"   {nm:44s} {b_all[k]:+7.3f}  [{lo:+.3f}, {hi:+.3f}]")
        rows.append({"term": nm, "coef": b_all[k], "lo": lo, "hi": hi})
    print(f"   (without author-region terms: population_n_{a.r1} {b_noaut[-2]:+.3f}, population_same {b_noaut[-1]:+.3f})")
    # contrasts: r1-r1 pair minus r1-r2 pair = b_n + b_same; r2-r2 pair minus r1-r2 pair = b_same - b_n
    CON = {f"{a.r1}{a.r1}_vs_{a.r1}{a.r2}": np.array([1., 1.]), f"{a.r2}{a.r2}_vs_{a.r1}{a.r2}": np.array([-1., 1.])}
    for nm, c in CON.items():
        lo, hi = np.percentile(BS[:, -2:] @ c, [2.5, 97.5])
        print(f"   {nm:44s} {b_all[-2:] @ c:+7.3f}  [{lo:+.3f}, {hi:+.3f}]")
        rows.append({"term": nm, "coef": b_all[-2:] @ c, "lo": lo, "hi": hi})
    # permutation tests for the population terms
    out = []
    for label, F, key in (("total", F_noaut, lambda i: (gene[i], corp[i])),
                          ("patient-carried", F_all, lambda i: (gene[i], corp[i], aut[i]))):
        plan = make_plan(cases, use, key, lambda i: grp[i], reg)
        b0 = fit(F, reg)[-2:]
        null = np.array([fit(F, draw(plan, reg, rng))[-2:] for _ in range(a.perms)])
        mu = null.mean(0)
        for k, nm in enumerate(pop_names):
            p = (1 + np.sum(np.abs(null[:, k] - mu[k]) >= abs(b0[k] - mu[k]))) / (1 + len(null))
            print(f"   [{label:15s}] {nm:24s} coef {b0[k]:+.3f}  null mean {mu[k]:+.3f}  excess {b0[k]-mu[k]:+.3f}  "
                  f"null 95% [{np.percentile(null[:,k],2.5):+.3f}, {np.percentile(null[:,k],97.5):+.3f}]  p={p:.4f}  units {sum(len(x[0]) for x in plan)}")
            out.append({"term": f"{label}:{nm}", "coef": b0[k], "lo": mu[k], "hi": p})
        for nm, c in CON.items():
            t0, tn = b0 @ c, null @ c; m = tn.mean()
            p = (1 + np.sum(np.abs(tn - m) >= abs(t0 - m))) / (1 + len(tn))
            print(f"   [{label:15s}] {nm:24s} coef {t0:+.3f}  null mean {m:+.3f}  excess {t0-m:+.3f}  "
                  f"null 95% [{np.percentile(tn,2.5):+.3f}, {np.percentile(tn,97.5):+.3f}]  p={p:.4f}")
            out.append({"term": f"{label}:{nm}", "coef": t0, "lo": m, "hi": p})
    with open(a.out + "_decomposition.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["term", "coef", "lo", "hi"]); w.writeheader(); w.writerows(rows + out)
    n_r = {r: int(((reg[I] == r).sum() + (reg[J] == r).sum())) for r in (a.r1, a.r2)}
    return [dict(r1=a.r1, r2=a.r2, n_pairs=len(I), n_genes=int(gi.max() + 1), **{f"n_case_slots_{k}": v for k, v in n_r.items()},
                 **{x["term"]: x["coef"] for x in rows}, **{x["term"] + "_lo": x["lo"] for x in rows},
                 **{x["term"] + "_hi": x["hi"] for x in rows}, **{x["term"] + "_p": x["hi"] for x in out},
                 **{x["term"] + "_coef": x["coef"] for x in out}, **{x["term"] + "_nullmean": x["lo"] for x in out})]


def _one_pair(rr):
    import contextlib, io
    r1, r2 = rr; a = _AP_ARGS
    a.r1, a.r2, a.out = r1, r2, f"{a.base}_{r1}_{r2}"
    if r2 == "REST":   # forked worker: relabelling is local to this process
        pops = set(a.pops.split(","))
        for c in _AP_PRE[1]:
            if c["region"] in pops and c["region"] != r1: c["region"] = "REST"
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print(f"\n########## {r1} vs {r2}")
        res = cmd_pairs(a, _AP_PRE)
    return buf.getvalue(), res


def cmd_allpairs(a):
    """Every population against every other one, loading the similarity matrix once; the pairs run
    in forked worker processes that share the matrix read-only."""
    import case_pairs as cp, itertools, multiprocessing as mp
    global _AP_ARGS, _AP_PRE
    _AP_PRE = cp.setup(a); a.base = a.out; _AP_ARGS = a
    combos = ([(p_, "REST") for p_ in a.pops.split(",")] if a.vs_rest else list(itertools.combinations(a.pops.split(","), 2)))
    res = []
    with mp.get_context("fork").Pool(max(1, a.pair_workers)) as pool:
        for txt, r in pool.imap_unordered(_one_pair, combos):
            print(txt, flush=True); res += r
            with open(a.base + "_allpairs.json", "w") as f: json.dump(res, f, indent=1, default=float)


# ---------------------------------------------------------------- terms
def cmd_terms(a):
    from rank_engine import parse_obo
    from phenotype_specificity import names, ancestors_fn, ROOT_PA
    rng = np.random.default_rng(19)
    where = [tuple(x.split("=", 1)) for x in (a.where or [])]
    G = load_groups(a.groups)
    terms, alt = parse_obo(a.obo); nm = names(a.obo)
    anc = ancestors_fn(terms, os.path.join(a.termdb, "graph_path.txt") if a.termdb else None)
    cases = []
    for path in a.cases:
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["region"] not in (a.r1, a.r2): continue
            if any(r.get(k_, "") != v_ for k_, v_ in where): continue
            if a.key_extra and not r.get(a.key_extra): continue
            ts = [alt.get(t, t) for t in r["terms"].split(";")]
            ts = [t for t in ts if t in terms and not terms[t]["obsolete"]]
            if ts: cases.append({"gene": r["gene"], "region": r["region"], "corpus": r.get("corpus", ""),
                                 "author": r.get("author_region") or "UNKNOWN", "group": group_of(G, r["pmid"]), "terms": ts,
                                 "extra": r.get(a.key_extra, "") if a.key_extra else ""})
    key = (lambda c: (c["gene"], c["corpus"], c["extra"])) if a.strata == "total" else (lambda c: (c["gene"], c["corpus"], c["author"], c["extra"]))
    sreg = collections.defaultdict(set)
    for c in cases: sreg[key(c)].add(c["region"])
    cases = [c for c in cases if len(sreg[key(c)]) == 2]
    print(f"{len(cases)} cases ({sum(c['region'] == a.r1 for c in cases)} {a.r1}) in {len({key(c) for c in cases})} strata; "
          f"corpora {dict(collections.Counter(c['corpus'] for c in cases))}; groups {len({c['group'] for c in cases})}")
    prop = [{t for x in c["terms"] for t in anc(x)} for c in cases]
    cnt = collections.Counter(t for p in prop for t in p)
    T = sorted(t for t, n in cnt.items() if n >= a.min_cases and t != ROOT_PA and ROOT_PA in anc(t))
    col = {t: j for j, t in enumerate(T)}
    X = np.zeros((len(cases), len(T)), np.float32)
    for i, p in enumerate(prop):
        for t in p:
            if t in col: X[i, col[t]] = 1
    depth = np.log([len(c["terms"]) for c in cases])
    extra = {}
    if a.termclass:
        tc = {x["id"]: x for x in map(json.loads, open(a.termclass))}
        TEST = {"laboratory", "imaging", "physiological_test", "pathology"}
        def cnt(c, pred): return sum(1 for t in set(c["terms"]) if t in tc and pred(tc[t]))
        nt = np.array([cnt(c, lambda x: x["class"] in TEST) for c in cases], float)
        nb = np.array([cnt(c, lambda x: x["class"] in ("clinical_exam", "symptom_history")) for c in cases], float)
        extra["share of findings needing a test"] = nt / np.maximum(nt + nb, 1)
        extra["log(1+ test-based findings)"] = np.log1p(nt)
        extra["log(1+ bedside/history findings)"] = np.log1p(nb)
        for k_ in ("laboratory", "imaging", "physiological_test", "pathology", "clinical_exam", "symptom_history"):
            extra[f"log(1+ {k_})"] = np.log1p(np.array([cnt(c, lambda x, k_=k_: x["class"] == k_) for c in cases], float))
        # per organ system (children of the phenotypic-abnormality root): findings in that system, and their test share
        BR = sorted(t for t in terms if ROOT_PA in terms[t]["is_a"] and not terms[t]["obsolete"])
        inb = [{b_ for t in set(c["terms"]) if t in tc for b_ in BR if b_ in anc(t)} for c in cases]
        for b_ in BR:
            nb_ = sum(b_ in x for x in inb)
            if nb_ < 200: continue
            bn = nm.get(b_, b_).replace("Abnormality of the ", "").replace("Abnormality of ", "")
            tt = np.array([sum(1 for t in set(c["terms"]) if t in tc and b_ in anc(t) and tc[t]["class"] in TEST) for c in cases], float)
            bb = np.array([sum(1 for t in set(c["terms"]) if t in tc and b_ in anc(t) and tc[t]["class"] in ("clinical_exam", "symptom_history")) for c in cases], float)
            extra[f"organ[{bn}] log(1+ findings)"] = np.log1p(tt + bb)
            extra[f"organ[{bn}] log(1+ test-based)"] = np.log1p(tt)
            extra[f"organ[{bn}] log(1+ bedside)"] = np.log1p(bb)
    sid = {s: k for k, s in enumerate(sorted({key(c) for c in cases}))}; S = np.array([sid[key(c)] for c in cases]); ns = len(sid)
    lab0 = np.array([c["region"] for c in cases])
    plan = make_plan(cases, range(len(cases)), lambda i: key(cases[i]), lambda i: cases[i]["group"], lab0)

    EX = np.column_stack(list(extra.values())) if extra else np.zeros((len(cases), 0))
    def stat(R):
        L1 = (R == a.r1).astype(float)
        n1 = np.bincount(S, weights=L1, minlength=ns); n = np.bincount(S, minlength=ns); n2 = n - n1
        ok = (n1 > 0) & (n2 > 0); w = np.where(ok, n1 * n2 / np.maximum(n, 1), 0)
        c_ = np.where(L1 == 1, w[S] / np.maximum(n1[S], 1), -w[S] / np.maximum(n2[S], 1)) / w.sum()
        return c_ @ X, c_ @ depth, c_ @ EX
    t0, d0, e0 = stat(lab0)
    null = np.empty((a.perms, len(T)), np.float32); dn = np.empty(a.perms); en = np.empty((a.perms, EX.shape[1]))
    for b in range(a.perms): null[b], dn[b], en[b] = stat(draw(plan, lab0, rng))
    if extra:
        print(f"case-level outcomes ({a.r1} minus {a.r2}, MH within strata; group-permutation p):")
        orows = []
        for k_, (nm_, v0) in enumerate(zip(extra, e0)):
            mu_ = en[:, k_].mean(); p_ = (1 + np.sum(np.abs(en[:, k_] - mu_) >= abs(v0 - mu_))) / (1 + a.perms)
            lo_, hi_ = np.percentile(en[:, k_], [2.5, 97.5])
            print(f"   {nm_:36s} {v0:+.4f}  null 95% [{lo_:+.4f}, {hi_:+.4f}]  p={p_:.4f}")
            orows.append({"outcome": nm_, "diff": v0, "null_mean": mu_, "null_lo": lo_, "null_hi": hi_, "p": p_,
                          "n_cases": len(cases), "n_r1": int((lab0 == a.r1).sum()), "strata": ns, "groups": len({c["group"] for c in cases})})
        with open(a.out + "_outcomes.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(orows[0])); w.writeheader(); w.writerows(orows)
    mu, sd = null.mean(0), null.std(0) + 1e-9; z = (t0 - mu) / sd
    p = (1 + (np.abs(null - mu) >= np.abs(t0 - mu)).sum(0)) / (1 + a.perms)
    maxz = np.abs((null - mu) / sd).max(1); fwer = (1 + (maxz[None, :] >= np.abs(z)[:, None]).sum(1)) / (1 + a.perms)
    o = np.argsort(p); q = np.empty_like(p); m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    zn = (null - mu) / sd; g0 = (z ** 2).sum(); gn = (zn ** 2).sum(1)
    print(f"depth diff {d0:+.3f} (null 95% [{np.percentile(dn,2.5):+.3f}, {np.percentile(dn,97.5):+.3f}])")
    print(f"terms {len(T)}; p<0.05 {int((p<0.05).sum())} (exp {0.05*len(T):.0f}); FDR<0.1 {int((q<0.1).sum())}; FWER<0.05 {int((fwer<0.05).sum())}; "
          f"global sum z^2 {g0:.0f} vs null {gn.mean():.0f}, p={(1 + (gn >= g0).sum()) / (1 + a.perms):.4f}")
    # organ-system branches (children of HP:0000118): mean term z, tested against the same group permutations
    br = sorted(t for t in terms if ROOT_PA in terms[t]["is_a"] and not terms[t]["obsolete"])
    bres = []
    for b_ in br:
        m = np.array([b_ in anc(t) for t in T])
        if m.sum() < 5: continue
        o0 = z[m].mean(); on = zn[:, m].mean(1)
        bres.append((b_, m.sum(), o0, on))
    if bres:
        sdn = np.array([x[3].std() for x in bres]); mx = np.abs(np.array([x[3] for x in bres]) / sdn[:, None]).max(0)
        print(f"organ-system branches (mean term z; group-permutation p; FWER over {len(bres)} branches):")
        with open(a.out + "_branches.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["hpo", "name", "n_terms", "mean_z", "p", "fwer"])
            for (b_, n_, o0, on), s_ in sorted(zip(bres, sdn), key=lambda x: -abs(x[0][2] / x[1])):
                pb = (1 + np.sum(np.abs(on) >= abs(o0))) / (1 + len(on)); fw = (1 + np.sum(mx >= abs(o0) / s_)) / (1 + len(mx))
                w.writerow([b_, nm.get(b_, ""), n_, f"{o0:.3f}", f"{pb:.4g}", f"{fw:.4g}"])
                if pb < 0.1: print(f"  {b_} {nm.get(b_,'')[:50]:50s} terms {n_:4d} mean z {o0:+.3f} p {pb:.4f} fwer {fw:.4f}")
    rows = sorted(zip(T, t0, z, p, q, fwer, X.sum(0)), key=lambda r: r[3])
    with open(a.out + "_terms.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["hpo", "name", "mh_diff", "z", "p", "q_bh", "fwer", "n_cases"])
        for r in rows: w.writerow([r[0], nm.get(r[0], ""), f"{r[1]:.4f}", f"{r[2]:.3f}", f"{r[3]:.4g}", f"{r[4]:.4g}", f"{r[5]:.4g}", int(r[6])])
    for r in rows[:a.show]:
        print(f"  {r[0]} {nm.get(r[0],'')[:48]:48s} diff {r[1]:+.3f} z {r[2]:+.2f} p {r[3]:.4f} q {r[4]:.3f} fwer {r[5]:.3f} n {int(r[6])}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["groups", "pairs", "allpairs", "terms"])
    ap.add_argument("--authors"); ap.add_argument("--groups"); ap.add_argument("--out", required=True)
    ap.add_argument("--cases", action="append"); ap.add_argument("--hpo-dir"); ap.add_argument("--cache")
    ap.add_argument("--sim", choices=["z", "raw"], default="z"); ap.add_argument("--obo"); ap.add_argument("--termdb")
    ap.add_argument("--pops", default="MENA,EUR,EASIA,SASIA,LATAM,AFR")
    ap.add_argument("--workers", type=int, default=1, help="processes for building the similarity matrix")
    ap.add_argument("--vs-rest", action="store_true", help="allpairs: each population against all other listed populations pooled")
    ap.add_argument("--pair-workers", type=int, default=1, help="allpairs: population pairs run in parallel")
    ap.add_argument("--r1", default="MENA"); ap.add_argument("--r2", default="EUR")
    ap.add_argument("--strata", choices=["total", "patient-carried"], default="total")
    ap.add_argument("--perms", type=int, default=1000); ap.add_argument("--boot", type=int, default=200)
    ap.add_argument("--r2-only", action="store_true", help="pairs/allpairs: only the variance decomposition (OUT_r2.csv)")
    ap.add_argument("--r2-boot", type=int, default=50, help="group-cluster bootstraps for the variance decomposition")
    ap.add_argument("--termclass", help="termclass.jsonl from qwen_tasks.py: adds case-level test-vs-bedside outcomes")
    ap.add_argument("--where", action="append", help="terms: keep only cases with COL=VALUE (repeatable)")
    ap.add_argument("--key-extra", help="terms: add this case-table column to the stratum key (cases without it dropped)")
    ap.add_argument("--min-cases", type=int, default=10); ap.add_argument("--show", type=int, default=15)
    a = ap.parse_args()
    a.authors_for_related = None
    {"groups": cmd_groups, "pairs": cmd_pairs, "allpairs": cmd_allpairs, "terms": cmd_terms}[a.cmd](a)


if __name__ == "__main__":
    main()
