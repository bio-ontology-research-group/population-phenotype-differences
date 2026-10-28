#!/usr/bin/env python3
"""Does the share of test-based findings follow the resources of the country where the patient was cared for?

Cases: PubMed abstracts and PMC-Patients of the six populations with a single place-of-care country (phenopackets).
Country covariates (World Bank, value in the publication year, nearest available year 2000-2023):
current health expenditure per capita (SH.XPD.CHEX.PC.CD), GDP per capita PPP (NY.GDP.PCAP.PP.CD),
physicians per 1,000 (SH.MED.PHYS.ZS); entered as log10.
Outcomes per case (finding classes from termclass.jsonl, as in population_model.py terms): share of findings
needing a test, log(1+ test-based), log(1+ bedside/history).
Model: outcome ~ covariate with gene x corpus fixed effects (within-stratum demeaning). Inference at country level:
the covariate values are permuted between countries (each country keeps its cases), 2,000 permutations; 95% CI by
bootstrap over countries. A second model adds population dummies (vs EUR) to ask how much of each population
difference the covariate absorbs.

  healthspend.py --index results/phenopackets/index.tsv --cases A.tsv [...] --termclass results/termclass.jsonl
                 --wb-dir results/v4/worldbank --years-abs abstracts.jsonl --years-tsv pmcp_pmid_country.tsv --out P
"""
import argparse, collections, csv, json, os
import numpy as np

ALIAS = {"United States of America": "United States", "USA": "United States", "Turkey": "Turkiye", "Türkiye": "Turkiye",
         "South Korea": "Korea, Rep.", "Korea": "Korea, Rep.", "Iran": "Iran, Islamic Rep.", "Czech Republic": "Czechia",
         "Russia": "Russian Federation", "Egypt": "Egypt, Arab Rep.", "Palestine": "West Bank and Gaza",
         "Slovakia": "Slovak Republic", "Syria": "Syrian Arab Republic", "Hong Kong": "Hong Kong SAR, China",
         "Puerto Rico": "Puerto Rico (US)", "Yemen": "Yemen, Rep.", "Macedonia": "North Macedonia",
         "Republic of Macedonia": "North Macedonia", "Venezuela": "Venezuela, RB", "Laos": "Lao PDR",
         "Gambia": "Gambia, The", "Democratic Republic of the Congo": "Congo, Dem. Rep.", "England": "United Kingdom",
         "The Netherlands": "Netherlands", "Bosnia and Hercegovina": "Bosnia and Herzegovina"}
IND = {"health_exp_pc": "SH.XPD.CHEX.PC.CD", "gdp_pc_ppp": "NY.GDP.PCAP.PP.CD", "physicians_per_1000": "SH.MED.PHYS.ZS"}
TEST = {"laboratory", "imaging", "physiological_test", "pathology"}


def load_wb(d):
    names = {x["name"]: x["id"] for x in json.load(open(os.path.join(d, "countries.json")))[1]}
    vals = {}
    for k, code in IND.items():
        v = collections.defaultdict(dict)
        for r in json.load(open(os.path.join(d, code + ".json")))[1]:
            if r["value"] is not None: v[r["countryiso3code"]][int(r["date"])] = float(r["value"])
        vals[k] = v
    return names, vals


def at_year(series, y):
    if not series: return None
    yrs = np.array(sorted(series)); y = min(max(y, yrs[0]), yrs[-1])
    return series[int(yrs[np.argmin(np.abs(yrs - y))])]


def fe_resid(y, S):
    m = np.bincount(S, weights=y) / np.bincount(S)
    return y - m[S]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True); ap.add_argument("--cases", action="append", required=True)
    ap.add_argument("--termclass", required=True); ap.add_argument("--wb-dir", required=True)
    ap.add_argument("--years-abs"); ap.add_argument("--years-tsv", action="append", default=[])
    ap.add_argument("--perms", type=int, default=2000); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(23)
    names, wb = load_wb(a.wb_dir)
    year = {}
    if a.years_abs:
        for l in open(a.years_abs):
            d = json.loads(l)
            if str(d.get("year", "")).isdigit(): year[d["pmid"]] = int(d["year"])
    for f in a.years_tsv:
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r.get("year", "").isdigit(): year.setdefault(r["pmid"], int(r["year"]))
    care = {r["id"]: r["place_of_care"] for r in csv.DictReader(open(a.index), delimiter="\t")}
    tc = {x["id"]: x["class"] for x in map(json.loads, open(a.termclass))}
    rows, miss = [], collections.Counter()
    for path in a.cases:
        for r in csv.DictReader(open(path), delimiter="\t"):
            c = care.get(r["id"], "")
            iso = names.get(ALIAS.get(c, c))
            if r["region"] not in ("EUR", "MENA", "EASIA", "SASIA", "LATAM", "AFR"): continue
            if not iso: miss[c] += 1; continue
            ts = set(r["terms"].split(";")); nt = sum(tc.get(t) in TEST for t in ts)
            nb = sum(tc.get(t) in ("clinical_exam", "symptom_history") for t in ts)
            if nt + nb == 0: continue
            y = year.get(r["pmid"], 2015)
            cov = {k: at_year(wb[k].get(iso), y) for k in IND}
            rows.append({"id": r["id"], "gene": r["gene"], "corpus": r["corpus"], "region": r["region"], "iso": iso,
                         "country": ALIAS.get(c, c), "year": y, "share": nt / (nt + nb), "logt": np.log1p(nt),
                         "logb": np.log1p(nb), **cov})
    print(f"cases {len(rows)}; countries {len({r['iso'] for r in rows})}; unmatched places of care: {sum(miss.values())} "
          f"(top {miss.most_common(6)})")
    out = []
    for k in IND:
        R = [r for r in rows if r[k] is not None and r[k] > 0]
        key = [(r["gene"], r["corpus"]) for r in R]
        sz = collections.Counter(key); keep = [i for i, kk in enumerate(key) if sz[kk] > 1]
        R = [R[i] for i in keep]
        sid = {s: j for j, s in enumerate(sorted({(r["gene"], r["corpus"]) for r in R}))}
        S = np.array([sid[(r["gene"], r["corpus"])] for r in R])
        cty = sorted({r["iso"] for r in R}); cid = {c: j for j, c in enumerate(cty)}
        C = np.array([cid[r["iso"]] for r in R])
        # country-level covariate: one value per country (mean of its case-level values), so permutation swaps countries
        xc = np.array([np.mean([np.log10(r[k]) for r in R if r["iso"] == c]) for c in cty])
        regs = sorted({r["region"] for r in R} & {"MENA", "EASIA", "SASIA", "LATAM", "AFR"})
        D = np.column_stack([[r["region"] == g for r in R] for g in regs]).astype(float)
        Dr = np.column_stack([fe_resid(D[:, j], S) for j in range(D.shape[1])])
        for yk, ylab in (("share", "share of findings needing a test"), ("logt", "log(1+ test-based)"), ("logb", "log(1+ bedside/history)")):
            y = fe_resid(np.array([r[yk] for r in R]), S)
            def slope(xcv):
                x = fe_resid(xcv[C], S); return (x @ y) / (x @ x)
            b0 = slope(xc)
            null = np.array([slope(rng.permutation(xc)) for _ in range(a.perms)])
            p = (1 + np.sum(np.abs(null) >= abs(b0))) / (1 + a.perms)
            bs = []
            for _ in range(500):
                pick = rng.integers(0, len(cty), len(cty)); w = np.bincount(pick, minlength=len(cty))[C].astype(float)
                if w.sum() == 0: continue
                x = fe_resid(xc[C], S); bs.append((w * x) @ y / ((w * x) @ x))
            lo, hi = np.percentile(bs, [2.5, 97.5])
            # population dummies with and without the covariate
            g0 = np.linalg.lstsq(Dr, y, rcond=None)[0]
            g1 = np.linalg.lstsq(np.column_stack([Dr, fe_resid(xc[C], S)]), y, rcond=None)[0]
            print(f"{k:20s} {ylab:34s} slope per 10-fold {b0:+.4f} [{lo:+.4f}, {hi:+.4f}] country-perm p={p:.4f}  "
                  f"(n={len(R)}, countries={len(cty)}, strata={len(sid)})")
            print("    population vs EUR without -> with covariate: " +
                  ", ".join(f"{g} {g0[j]:+.4f}->{g1[j]:+.4f}" for j, g in enumerate(regs)))
            out.append({"covariate": k, "outcome": ylab, "slope_per_10fold": b0, "lo": lo, "hi": hi, "p": p, "n": len(R),
                        "countries": len(cty), "strata": len(sid),
                        **{f"{g}_without": g0[j] for j, g in enumerate(regs)}, **{f"{g}_with": g1[j] for j, g in enumerate(regs)}})
        if k == "health_exp_pc":
            ys = fe_resid(np.array([r["share"] for r in R]), S)
            with open(a.out + "_country_points.csv", "w", newline="") as f:
                w = csv.writer(f); w.writerow(["iso", "country", "n", "log10_health_exp_pc", "mean_resid_test_share", "region_mode"])
                for c in cty:
                    m = C == cid[c]
                    reg = collections.Counter(r["region"] for r, mm in zip(R, m) if mm).most_common(1)[0][0]
                    w.writerow([c, next(r["country"] for r in R if r["iso"] == c), int(m.sum()), f"{xc[cid[c]]:.4f}", f"{ys[m].mean():.5f}", reg])
    with open(a.out + "_slopes.csv", "w", newline="") as f:
        fields = sorted({k for o in out for k in o}, key=lambda x: list(out[0]).index(x) if x in out[0] else 99)
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(out)


if __name__ == "__main__":
    main()
