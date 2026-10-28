#!/usr/bin/env python3
"""Pooled within-gene test of patient-origin effects in the Phenopacket Store.

Outcome: log2 rank (and Hits@1) of the causative gene under the leave-publication-out (LPO)
reference. Gene fixed effects (only within-gene contrasts identify region effects), covariates
log(n_terms), mean IC and HPOA citation of the publication, standard errors clustered by publication.
Model 2 adds the region of the first author's institution (affiliation proxy), so that patient
origin and recording institution are estimated jointly where they are crossed.
Usage: pooled_origin_model.py --results results
"""
import argparse, csv, os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf

ap = argparse.ArgumentParser(); ap.add_argument("--results", required=True); a = ap.parse_args()
R = a.results
lab = {r["pmid"]: r["region"] for r in csv.DictReader(open(os.path.join(R, "labels_llm.tsv")), delimiter="\t")}
aut = {r["pmid"]: r["region"] for r in csv.DictReader(open(os.path.join(R, "labels_proxy.tsv")), delimiter="\t")}
base = {r["id"]: r for r in csv.DictReader(open(os.path.join(R, "baseline_percase.csv"))) if r["ic_type"] == "extrinsic" and r["source"] == "Store"}
rows = []
for r in csv.DictReader(open(os.path.join(R, "lpo_percase.csv"))):
    b = base.get(r["id"])
    if b is None:
        continue
    rows.append(dict(id=r["id"], pmid=r["pmid"], gene=r["gene"], patient=lab.get(r["pmid"], "UNKNOWN"),
                     author=aut.get(r["pmid"], "UNKNOWN"), circ=int(r["pmid_in_hpoa"]),
                     lrank_lpo=np.log2(int(r["rank_lpo"])), lrank_full=np.log2(int(r["rank_full"])),
                     hit_lpo=int(int(r["rank_lpo"]) == 1), hit_full=int(int(r["rank_full"]) == 1),
                     lterms=np.log(int(b["n_terms"])), ic=float(b["mean_ic"])))
df = pd.DataFrame(rows)
out = []

def fit(name, d, formula):
    # keep genes with >=2 distinct values of the contrast variable (others carry no within-gene information)
    m = smf.ols(formula, data=d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d["pmid"])[0]})
    for k in m.params.index:
        if k.startswith("C(patient") or k.startswith("C(author") or k.startswith("C(cross"):
            lo, hi = m.conf_int().loc[k]
            print(f"  {name:38s} {k:60s} {m.params[k]:+.3f} [{lo:+.3f}, {hi:+.3f}] p={m.pvalues[k]:.3g}")
            out.append(dict(model=name, term=k, coef=round(m.params[k], 4), lo=round(lo, 4), hi=round(hi, 4), p=m.pvalues[k], n=int(m.nobs), papers=d["pmid"].nunique(), genes=d["gene"].nunique()))
    return m

def informative(d, col):
    g = d.groupby("gene")[col].nunique()
    return d[d["gene"].isin(g[g >= 2].index)]

known = df[df["patient"].isin(["EUR", "MENA", "EASIA", "OTHER"])]
d1 = informative(known, "patient")
print(f"Model 1 (patient region, within gene): n={len(d1)} papers={d1.pmid.nunique()} genes={d1.gene.nunique()} "
      f"MENA={sum(d1.patient=='MENA')} EASIA={sum(d1.patient=='EASIA')} EUR={sum(d1.patient=='EUR')}")
cov = "lterms + ic + circ + C(gene)"
for y in ("lrank_lpo", "hit_lpo", "lrank_full"):
    fit(f"M1 {y}", d1, f"{y} ~ C(patient, Treatment('EUR')) + {cov}")

# Model 2: patient and author region jointly (crossed where European authors report MENA patients)
k2 = known[known["author"].isin(["EUR", "MENA", "EASIA", "OTHER"])]
d2 = informative(k2, "patient")
ct = pd.crosstab(d2.author, d2.patient)
print("Model 2 author x patient cases (informative genes):\n" + ct.to_string())
for y in ("lrank_lpo", "hit_lpo"):
    fit(f"M2 {y}", d2, f"{y} ~ C(patient, Treatment('EUR')) + C(author, Treatment('EUR')) + {cov}")

# Model 3: MENA patients only, European vs MENA authors (recorder effect with population fixed), no gene FE needed but kept
m3 = k2[(k2.patient == "MENA") & k2.author.isin(["EUR", "MENA"])].copy()
m3["cross"] = m3["author"]
print(f"Model 3 MENA patients: EUR authors {sum(m3.author=='EUR')} cases/{m3[m3.author=='EUR'].pmid.nunique()} papers, MENA authors {sum(m3.author=='MENA')}/{m3[m3.author=='MENA'].pmid.nunique()}")
for y in ("lrank_lpo", "hit_lpo"):
    fit(f"M3 {y} (no gene FE)", m3, f"{y} ~ C(cross, Treatment('MENA')) + lterms + ic + circ")
pd.DataFrame(out).to_csv(os.path.join(R, "origin_pooled_model.csv"), index=False)
