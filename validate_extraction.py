#!/usr/bin/env python3
"""Validate Qwen phenotype extraction (PMC-Patients) against Phenopacket Store curation.

For publications present in both, compare the union of observed HPO terms per publication and
gene (Store curators vs Qwen+mapping). Reports exact precision/recall, 'lenient' precision/recall
(a term counts if it equals, or is an ancestor or descendant of, a term on the other side), and
the mean case-level depth on both sides, overall and by patient region (to check that extraction
quality does not differ by region). Usage:
  validate_extraction.py --store-cases cases_store.tsv --pmcp-cases cases_pmcp.tsv --hpo-dir H --out PREFIX
"""
import argparse, csv, os, sys, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rank_engine import parse_obo

ap = argparse.ArgumentParser()
ap.add_argument("--store-cases", required=True); ap.add_argument("--pmcp-cases", required=True)
ap.add_argument("--hpo-dir", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
terms, alt = parse_obo(os.path.join(a.hpo_dir, "hp.obo"))
anc = {}
def get(t):
    if t in anc: return anc[t]
    s = {t}
    for p in terms[t]["is_a"]:
        if p in terms: s |= get(p)
    anc[t] = s; return s
sys.setrecursionlimit(100000)

def load(p):
    by = collections.defaultdict(lambda: {"terms": set(), "n": 0, "depth": [], "region": None})
    for r in csv.DictReader(open(p), delimiter="\t"):
        ts = [t for t in r["terms"].split(";") if t in terms]
        k = (r["pmid"], r["gene"]); by[k]["terms"] |= set(ts); by[k]["n"] += 1; by[k]["depth"].append(len(ts)); by[k]["region"] = r["region"]
    return by

S, Q = load(a.store_cases), load(a.pmcp_cases)
shared = sorted(set(S) & set(Q))
print(f"publication-gene units in both: {len(shared)} ({len({k[0] for k in shared})} publications)")
rows = []
for k in shared:
    s, q = S[k]["terms"], Q[k]["terms"]
    def lenient_hits(A, B):
        return sum(1 for t in A if any(t == u or t in get(u) or u in get(t) for u in B))
    rows.append({"pmid": k[0], "gene": k[1], "region": S[k]["region"], "store_cases": S[k]["n"], "qwen_cases": Q[k]["n"],
                 "store_terms": len(s), "qwen_terms": len(q),
                 "prec_exact": len(s & q) / len(q) if q else np.nan, "rec_exact": len(s & q) / len(s) if s else np.nan,
                 "prec_lenient": lenient_hits(q, s) / len(q) if q else np.nan, "rec_lenient": lenient_hits(s, q) / len(s) if s else np.nan,
                 "store_depth": np.mean(S[k]["depth"]), "qwen_depth": np.mean(Q[k]["depth"])})
with open(a.out + "_validation.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
def summ(name, R):
    if not R: return
    m = lambda c: np.nanmean([r[c] for r in R])
    print(f"  {name:10s} units {len(R):4d}  precision exact {m('prec_exact'):.2f} lenient {m('prec_lenient'):.2f}  "
          f"recall exact {m('rec_exact'):.2f} lenient {m('rec_lenient'):.2f}  depth store {m('store_depth'):.1f} qwen {m('qwen_depth'):.1f}  "
          f"cases store {sum(r['store_cases'] for r in R)} qwen {sum(r['qwen_cases'] for r in R)}")
summ("all", rows)
for reg in ("EUR", "EASIA", "MENA", "OTHER", "MIXED", "UNKNOWN"):
    summ(reg, [r for r in rows if r["region"] == reg])
