#!/usr/bin/env python3
"""Normalise a reviewer's returned sheets (prepare_human_review.py packets) before scoring.

Trims and lower-cases document ids, status and observation basis; writes one row per HPO identifier when a cell
lists several; maps a basis cell with several values to the single stated value if the others are "unclear",
else to "mixed"; drops empty rows.  The returned file is not changed.

  python3 clinician_review_normalize.py --in findings_MA.csv --out findings_MA.norm.csv
"""
import argparse, csv, re

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--in", dest="inp", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
rows = list(csv.DictReader(open(a.inp, encoding="utf-8-sig"))); out = []
for r in rows:
    r = {k: (v or "").strip() for k, v in r.items() if k}
    if not any(r.values()): continue
    r["document_id"] = r["document_id"].upper()
    if "status" in r: r["status"] = r["status"].lower().replace(" ", "_")
    if "observation_basis" in r:
        b = [x.strip().replace(" ", "_") for x in re.split(r"[,;/|]", r["observation_basis"].lower()) if x.strip()]
        known = sorted({x for x in b if x != "unclear"})
        r["observation_basis"] = known[0] if len(known) == 1 else ("mixed" if known else ("unclear" if b else ""))
    ids = re.findall(r"HP[:_]\d{7}", r.get("hpo_id", ""))
    if "hpo_id" in r and len(ids) > 0:
        for i in ids: out.append(dict(r, hpo_id=i.replace("_", ":")))
    else: out.append(r)
with open(a.out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
print(len(rows), "rows in,", len(out), "rows out")
