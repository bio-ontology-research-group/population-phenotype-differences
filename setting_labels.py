#!/usr/bin/env python3
"""Ancestry x place-of-care labels for every population, from the phenopackets.
Ancestry region: patient origin country or ancestry text (never the place of care); care region: placeOfCare.
Regions as in origin_analyses.region_of (MENA, EASIA, EUR) and fine_regions (SASIA, LATAM, AFR).
setting = <anc>anc_<care>care for known ancestry, EURorlocal_EURcare for EUR ancestry or unstated ancestry
cared for in EUR (the reference used before), else NA. Writes the case tables with region := setting.

  setting_labels.py --packets results/phenopackets --cases A.tsv [...] --out-dir D
"""
import argparse, collections, csv, json, os, re
from origin_analyses import region_of
from fine_regions import classify
from pmcp_to_cases import ETH


def region(text, words=True):
    t = (text or "").strip()
    if not t or t.lower() == "null": return None
    r = region_of(t)
    if r not in ("OTHER", "UNKNOWN"): return r
    r = classify(t)
    if r: return r
    if words:
        lo = t.lower()
        hits = {reg for reg, ws in ETH.items() if any(re.search(r"\b" + w + r"\b", lo) for w in ws)}
        if len(hits) == 1: return hits.pop()
    return None


def val(v):
    return (v.get("label") or v.get("text") or "") if isinstance(v, dict) else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--packets", required=True); ap.add_argument("--cases", action="append", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    idx = {r["id"]: r for r in csv.DictReader(open(os.path.join(a.packets, "index.tsv")), delimiter="\t")}
    st = collections.Counter()
    for path in a.cases:
        rows = list(csv.DictReader(open(path), delimiter="\t"))
        for r in rows:
            e = {x["name"]: x["value"] for x in json.load(open(os.path.join(a.packets, idx[r["id"]]["file"])))["subject"].get("extensions", [])}
            anc = region(val(e.get("geographicOrigin")), words=False) or region(val(e.get("ancestry")))
            care = region(val(e.get("placeOfCare")), words=False)
            if care and anc: s = "EURorlocal_EURcare" if (anc == "EUR" and care == "EUR") else f"{anc}anc_{care}care"
            elif care == "EUR" and not anc: s = "EURorlocal_EURcare"
            else: s = "NA"
            r["region"] = s; st[s] += 1
        with open(os.path.join(a.out_dir, os.path.basename(path)), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    for k, n in st.most_common(30): print(k, n, sep="\t")


if __name__ == "__main__":
    main()
