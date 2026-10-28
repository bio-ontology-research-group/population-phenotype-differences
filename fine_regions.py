#!/usr/bin/env python3
"""Split the OTHER population bin into South Asia (SASIA), Latin America (LATAM) and sub-Saharan
Africa (AFR). MENA, EUR and EASIA labels are kept as they are (Pakistan and Afghanistan stay MENA,
as in origin_analyses.region_of). The evidence for an OTHER case is read from its phenopacket, in
the order the original label used it: patient origin country, ancestry text, then place of care
(only when the label was assigned from the place of care), then paper-level origins (Store).
Cases that fit none of the three stay OTHER. The author region (first-author affiliation country)
is split the same way, from the pmid -> first_country tables given with --author-country.

  fine_regions.py --packets results/phenopackets --cases A.tsv [--cases B.tsv ...] --out-dir D
writes D/<basename of each case table> with the region column relabelled.
"""
import argparse, collections, csv, json, os, re

SASIA = {"india", "sri lanka", "bangladesh", "nepal", "bhutan", "maldives", "south asia"}
LATAM = {"brazil", "mexico", "colombia", "puerto rico", "argentina", "ecuador", "chile", "peru", "honduras",
         "costa rica", "dominican republic", "venezuela", "guatemala", "el salvador", "uruguay", "cuba", "haiti",
         "guyana", "suriname", "bolivia", "paraguay", "panama", "nicaragua", "latin america", "south america",
         "central america", "caribbean", "jamaica", "trinidad and tobago", "french guiana", "martinique", "guadeloupe"}
AFR = {"south africa", "somalia", "mali", "senegal", "ethiopia", "nigeria", "tanzania", "cameroon", "ghana", "congo",
       "democratic republic of the congo", "africa", "sub-saharan africa", "gambia", "cape verde", "kenya", "uganda",
       "rwanda", "burundi", "zambia", "zimbabwe", "malawi", "mozambique", "angola", "namibia", "botswana", "benin",
       "togo", "burkina faso", "niger", "chad", "guinea", "sierra leone", "liberia", "ivory coast", "cote d'ivoire",
       "côte d'ivoire", "madagascar", "eritrea", "djibouti", "gabon", "lesotho", "eswatini", "swaziland", "mauritania"}
WORDS = {"SASIA": r"\b(indian|south[ -]asian|sri lankan|bangladeshi|nepali|nepalese|tamil|gujarati|bengali)\b",
         "LATAM": r"\b(hispanic|latino|latina|latin american|mexican|brazilian|colombian|puerto rican|argentin\w*|chilean|"
                  r"peruvian|ecuadorian|cuban|venezuelan|dominican|guatemalan|mestizo|amerindian)\b",
         "AFR": r"\b(african|black|nigerian|ghanaian|somali|ethiopian|kenyan|senegalese|malian|cameroonian|yoruba|igbo|bantu|xhosa|zulu)\b"}
SETS = {"SASIA": SASIA, "LATAM": LATAM, "AFR": AFR}


def classify(text):
    t = (text or "").strip().lower()
    if not t or t == "null": return None
    for r, s in SETS.items():
        if t in s: return r
    hits = {r for r, rx in WORDS.items() if re.search(rx, t)}
    hits |= {r for r, s in SETS.items() if any(re.search(r"\b" + re.escape(c) + r"\b", t) for c in s)}
    if "AFR" in hits and re.search(r"\bnorth african\b", t) and not re.search(r"\b(sub-saharan|black)\b", t):
        hits.discard("AFR")
    return hits.pop() if len(hits) == 1 else None


def val(v):
    if isinstance(v, dict): return v.get("label") or v.get("text") or ""
    return ""


def fine_label(ext):
    order = ["geographicOrigin", "ancestry"]
    if (val(ext.get("populationRegionSource")) or "") != "patient": order.append("placeOfCare")
    order += ["publicationPatientOrigins", "recruitmentCountries"]
    for k in order:
        v = val(ext.get(k))
        if k in ("publicationPatientOrigins", "recruitmentCountries") and v.startswith("["):
            try:
                L = json.loads(v)
                labs = {classify(o.get("country") or "") or classify(o.get("ethnicity_or_ancestry") or "") for o in L}
                labs.discard(None)
                if len(labs) == 1: return labs.pop()
            except (ValueError, AttributeError):
                pass
            continue
        r = classify(v)
        if r: return r
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--packets", required=True); ap.add_argument("--cases", action="append", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--author-country", action="append", default=[], help="TSV or JSONL with pmid, first_country")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    idx = {r["id"]: r for r in csv.DictReader(open(os.path.join(a.packets, "index.tsv")), delimiter="\t")}
    fc = {}
    for path in a.author_country:
        if path.endswith(".jsonl"):
            for l in open(path):
                d = json.loads(l); fc.setdefault(d["pmid"], d.get("first_country") or "")
        else:
            for r in csv.DictReader(open(path), delimiter="\t"): fc.setdefault(r["pmid"], r.get("first_country") or "")
    st = collections.Counter(); sa = collections.Counter()
    for path in a.cases:
        rows = list(csv.DictReader(open(path), delimiter="\t"))
        for row in rows:
            if row.get("author_region") == "OTHER":
                r = classify(fc.get(row["pmid"], ""))
                if r: row["author_region"] = r
                sa[(row["corpus"], row["author_region"])] += 1
            if row["region"] != "OTHER": continue
            pp = json.load(open(os.path.join(a.packets, idx[row["id"]]["file"])))
            ext = {e["name"]: e["value"] for e in pp["subject"].get("extensions", [])}
            r = fine_label(ext)
            if r: row["region"] = r
            st[(row["corpus"], row["region"])] += 1
        with open(os.path.join(a.out_dir, os.path.basename(path)), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
    for k, n in sorted(st.items()): print("population", *k, n, sep="\t")
    for k, n in sorted(sa.items()): print("author", *k, n, sep="\t")


if __name__ == "__main__":
    main()
