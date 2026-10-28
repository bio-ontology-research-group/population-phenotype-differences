#!/usr/bin/env python3
"""Add age_band and sex columns to case tables from the phenopacket index.
Age is the age at report (ISO 8601 duration) or, when only an HPO onset class is known, that class.
Bands: infant (<1 y, incl. neonatal), child (1-11), adolescent (12-17), adult (>=18); prenatal and unknown left empty.

  add_age_sex.py --index results/phenopackets/index.tsv --cases A.tsv [...] --out-dir D
"""
import argparse, csv, os, re

ONSET = {"Neonatal onset": "infant", "Infantile onset": "infant", "Congenital onset": "infant",
         "Childhood onset": "child", "Juvenile onset": "adolescent", "Adult onset": "adult",
         "Young adult onset": "adult", "Middle age onset": "adult", "Late onset": "adult"}


def years(iso):
    m = re.fullmatch(r"P(?:(\d+)Y)?(?:(\d+)M)?(?:(\d+)W)?(?:(\d+)D)?(?:T(\d+)H)?", iso or "")
    if not m or not any(m.groups()): return None
    y, mo, w, d, h = (int(x) if x else 0 for x in m.groups())
    return y + mo / 12 + w / 52.18 + d / 365.25 + h / 8766


def band(age):
    y = years(age)
    if y is None: return ONSET.get(age, "")
    return "infant" if y < 1 else "child" if y < 12 else "adolescent" if y < 18 else "adult"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True); ap.add_argument("--cases", action="append", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    idx = {r["id"]: r for r in csv.DictReader(open(a.index), delimiter="\t")}
    for path in a.cases:
        rows = list(csv.DictReader(open(path), delimiter="\t"))
        for r in rows:
            x = idx.get(r["id"], {})
            r["age_band"] = band(x.get("age", ""))
            r["sex"] = {"MALE": "male", "FEMALE": "female"}.get(x.get("sex", ""), "")
        fields = [k for k in rows[0] if k not in ("age_band", "sex")] + ["age_band", "sex"]
        with open(os.path.join(a.out_dir, os.path.basename(path)), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
