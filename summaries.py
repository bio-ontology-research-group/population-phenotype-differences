#!/usr/bin/env python3
"""Descriptive numbers quoted in RESULTS.md and the paper that have no other result file.

Writes results/v4/summaries.txt: cases per corpus and population (fine2 tables and the PAVS table), age-band and sex
composition per population, finding classes per HPO class (termclass.jsonl), and the annotation-free gene recovery
Hits@1 per corpus and population (recovery/*_percase.csv; a case counts as a hit when the causative gene ranks first).

  summaries.py [--res results/v4] [--termclass results/termclass.jsonl]
"""
import argparse, collections, csv, glob, json, os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", default="results/v4"); ap.add_argument("--termclass", default="results/termclass.jsonl")
    a = ap.parse_args()
    out = []
    P = ["EUR", "MENA", "SASIA", "LATAM", "AFR", "EASIA"]
    cases = []
    for f in ("cases_store.tsv", "cases_pmcp_map2.tsv", "cases_abs.tsv"):
        cases += list(csv.DictReader(open(os.path.join(a.res, "fine2", f)), delimiter="\t"))
    pavs = list(csv.DictReader(open(os.path.join(a.res, "fine", "cases_pavs.tsv")), delimiter="\t"))
    out.append("# Cases by corpus and population (fine2 tables; PAVS from fine/cases_pavs.tsv)")
    tab = collections.Counter((c["corpus"], c["region"]) for c in cases + pavs)
    corp = sorted({c for c, _ in tab})
    out.append("corpus\t" + "\t".join(P) + "\tother\ttotal")
    for c in corp:
        row = [tab[(c, p)] for p in P]; oth = sum(v for (cc, r), v in tab.items() if cc == c and r not in P)
        out.append(c + "\t" + "\t".join(map(str, row)) + f"\t{oth}\t{sum(row) + oth}")
    out.append(f"total cases: {len(cases) + len(pavs)}")
    out.append("\n# Age band and sex per population (fine2; share of cases with the attribute stated)")
    for var, levels in (("age_band", ["infant", "child", "adolescent", "adult"]), ("sex", ["male", "female"])):
        out.append(f"{var}\t" + "\t".join(levels) + "\tstated\tcases")
        for p in P:
            sub = [c for c in cases if c["region"] == p]; st = [c for c in sub if c.get(var) in levels]
            shares = [f"{100 * sum(c[var] == l for c in st) / max(len(st), 1):.1f}" for l in levels]
            out.append(p + "\t" + "\t".join(shares) + f"\t{len(st)}\t{len(sub)}")
    out.append("\n# Finding classes per HPO class (termclass.jsonl)")
    tc = collections.Counter(x["class"] for x in map(json.loads, open(a.termclass)))
    for k, v in sorted(tc.items(), key=lambda kv: -kv[1]): out.append(f"{k}\t{v}")
    test = sum(v for k, v in tc.items() if k in ("laboratory", "imaging", "physiological_test", "pathology"))
    out.append(f"test-based {test}; bedside {tc['clinical_exam'] + tc['symptom_history']} (signs {tc['clinical_exam']}, symptoms {tc['symptom_history']}); total {sum(tc.values())}")
    out.append("\n# Annotation-free gene recovery, Hits@1 by corpus and population (recovery/*_percase.csv)")
    hits = collections.defaultdict(lambda: [0, 0]); seen = set()
    for f in sorted(glob.glob(os.path.join(a.res, "recovery", "recovery_*_percase.csv"))):
        for r in csv.DictReader(open(f)):
            key = (r["id"], r["corpus"])
            if key in seen: continue
            seen.add(key)
            h = hits[(r["corpus"], r["region"])]; h[1] += 1; h[0] += float(r["rank"]) <= 1
    out.append("corpus\tpopulation\tqueries\thits@1 (%)")
    for (c, p), (h, n) in sorted(hits.items()): out.append(f"{c}\t{p}\t{n}\t{100 * h / n:.1f}")
    path = os.path.join(a.res, "summaries.txt")
    open(path, "w").write("\n".join(out) + "\n"); print("\n".join(out)); print("wrote", path)


if __name__ == "__main__":
    main()
