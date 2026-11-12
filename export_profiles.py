#!/usr/bin/env python3
"""Export the three HPO profiles (all, bedside-only, random same-size subset) of every paired case, for ranking with
external tools (tool_restriction).

Uses the code path of `gda_impact.py rank --random-control`: same case files (read_cases), same term mapping
(Engine.map_terms over hp.obo), same restriction to cases whose gene is in genes_to_phenotype.txt, same bedside/test
definition (termclass.jsonl, BED/TEST), same seeded random draw. Writes one JSONL line per case and validates the
counts against the existing random-control output (gdarand_chunk*.csv): n_terms, n_bed, n_rand_test must match.

  python3 export_profiles.py --cases A.tsv [...] --hpo-dir hpo --termclass T --rand 'gda_rand/gdarand_chunk*.csv' --out P.jsonl
"""
import argparse, csv, glob, hashlib, json, os, sys
import numpy as np

from gda_impact import read_cases, BED, TEST
from rank_engine import Engine


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", action="append", required=True); ap.add_argument("--hpo-dir", required=True)
    ap.add_argument("--termclass", required=True); ap.add_argument("--rand", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    eng = Engine(os.path.join(a.hpo_dir, "hp.obo"), os.path.join(a.hpo_dir, "genes_to_phenotype.txt"))
    tc = {x["id"]: x["class"] for x in map(json.loads, open(a.termclass))}
    ref = {r["id"]: r for f in sorted(glob.glob(a.rand)) for r in csv.DictReader(open(f))}
    cases = [r for r in read_cases(a.cases) if r["gene"] in eng.gid]
    n_out = mism = 0; seen = set()
    with open(a.out, "w") as fo:
        for r in cases:
            if r["id"] not in ref: continue
            terms = eng.map_terms(r["terms"].split(";"))
            if not terms: continue
            bed = [t for t in terms if tc.get(t) in BED]; tst = [t for t in terms if tc.get(t) in TEST]
            if not bed or not tst: continue
            seed = int(hashlib.sha1(r["id"].encode()).hexdigest()[:8], 16)
            keep = sorted(np.random.default_rng(seed).choice(len(terms), size=len(bed), replace=False))
            rnd = [terms[k] for k in keep]
            g = ref[r["id"]]
            ok = (int(g["n_terms"]) == len(terms) and int(g["n_bed"]) == len(bed) and int(g["n_test"]) == len(tst)
                  and int(g["n_rand_test"]) == sum(tc.get(t) in TEST for t in rnd) and int(g["seed"]) == seed
                  and g["gene"] == r["gene"] and g["region"] == r["region"])
            if not ok:
                mism += 1
                if mism <= 10: print("MISMATCH", r["id"], len(terms), len(bed), len(tst), g, file=sys.stderr)
            if r["id"] in seen: print("DUPLICATE", r["id"], file=sys.stderr)
            seen.add(r["id"])
            fo.write(json.dumps({"id": r["id"], "gene": r["gene"], "region": r["region"], "corpus": r["corpus"],
                                 "n_terms": len(terms), "n_bed": len(bed), "n_test": len(tst),
                                 "terms_all": terms, "terms_bed": bed, "terms_rand": rnd}) + "\n")
            n_out += 1
    missing = set(ref) - seen
    print(f"exported {n_out} cases; reference {len(ref)}; reference ids not exported {len(missing)}; count mismatches {mism}")
    if mism or missing: sys.exit(1)


if __name__ == "__main__":
    main()
