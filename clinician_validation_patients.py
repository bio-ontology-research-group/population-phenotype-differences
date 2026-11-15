#!/usr/bin/env python3
"""Score the pipeline's patient records against a reviewer's patients.csv (prepare_human_review.py packets).

Per source document, the reviewer's rows (unit, causal gene, allelic state, ancestry, care country) are compared
with the common-pipeline extraction of the same document (data/v5_private/paired_extraction.jsonl).  Ancestry and
care entries are mapped to regional groups with setting_labels.region, the rule used for the case tables
(ancestry region from origin country or ancestry text, care region from the place of care).  Reported:

  units      individual patients per document, and aggregate rows on either side
  gene       reviewer's causal gene against the target gene of the document
  zygosity   allelic state of matched individuals (hom, het, comp_het, hemi), where both sides state one
  ancestry   regional group of the ancestry per document: same, stated by one side only, or conflicting
  care       regional group of the place of care per document, likewise
  region     per publication, the reviewer's ancestry region (else care region) against the assigned region

  python3 clinician_validation_patients.py --patients data/v5_private/human_review/returned/Paul/patients_PNS.csv \
      --key data/v5_private/human_review/PRIVATE_alignment_key.csv \
      --extraction data/v5_private/paired_extraction.jsonl --out results/v5/clinician_validation/PNS_patients
"""
import argparse, collections, csv, json, re
from setting_labels import region


def zyg(t):
    t = (t or "").lower()
    if "compound" in t: return "comp_het"
    if "hemi" in t: return "hemi"
    if "homo" in t: return "hom"
    if "hetero" in t: return "het"
    return None


def norm(t):
    return re.sub(r"[^a-z0-9]", "", (t or "").lower())


def one(xs):
    s = {x for x in xs if x}
    return "/".join(sorted(s)) if s else None


def cmp(a, b):
    if a == b: return "same" if a else "both_unstated"
    if a and b: return "conflict"
    return "reviewer_only" if a else "pipeline_only"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("patients", "key", "extraction", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    key = {r["document_id"]: r for r in csv.DictReader(open(a.key))}
    ex = {}
    for l in open(a.extraction):
        d = json.loads(l)
        if d.get("result"): ex[(str(d["pmid"]), d["gene"], d["source"])] = d["result"]
    rev = collections.defaultdict(list)
    for r in csv.DictReader(open(a.patients, encoding="utf-8-sig")): rev[r["document_id"].strip()].append(r)
    rows, zy = [], collections.Counter()
    for d in sorted(key):
        k = key[d]; pp = ex.get((k["pmid"], k["gene"], k["source"]), {}).get("patients", [])
        ri = [r for r in rev[d] if r["unit"].strip() == "individual"]; pi = [p for p in pp if p.get("unit") == "individual"]
        # match individuals by label; a single individual on each side is the same patient
        pairs = [(ri[0], pi[0])] if len(ri) == len(pi) == 1 else [(r, p) for r in ri for p in pi if norm(r["patient_label"]) and norm(r["patient_label"]) in norm(p.get("local_id"))]
        for r, p in pairs:
            x, y = zyg(r["allelic_state"]), zyg(p.get("zygosity")); zy[cmp(x, y)] += 1
        ra, rc = one(region(r["ancestry"]) for r in ri), one(region(r["care_country"]) for r in ri)
        pa, pc = one(region(p.get("origin_country")) or region(p.get("ancestry")) for p in pi), one(region(p.get("care_country")) for p in pi)
        rows.append({"document_id": d, "pmid": k["pmid"], "source": k["source"], "assigned_region": k["region"],
                     "individuals_reviewer": len(ri), "individuals_pipeline": len(pi), "matched": len(pairs),
                     "aggregates_reviewer": len(rev[d]) - len(ri), "aggregates_pipeline": len(pp) - len(pi),
                     "gene_ok": all(norm(k["gene"]) in norm(r["causal_gene"]) for r in rev[d]) and bool(rev[d]),
                     "ancestry_reviewer": ra or "", "ancestry_pipeline": pa or "", "ancestry": cmp(ra, pa),
                     "care_reviewer": rc or "", "care_pipeline": pc or "", "care": cmp(rc, pc)})
    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    pub = collections.defaultdict(list)
    for r in rows: pub[r["pmid"]].append(r)
    reg = collections.Counter()
    for v in pub.values():
        got = one(r["ancestry_reviewer"] for r in v) or one(r["care_reviewer"] for r in v)
        reg["unstated" if not got else "same" if got == v[0]["assigned_region"] else "conflict"] += 1
    rep = {"documents": len(rows), "reviewer_rows": sum(len(v) for v in rev.values()),
           "individuals": {"reviewer": sum(r["individuals_reviewer"] for r in rows), "pipeline": sum(r["individuals_pipeline"] for r in rows),
                           "documents_same_count": sum(r["individuals_reviewer"] == r["individuals_pipeline"] for r in rows), "matched": sum(r["matched"] for r in rows)},
           "aggregates": {"reviewer": sum(r["aggregates_reviewer"] for r in rows), "pipeline": sum(r["aggregates_pipeline"] for r in rows)},
           "gene_ok_documents": sum(r["gene_ok"] for r in rows), "zygosity_matched_individuals": dict(zy),
           "ancestry_region_documents": dict(collections.Counter(r["ancestry"] for r in rows)),
           "care_region_documents": dict(collections.Counter(r["care"] for r in rows)),
           "publication_region_vs_assigned": dict(reg)}
    json.dump(rep, open(a.out + ".json", "w"), indent=1); print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
