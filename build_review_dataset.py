#!/usr/bin/env python3
"""Assemble the independent-review dataset for archiving (results/v5/review_dataset).

Joins the two reviewers' returned sheets with the document key and the publication metadata, adds the pipeline
output for the same documents, and writes tab-separated tables with a Frictionless data package descriptor
(datapackage.json), a README and the reviewer protocol.  Source documents are not copied: each document is
identified by PMID, PMCID and DOI, and only the phrases the reviewers cited are included.

  python3 build_review_dataset.py --review data/v5_private/human_review --sources data/v5_private/paired_sources \
      --mapping data/v5_private/paired_mapping_final --termclass results/termclass.jsonl --obo data/hpo/hp.obo \
      --out results/v5/review_dataset
"""
import argparse, csv, hashlib, json, re, zipfile
from pathlib import Path
from rank_engine import parse_obo

REV = {"PNS": ("Paul", "Paul N. Schofield"), "MA": ("Marwa", "Marwa Abdelhakim")}
BASIS = ["history", "examination", "laboratory", "imaging", "physiological_test", "pathology", "mixed", "unclear"]
STATUS = ["present", "explicitly_absent", "uncertain"]
REGIONS = ["EUR", "MENA", "SASIA", "LATAM", "AFR", "EASIA"]


def clean(v):
    return re.sub(r"\s+", " ", (v or "").replace("\r", " ")).strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("review", "sources", "mapping", "termclass", "obo", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(); R = Path(a.review); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    terms, alt = parse_obo(a.obo); ver = next(l.split("releases/")[1].strip() for l in open(a.obo) if l.startswith("data-version"))
    names = {}; cur = None
    for l in open(a.obo):
        if l.startswith("id: "): cur = l[4:].strip()
        elif l.startswith("name: ") and cur: names.setdefault(cur, l[6:].strip())
    tc = {}
    for f in (Path(a.mapping) / "termclass_extra.jsonl", Path(a.termclass)):
        for x in map(json.loads, open(f)):
            if x.get("class"): tc[x["id"]] = x["class"]
    key = {r["document_id"]: r for r in csv.DictReader(open(R / "PRIVATE_alignment_key.csv"))}

    def hpo(h):
        h = alt.get(h, h)
        return (h, names.get(h, "")) if h in terms and not terms[h]["obsolete"] else ("", "")
    tables = {}
    # documents
    rows = []
    for d in sorted(key):
        k = key[d]; m = json.load(open(Path(a.sources) / f"{k['pmid']}.metadata.json"))
        rows.append({"document_id": d, "source": k["source"], "pmid": k["pmid"], "pmcid": m.get("pmcid", ""), "doi": m.get("doi", ""),
                     "publication_year": m.get("pubYear", ""), "title": clean(m.get("title")), "licence": (m.get("license") or "").upper(),
                     "target_gene": k["gene"], "regional_group": k["region"]})
    tables["documents"] = rows
    # reviewer sheets
    F, P, Q, A = [], [], [], []
    for code, (folder, _) in REV.items():
        f = R / "returned" / folder / f"findings_{code}.csv"
        for n, r in enumerate(csv.DictReader(open(f, encoding="utf-8-sig")), 1):
            r = {k: clean(v) for k, v in r.items() if k}
            if not any(r.values()): continue
            b = [x.strip().replace(" ", "_") for x in re.split(r"[,;/|]", r["observation_basis"].lower()) if x.strip()]
            known = sorted({x for x in b if x != "unclear"}); basis = known[0] if len(known) == 1 else ("mixed" if known else ("unclear" if b else ""))
            ids = [i.replace("_", ":") for i in re.findall(r"HP[:_]\d{7}", r["hpo_id"])] or [""]
            for i in ids:
                h, lab = hpo(i)
                F.append({"reviewer": code, "row": n, "document_id": r["document_id"].upper(), "patient_label": r["patient_label"], "source_phrase": r["source_phrase"],
                          "status": r["status"].lower().replace(" ", "_"), "hpo_id_as_entered": r["hpo_id"], "hpo_label_as_entered": r["hpo_label"], "hpo_id": h, "hpo_label": lab,
                          "observation_basis_as_entered": r["observation_basis"], "observation_basis": basis if basis in BASIS else "", "source_locator": r["source_locator"], "notes": r.get("notes", "")})
        f = R / "returned" / folder / f"patients_{code}.csv"
        for n, r in enumerate(csv.DictReader(open(f, encoding="utf-8-sig")), 1):
            r = {k: clean(v) for k, v in r.items() if k and not k.startswith("Column")}
            if any(r.values()): P.append({"reviewer": code, "row": n, **r, "document_id": r["document_id"].upper()})
        f = R / "returned" / folder / f"frequencies_{code}.csv"
        if f.exists():
            for n, r in enumerate(csv.DictReader(open(f, encoding="utf-8-sig")), 1):
                r = {k: clean(v) for k, v in r.items() if k and not k.startswith("Column")}
                if any(r.values()):
                    h, lab = hpo(r["hpo_id"].replace("_", ":")); Q.append({"reviewer": code, "row": n, **r, "document_id": r["document_id"].upper(), "hpo_label": lab if h else ""})
        f = R / "returned" / folder / f"documents_{code}.csv"
        if f.exists():
            for r in csv.DictReader(open(f, encoding="utf-8-sig")):
                r = {k: clean(v) for k, v in r.items() if k}
                if any(r.values()): A.append({"reviewer": code, **r})
    tables["findings"], tables["patients"], tables["frequencies"], tables["document_assessments"] = F, P, Q, A
    # pipeline output for the same documents
    sel = {x["key"]: x for x in map(json.loads, open(Path(a.mapping) / "selected.jsonl"))}
    phr = {x["key"]: x for x in map(json.loads, open(Path(a.mapping) / "phrases.jsonl"))} if (Path(a.mapping) / "phrases.jsonl").exists() else {}
    doc_of = {f"{k['pmid']}:{k['gene']}:{k['source']}": d for d, k in key.items()}; G = []
    for m in map(json.loads, open(Path(a.mapping) / "mentions.jsonl")):
        d = doc_of.get(m["document_key"])
        if not d: continue
        s = sel.get(m["phrase_key"], {}); h, lab = hpo(s.get("hpo") or "")
        G.append({"document_id": d, "patient_label": clean(str(m.get("patient_label") or "")), "unit": m.get("unit") or "", "status": m["status"],
                  "phrase": clean(phr.get(m["phrase_key"], {}).get("text") or phr.get(m["phrase_key"], {}).get("phrase") or ""), "evidence": clean(m.get("evidence")),
                  "hpo_id": h, "hpo_label": lab, "finding_class": tc.get(h, "")})
    tables["pipeline_findings"] = G
    C = [{"hpo_id": h, "hpo_label": names.get(h, ""), "finding_class": tc[h]} for h in sorted({r["hpo_id"] for r in F + G if r["hpo_id"] and r["hpo_id"] in tc})]
    tables["finding_classes"] = C

    desc = {
     "documents": {"document_id": "Neutral identifier shown to the reviewers (D001 to D040).", "source": "Text presented: abstract or fulltext of the publication.", "pmid": "PubMed identifier of the publication.", "pmcid": "PubMed Central identifier.", "doi": "Digital object identifier.", "publication_year": "Year of publication.", "title": "Title of the publication.", "licence": "Licence of the publication as reported by Europe PMC; applies to quoted phrases.", "target_gene": "HGNC symbol of the causative gene for which the publication was sampled.", "regional_group": "Regional group assigned to the publication by the study (withheld from reviewers)."},
     "findings": {"reviewer": "Reviewer code (PNS or MA).", "row": "Row number in the reviewer's returned sheet.", "document_id": "Document reviewed.", "patient_label": "Label of the patient or aggregate as written by the reviewer.", "source_phrase": "Phrase or sentence in the document that supports the finding.", "status": "present, explicitly_absent or uncertain; empty if not entered.", "hpo_id_as_entered": "HPO identifier(s) exactly as entered by the reviewer.", "hpo_label_as_entered": "Label as entered by the reviewer.", "hpo_id": "Identifier resolved in the HPO release named in the README (one row per identifier; empty if absent from that release).", "hpo_label": "Primary label of hpo_id in that release.", "observation_basis_as_entered": "How the finding was established in this report, as entered.", "observation_basis": "Normalised basis: a single stated value, mixed if several, unclear, or empty.", "source_locator": "Section or location in the document.", "notes": "Reviewer's notes."},
     "patients": {"reviewer": "Reviewer code.", "row": "Row number in the returned sheet.", "document_id": "Document reviewed.", "patient_label": "Label of the patient or aggregate.", "unit": "individual or aggregate.", "number_of_patients": "Number of patients the row describes.", "diagnosis": "Diagnosis as stated.", "causal_gene": "Causal gene as entered.", "variant": "Variant as reported.", "allelic_state": "Allelic state as entered.", "ancestry": "Stated ancestry or origin.", "care_country": "Country in which the patient was cared for.", "care_evidence": "Text supporting the care country.", "identity_evidence": "Text identifying the patient.", "notes": "Reviewer's notes."},
     "frequencies": {"reviewer": "Reviewer code.", "row": "Row number in the returned sheet.", "document_id": "Document reviewed.", "phenotype": "Phenotype as stated.", "hpo_id": "HPO identifier as entered.", "hpo_label": "Primary label in the HPO release, if the identifier resolves.", "affected": "Number of patients with the phenotype, as entered.", "assessed": "Number of patients assessed, as entered.", "cohort_definition": "Cohort to which the frequency refers.", "age_range": "Age range, as entered.", "source_locator": "Location in the document.", "overlap_notes": "Notes on overlap with other cohorts.", "notes": "Reviewer's notes."},
     "document_assessments": {"reviewer": "Reviewer code.", "document_id": "Document reviewed.", "target_gene": "Gene as entered.", "assessment_complete": "Whether the reviewer completed the document.", "source_incomplete": "Whether the reviewer judged the supplied text incomplete.", "notes": "Reviewer's notes."},
     "pipeline_findings": {"document_id": "Document processed.", "patient_label": "Patient label in the model output.", "unit": "individual or aggregate.", "status": "present or absent (explicitly excluded), as returned by the model.", "phrase": "Extracted phrase.", "evidence": "Text the model cited as evidence.", "hpo_id": "HPO class selected by the mapping procedure (empty if unmapped).", "hpo_label": "Primary label of hpo_id.", "finding_class": "Finding class of hpo_id used in the study."},
     "finding_classes": {"hpo_id": "HPO class occurring in the findings or pipeline_findings table.", "hpo_label": "Primary label.", "finding_class": "Class assigned once per HPO class: symptom_history, clinical_exam, laboratory, imaging, physiological_test, pathology or mixed."}}
    enum = {"source": ["abstract", "fulltext"], "regional_group": REGIONS, "reviewer": list(REV)}
    res = []
    for name, rows in tables.items():
        cols = list(desc[name]); p = out / f"{name}.tsv"
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n"); w.writeheader()
            for r in rows: w.writerow({c: r.get(c, "") for c in cols})
        fields = []
        for c in cols:
            fd = {"name": c, "type": "integer" if c == "row" else "string", "description": desc[name][c]}
            if c in enum: fd["constraints"] = {"enum": enum[c]}
            if c == "hpo_id" and name != "frequencies": fd["rdfType"] = "http://purl.obolibrary.org/obo/HP_0000001"; fd["constraints"] = {"pattern": "(HP:[0-9]{7})?"}
            fields.append(fd)
        sch = {"fields": fields, "missingValues": [""]}
        if name != "documents" and "document_id" in cols: sch["foreignKeys"] = [{"fields": "document_id", "reference": {"resource": "documents", "fields": "document_id"}}]
        if name == "documents": sch["primaryKey"] = "document_id"
        res.append({"name": name.replace("_", "-"), "path": p.name, "format": "tsv", "mediatype": "text/tab-separated-values", "encoding": "utf-8", "dialect": {"delimiter": "\t"},
                    "rows": len(rows), "hash": "sha256:" + hashlib.sha256(p.read_bytes()).hexdigest(), "schema": sch})
        print(name, len(rows))
    z = zipfile.ZipFile(R / "Paul_source_validation.zip"); (out / "protocol.txt").write_bytes(z.read("Paul/protocol.txt"))
    pkg = {"name": "independent-review-of-phenotype-extraction", "title": "Independent review of phenotype extraction from 20 rare disease case reports",
           "description": "Two reviewers independently recorded patients, phenotypic findings (HPO), observation basis and frequency statements in the abstract and full text of 20 case reports, blind to model output and regional assignment; with the automated pipeline output for the same documents.",
           "version": "1.0.0", "licenses": [{"name": "CC-BY-4.0", "path": "https://creativecommons.org/licenses/by/4.0/", "title": "Creative Commons Attribution 4.0"}],
           "contributors": [{"title": "Paul N. Schofield", "role": "author", "path": "https://orcid.org/0000-0002-5111-7263", "organization": "Department of Physiology, Development and Neuroscience, University of Cambridge, United Kingdom"},
                            {"title": "Marwa Abdelhakim", "role": "author", "path": "https://orcid.org/0000-0001-6816-2119", "organization": "King Abdullah University of Science and Technology, Thuwal, Saudi Arabia"}],
           "keywords": ["Human Phenotype Ontology", "rare disease", "case reports", "phenotype extraction", "inter-annotator agreement", "gold standard"],
           "sources": [{"title": "Human Phenotype Ontology release " + ver, "path": "http://purl.obolibrary.org/obo/hp/releases/" + ver + "/hp.obo"}, {"title": "Europe PMC", "path": "https://europepmc.org"}],
           "homepage": "https://github.com/bio-ontology-research-group/population-phenotype-differences", "resources": res}
    json.dump(pkg, open(out / "datapackage.json", "w"), indent=1)
    json.dump({"hpo_release": ver, "counts": {k: len(v) for k, v in tables.items()}}, open(out / "summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()
