#!/usr/bin/env python3
"""Write every analysed patient (Store, PMC-Patients, PubMed abstracts, PAVS) as a GA4GH Phenopacket v2 JSON.

Core v2 fields: subject (sex, timeAtLastEncounter), phenotypicFeatures (HPO; verbatim phrase in
`description`; excluded where extracted), interpretations (causative gene, allelic state when
reported), metaData (resources, source PMID). Store and PAVS packets are copied unchanged and
only augmented. Fields that Phenopacket v2 has no slot for (origin, place of care, migrant
status, onset as reported, analysis labels, provenance) go into `subject.extensions` as
name/value pairs, the convention PAVS already uses (geographicOrigin, consanguinity).
Countries are GAZ terms when the name maps to a GAZ country (OLS-verified table), text otherwise.

Usage: to_phenopackets.py --sp SCRATCH --gaz gaz_map.json --out results/phenopackets
"""
import argparse, collections, copy, csv, datetime, glob, hashlib, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from abs_to_cases import parse_age
from origin_analyses import region_of

TODAY = datetime.date.today().isoformat()
ONSET = {"prenatal": ("HP:0030674", "Antenatal onset"), "neonate": ("HP:0003623", "Neonatal onset"),
         "infant": ("HP:0003593", "Infantile onset"), "child": ("HP:0011463", "Childhood onset"),
         "adolescent": ("HP:0003621", "Juvenile onset"), "adult": ("HP:0003581", "Adult onset")}
GENO = {"heterozygous": ("GENO:0000135", "heterozygous"), "homozygous": ("GENO:0000136", "homozygous"),
        "hemizygous": ("GENO:0000134", "hemizygous"), "compound_heterozygous": ("GENO:0000402", "compound heterozygous")}
CARE_SOURCE = {"A": "Saudi Arabia", "B": "Saudi Arabia", "F": "Saudi Arabia", "P": "Saudi Arabia", "Q": "USA"}
RES = [{"id": "hp", "name": "Human Phenotype Ontology", "url": "http://purl.obolibrary.org/obo/hp.owl", "version": "2026-01-08",
        "namespacePrefix": "HP", "iriPrefix": "http://purl.obolibrary.org/obo/HP_"},
       {"id": "hgnc", "name": "HUGO Gene Nomenclature Committee", "url": "https://www.genenames.org", "version": "",
        "namespacePrefix": "HGNC", "iriPrefix": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/"},
       {"id": "geno", "name": "Genotype Ontology", "url": "http://purl.obolibrary.org/obo/geno.owl", "version": "",
        "namespacePrefix": "GENO", "iriPrefix": "http://purl.obolibrary.org/obo/GENO_"},
       {"id": "gaz", "name": "Gazetteer", "url": "http://purl.obolibrary.org/obo/gaz.obo", "version": "",
        "namespacePrefix": "GAZ", "iriPrefix": "http://purl.obolibrary.org/obo/GAZ_"}]
CORPUS_DESC = {"store": "Phenopacket Store 0.1.26 (curated); origin labels extracted from the publication by Qwen3.8-27B",
               "pmcp": "PMC-Patients V2 patient summary, extracted by Qwen3.8-27B, HPO by RAG-style candidate selection",
               "pubmed_abs": "PubMed abstract, extracted by Qwen3.8-27B, HPO by RAG-style candidate selection",
               "pavs_clin": "PAVS Saudi clinical records", "pavs_lit": "PAVS Saudi literature curation",
               "pavs_us": "PAVS US clinical records"}


def ext(name, value):
    return {"name": name, "value": value}


def oc(i, l):
    return {"id": i, "label": l}


def load_hpo(obo):
    lab, alt, cur = {}, {}, None
    for line in open(obo):
        line = line.rstrip("\n")
        if line == "[Term]": cur = {}
        elif cur is not None:
            if line.startswith("id: "): cur["id"] = line[4:]
            elif line.startswith("name: "): lab[cur["id"]] = line[6:]
            elif line.startswith("alt_id: "): alt[line[8:]] = cur["id"]
            elif line == "": cur = None
    return lab, alt


def time_element(s):
    y, b = parse_age(s)
    if y is not None:
        if y >= 0.25:
            mo = round(y * 12); yy, mm = divmod(mo, 12)
            dur = "P" + (f"{yy}Y" if yy else "") + (f"{mm}M" if mm or not yy else "")
        elif y >= 2 / 365.25:
            dur = f"P{max(1, round(y * 365.25))}D"
        else:
            dur = f"PT{max(1, round(y * 8766))}H"
        return {"age": {"iso8601duration": dur}}
    if b: return {"ontologyClass": oc(*ONSET[b])}
    return None


def mkkey(t, l):
    return hashlib.md5(((t or "").strip().lower() + "||" + (l or "").strip().lower()).encode()).hexdigest()[:16]


class Builder:
    def __init__(self, a):
        self.lab, self.alt = load_hpo(os.path.join(a.sp, "hpo", "hp.obo"))
        self.gaz = json.load(open(a.gaz))
        self.hgnc = {}
        for r in csv.DictReader(open(os.path.join(a.sp, "hgnc", "hgnc_complete_set.txt")), delimiter="\t"):
            if r["status"] == "Approved": self.hgnc[r["symbol"]] = r["hgnc_id"]
        self.m2 = {}
        for f in glob.glob("results/abs/sel/sel_*.jsonl") + ["results/pmcp/hpo_map2_select.jsonl"]:
            for l in open(f):
                d = json.loads(l)
                if d.get("ok") and d.get("hpo"): self.m2[d["key"]] = d["hpo"]
        self.grp = {r[0]: r[1] for r in csv.reader(open("results/groups_lastauthor.tsv"), delimiter="\t")}
        self.authors = {}
        for l in open("results/all_authors.tsv"):
            p, _, au = l.rstrip("\n").partition("\t"); self.authors[p] = au

    def country(self, name):
        if not name: return None
        g = self.gaz.get(name)
        return oc(g["id"], g["label"]) if g else {"text": name}

    def term(self, h):
        h = self.alt.get(h, h)
        return oc(h, self.lab[h]) if h in self.lab else None

    def features(self, obs, excl=()):
        out, unm, seen = [], [], set()
        for ex, lst in ((False, obs), (True, excl)):
            for q in lst or []:
                if isinstance(q, str): q = {"text": q, "hpo_label": q}
                if not isinstance(q, dict): continue
                h = self.m2.get(mkkey(q.get("text"), q.get("hpo_label")))
                t = self.term(h) if h else None
                if not t: unm.append({"text": q.get("text"), "excluded": ex}); continue
                if (t["id"], ex) in seen: continue
                seen.add((t["id"], ex))
                f = {"type": t, "description": q.get("text") or ""}
                if ex: f["excluded"] = True
                out.append(f)
        return out, unm

    def interp(self, pid, gene, zyg):
        gd = {"valueId": self.hgnc.get(gene, ""), "symbol": gene}
        gi = {"subjectOrBiosampleId": pid, "interpretationStatus": "CAUSATIVE"}
        if zyg in GENO:
            gi["variantInterpretation"] = {"variationDescriptor": {"id": f"{pid}-{gene}", "geneContext": gd,
                                                                  "allelicState": oc(*GENO[zyg])}}
        else:
            gi["gene"] = gd
        return [{"id": f"interp-{pid}", "progressStatus": "SOLVED", "diagnosis": {"genomicInterpretations": [gi]}}]

    def meta(self, pmid, corpus):
        refs = []
        if pmid and pmid.isdigit(): refs.append({"id": f"PMID:{pmid}", "reference": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        refs.append({"id": f"pavs-population:{corpus}", "description": CORPUS_DESC[corpus]})
        return {"created": TODAY + "T00:00:00Z", "createdBy": "pavs-saudi-specific-phenotypes/to_phenopackets.py",
                "resources": RES, "phenopacketSchemaVersion": "2.0", "externalReferences": refs}

    def analysis_ext(self, row):
        e = [ext("populationRegion", {"text": row["region"]}), ext("populationRegionSource", {"text": row.get("region_source", "")}),
             ext("reportingAuthorRegion", {"text": row.get("author_region", "")}), ext("curationCorpus", {"text": row["corpus"]})]
        if row["pmid"] in self.grp: e.append(ext("reportingGroup", {"text": self.grp[row["pmid"]]}))
        if row.get("setting") and row["setting"] != "NA": e.append(ext("ancestryCareSetting", {"text": row["setting"]}))
        return e

    def origin_ext(self, anc, anc_c, care_c, care_basis, migrant, cons, anc_ev=None, care_ev=None):
        e = []
        if anc_c: e.append(ext("geographicOrigin", self.country(anc_c)))
        if anc: e.append(ext("ancestry", {"text": anc}))
        # supporting quotes (anc_ev, care_ev) are article text and are not written to the packets
        if care_c: e.append(ext("placeOfCare", self.country(care_c)))
        if care_basis: e.append(ext("placeOfCareBasis", {"text": care_basis}))
        if migrant: e.append(ext("migrant", {"text": migrant}))
        if cons == "yes": e.append(ext("consanguinity", {"ontologyClass": oc("HP:0025820", self.lab.get("HP:0025820", ""))}))
        elif cons: e.append(ext("consanguinity", {"text": cons}))
        return e


def sex_of(s):
    return {"male": "MALE", "female": "FEMALE"}.get((s or "").lower(), "UNKNOWN_SEX")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sp", required=True); ap.add_argument("--gaz", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--pavs-json", default=os.path.expanduser("~/Public/software/pavs/data/PAVS_phenopackets.json"),
                    help="PAVS_phenopackets.json from the PAVS database export")
    a = ap.parse_args()
    B = Builder(a); st = collections.Counter(); index = []
    def write(corpus, pp, row, extra):
        d = os.path.join(a.out, corpus); os.makedirs(d, exist_ok=True)
        fn = re.sub(r"[^A-Za-z0-9_.-]", "_", pp["id"]) + ".json"
        json.dump(pp, open(os.path.join(d, fn), "w"), indent=1)
        exts = {e["name"]: e["value"] for e in pp["subject"].get("extensions", [])}
        tle = pp["subject"].get("timeAtLastEncounter") or {}
        index.append({"id": pp["id"], "file": f"{corpus}/{fn}", "corpus": row["corpus"], "pmid": row["pmid"], "gene": row["gene"],
                      "region": row["region"], "setting": row.get("setting", ""), "sex": pp["subject"].get("sex", ""),
                      "age": (tle.get("age") or {}).get("iso8601duration") or (tle.get("ontologyClass") or {}).get("label", ""),
                      "origin": (exts.get("geographicOrigin") or {}).get("label") or (exts.get("geographicOrigin") or {}).get("text", ""),
                      "place_of_care": (exts.get("placeOfCare") or {}).get("label") or (exts.get("placeOfCare") or {}).get("text", ""),
                      "n_present": sum(1 for f in pp["phenotypicFeatures"] if not f.get("excluded")),
                      "n_excluded": sum(1 for f in pp["phenotypicFeatures"] if f.get("excluded")), **extra})
        st[corpus] += 1

    # ---- PubMed abstracts
    meta = {}
    for l in open(os.path.join(a.sp, "abs", "abstracts.jsonl")):
        d = json.loads(l); meta[d["pmid"]] = d
    ex = {}
    for f in sorted(glob.glob("results/abs/abs_*.jsonl")):
        for l in open(f):
            d = json.loads(l)
            if d["pmid"] not in ex or d.get("ok"): ex[d["pmid"]] = d
    for row in csv.DictReader(open("results/cases_abs.tsv"), delimiter="\t"):
        pmid, i = row["id"][4:].split("-"); p = ex[pmid]["result"]["patients"][int(i) - 1]
        pid = row["id"]
        zyg = next((g.get("zygosity") for g in p.get("causal_genes") or [] if isinstance(g, dict)), None)
        feats, unm = B.features(p.get("observed_phenotypes"))
        subj = {"id": pid, "sex": sex_of(p.get("sex"))}
        te = time_element(p.get("age"))
        if te: subj["timeAtLastEncounter"] = te
        e = B.origin_ext(p.get("ancestry"), p.get("ancestry_country"), p.get("care_country"), p.get("care_basis"), p.get("migrant"), p.get("consanguinity"))
        if p.get("age"): e.append(ext("ageAsReported", {"text": str(p["age"])}))
        if p.get("age_at_onset"):
            e.append(ext("ageAtOnsetAsReported", {"text": str(p["age_at_onset"])}))
            to = time_element(p["age_at_onset"])
            if to: e.append(ext("ageAtOnset", to))
        if zyg: e.append(ext("zygosityAsReported", {"text": zyg}))
        if unm: e.append(ext("unmappedPhenotypeText", {"text": json.dumps(unm, ensure_ascii=False)}))
        m = meta.get(pmid, {})
        if m.get("countries"): e.append(ext("authorAffiliationCountries", {"text": "; ".join(m["countries"])}))
        subj["extensions"] = e + B.analysis_ext(row)
        pp = {"id": pid, "subject": subj, "phenotypicFeatures": feats, "interpretations": B.interp(pid, row["gene"], zyg),
              "metaData": B.meta(pmid, "pubmed_abs")}
        write("pubmed_abs", pp, row, {})

    # ---- PMC-Patients
    pat = {}
    for l in open("results/pmcp/patients.jsonl"):
        d = json.loads(l)
        if d.get("ok"): pat[d["uid"]] = d["result"]
    sett = {}
    for l in open("results/pmcp/setting.jsonl"):
        d = json.loads(l)
        if d.get("ok") and d.get("result"): sett[d["uid"]] = d["result"]
    for row in csv.DictReader(open("results/cases_pmcp_map2.tsv"), delimiter="\t"):
        uid = row["id"][5:]; p = pat[uid]; s = sett.get(uid, {}); pid = row["id"]
        zyg = next((g.get("zygosity") for g in p.get("causal_genes") or [] if isinstance(g, dict)), None)
        feats, unm = B.features(p.get("observed_phenotypes"), p.get("excluded_phenotypes"))
        subj = {"id": pid, "sex": sex_of(p.get("sex"))}
        te = time_element(p.get("age"))
        if te: subj["timeAtLastEncounter"] = te
        o = p.get("origin") or {}
        anc_c = s.get("ancestry_country") or o.get("country")
        e = B.origin_ext(s.get("ancestry") or o.get("ethnicity_or_ancestry"), anc_c, s.get("care_country"), s.get("care_basis"),
                         s.get("migrant"), p.get("consanguinity"), s.get("ancestry_evidence") or o.get("evidence"), s.get("care_evidence"))
        if p.get("age"): e.append(ext("ageAsReported", {"text": str(p["age"])}))
        if p.get("diagnosis"): e.append(ext("diagnosisAsReported", {"text": str(p["diagnosis"])}))
        if zyg: e.append(ext("zygosityAsReported", {"text": zyg}))
        if unm: e.append(ext("unmappedPhenotypeText", {"text": json.dumps(unm, ensure_ascii=False)}))
        e.append(ext("pmcPatientsUid", {"text": uid}))
        subj["extensions"] = e + B.analysis_ext(row)
        pp = {"id": pid, "subject": subj, "phenotypicFeatures": feats, "interpretations": B.interp(pid, row["gene"], zyg),
              "metaData": B.meta(row["pmid"], "pmcp")}
        write("pmcp", pp, row, {})

    # ---- Phenopacket Store (copied, augmented)
    orig = {}
    for l in open("results/origin_merged_v3.jsonl"):
        d = json.loads(l)
        if d.get("ok"): orig[d["pmid"]] = d.get("result") or {}
    rows = {r["id"]: r for r in csv.DictReader(open("results/cases_store.tsv"), delimiter="\t")}
    for fpath in glob.glob(os.path.join(a.sp, "ppktstore", "0.1.26", "*", "*.json")):
        try: pp = json.load(open(fpath))
        except Exception: continue
        row = rows.get(pp.get("id"))
        if not row: continue
        pp = copy.deepcopy(pp); r = orig.get(row["pmid"], {})
        e = list(pp["subject"].get("extensions", []))
        po = [o for o in r.get("patient_origins") or [] if o.get("country") or o.get("ethnicity_or_ancestry")]
        if len(po) == 1:
            if po[0].get("country"): e.append(ext("geographicOrigin", B.country(po[0]["country"])))
            if po[0].get("ethnicity_or_ancestry"): e.append(ext("ancestry", {"text": po[0]["ethnicity_or_ancestry"]}))
            e.append(ext("originLevel", {"text": "publication (single origin stated for the paper's patients)"}))
        elif po:
            e.append(ext("publicationPatientOrigins", {"text": json.dumps(po, ensure_ascii=False)}))
            e.append(ext("originLevel", {"text": "publication (several origins; not resolved to this patient)"}))
        rc = [x.get("country") for x in r.get("recruitment_countries") or [] if x.get("country")]
        if rc: e.append(ext("recruitmentCountries", {"text": "; ".join(rc)}))
        if r.get("consanguinity") == "yes": e.append(ext("consanguinityInPublication", {"text": r.get("consanguinity_evidence") or "yes"}))
        pp["subject"]["extensions"] = e + B.analysis_ext(row)
        pp["metaData"].setdefault("externalReferences", []).append({"id": "pavs-population:store", "description": CORPUS_DESC["store"]})
        write("store", pp, row, {})

    # ---- PAVS (copied, augmented)
    rows = {r["id"]: r for r in csv.DictReader(open("results/cases_pavs.tsv"), delimiter="\t")}
    for pp in json.load(open(a.pavs_json)):
        row = rows.get(pp["id"])
        if not row: continue
        pp = copy.deepcopy(pp); e = list(pp["subject"].get("extensions", []))
        cr = pp["metaData"].get("created", "")
        if re.fullmatch(r"\d{4}-\d{2}-\d{2}", cr):  # PAVS writes a bare date; v2 needs a timestamp
            pp["metaData"]["created"] = cr + "T00:00:00Z"
        for dz in pp.get("diseases", []):  # PAVS puts modeOfInheritance on Disease; v2 has no such field there
            if "modeOfInheritance" in dz: e.append(ext("modeOfInheritance", {"ontologyClass": dz.pop("modeOfInheritance")}))
        src = pp["id"][5]
        if src in CARE_SOURCE:
            e.append(ext("placeOfCare", B.country(CARE_SOURCE[src])))
            e.append(ext("placeOfCareBasis", {"text": "source dataset (clinical records of that institution)"}))
        # PAVS records ancestry as a HANCESTRO "phenotypic feature"; move it to the subject
        for f in pp["phenotypicFeatures"]:
            if f["type"]["id"].startswith("HANCESTRO:"): e.append(ext("ancestry", {"ontologyClass": f["type"]}))
        pp["phenotypicFeatures"] = [f for f in pp["phenotypicFeatures"] if not f["type"]["id"].startswith("HANCESTRO:")]
        pp["subject"]["extensions"] = e + B.analysis_ext(row)
        pp["metaData"].setdefault("externalReferences", []).append({"id": f"pavs-population:{row['corpus']}", "description": CORPUS_DESC[row["corpus"]]})
        write(row["corpus"], pp, row, {})

    with open(os.path.join(a.out, "index.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(index[0]), delimiter="\t"); w.writeheader(); w.writerows(index)
    print(dict(st), "total", sum(st.values()))


if __name__ == "__main__":
    main()
