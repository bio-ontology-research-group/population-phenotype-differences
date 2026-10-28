#!/usr/bin/env python3
"""Additional Qwen passes for the provenance/population analysis (BORG endpoint on unimatrix).

Stage 'setting': for each Mendelian PMC-Patients case, separate ancestry from place of care.
  --cands genetic_candidates.json --patients patients.jsonl --affil pmcp_pmid_country.tsv --out setting.jsonl
Stage 'termclass': classify HPO terms by how the finding is obtained (bedside vs a test somebody
  had to decide to order). --terms terms.tsv (id, name, definition) --out termclass.jsonl
Resumable; the API key is read from a file and never printed.
"""
import argparse, csv, json, os, queue, re, threading, time, urllib.request

SYSTEM = "You are a careful clinical geneticist. Use only what the text states. Answer with JSON only."

SETTING_PROMPT = """Below is the description of ONE patient from a published case report, with the article title and the countries of the authors' affiliations.

Separate the patient's ANCESTRY from the PLACE where the patient was examined and treated.

- ancestry: ethnicity, ancestry or country of origin of the patient or family as stated in the text (e.g. "Saudi", "Turkish", "of Moroccan descent", "Han Chinese"). If the text says nothing about origin or ethnicity, give null. Do not infer ancestry from the hospital location alone.
- care_country: country where the patient was examined, diagnosed or treated. Use "explicit" if the text names the country, city or hospital; use "our_institution" if the text only says "our hospital/clinic/department" or similar, in which case take the country of the authors' affiliation that the text refers to (usually the first author's); use "none" if neither applies.
- migrant: "yes" if the text indicates the patient or family lives in, or migrated to, a country different from their stated origin (e.g. "a Turkish family living in Germany", "refugee", "immigrant", "of Pakistani descent born in the UK"); "no" if origin and place of care are the same country or the text implies a local patient; "unknown" otherwise.

Return:
{"ancestry": "<as stated or null>", "ancestry_country": "<country implied by the stated ancestry or null>", "ancestry_evidence": "<quote, max 20 words, or null>",
 "care_country": "<country or null>", "care_basis": "explicit|our_institution|none", "care_evidence": "<quote, max 20 words, or null>",
 "migrant": "yes|no|unknown", "migrant_evidence": "<quote or null>"}

Article title: {title}
Author affiliation countries (first author first): {affil}
PATIENT DESCRIPTION:
{text}
"""

TERMCLASS_PROMPT = """For each Human Phenotype Ontology term below, say how a clinician would normally establish that a patient has it.

Classes:
- symptom_history: reported by the patient/family or taken from the history (pain, fatigue, developmental milestones, seizures witnessed, feeding difficulty).
- clinical_exam: seen or measured at the bedside without special equipment beyond simple tools (dysmorphic features, skin findings, hypotonia, ataxia on examination, hepatomegaly by palpation, height/weight, blood pressure, fever).
- laboratory: requires a laboratory test on blood, urine, CSF or other fluid (chemistry, metabolites, hormones, enzymes, cell counts, immunoglobulins, lymphocyte subsets, antibodies).
- imaging: requires radiography, ultrasound, echocardiography, CT, MRI or other imaging.
- physiological_test: requires an instrument-based functional test (ECG, Holter, EEG, EMG/nerve conduction, audiometry, ERG, visual fields, spirometry, polysomnography, slit-lamp or fundus examination by an ophthalmologist).
- pathology: requires biopsy, histology, microscopy of tissue or cells, or muscle/skin/liver biopsy findings.
- mixed: commonly established in more than one of these ways, or a broad grouping term.

Also say whether establishing the term usually requires a DECISION TO ORDER a test (true for laboratory, imaging, physiological_test, pathology; false for symptom_history and clinical_exam; for mixed decide by the most common route).

Return a JSON list with one object per term, in the same order: [{"id": "HP:...", "class": "...", "needs_test": true|false}]

Terms:
{terms}
"""


SELECT_PROMPT = """You map clinical findings from case reports to the Human Phenotype Ontology (HPO).
For each numbered finding choose exactly ONE term from that finding's candidate list:
- Choose the candidate that means the same as the finding. If none means the same, choose the closest candidate that is MORE GENERAL and still correct (a parent-like term).
- Never choose a term that asserts more than the finding (a more specific type, a particular body part, a cause, or a diagnosis the text does not state).
- Match the anatomy and the kind of abnormality exactly; similar-sounding words in different body systems are wrong (e.g. limb vs limbic, hand vs foot).
- Laboratory values: choose an increased/decreased term only if the finding says the value is abnormal in that direction.
- The proposed label is only a hint from another model and may be wrong.
- Answer null if no candidate is correct, or if the finding is not an abnormality.
Return a JSON list: [{"n": 1, "hpo": "HP:..." or null}, ...] covering every finding.

{items}
"""


ABSTRACT_PROMPT = """Below is the abstract of a published case report, with its title and the countries of the authors' affiliations.

Task 1. List each PATIENT described who has a GERMLINE monogenic (Mendelian) disorder whose causal gene is identified (by sequencing or another genetic test). Somatic tumour mutations, pharmacogenetic findings, risk alleles, aneuploidies and large copy-number syndromes do not count; a single-gene deletion or intragenic CNV does count. If several affected relatives are described only together, list them as one entry. If no patient qualifies, return an empty list.

Task 2. For each patient give:
- causal gene(s) and zygosity;
- observed phenotypes: every distinct clinical abnormality of that patient (signs, symptoms, examination, laboratory, imaging, electrophysiology, histology), each as the verbatim phrase plus the best matching HPO term label; do not list the disease name, genetic results, treatments or normal findings;
- ancestry: ethnicity, ancestry or country of origin as stated (null if not stated; do not infer it from the hospital location);
- care_country: where the patient was examined/treated: "explicit" if the abstract names the country, city or hospital; "our_institution" if it only says "our hospital/clinic" (then use the authors' affiliation country); "none" otherwise;
- migrant: "yes" if the patient/family lives in a country different from their stated origin, "no" if origin and place of care coincide, "unknown" otherwise;
- consanguinity: "reported", "denied" or "not_mentioned".

Return one JSON object:
{"patients": [{"causal_genes": [{"symbol": "<HGNC>", "zygosity": "homozygous|compound_heterozygous|heterozygous|hemizygous|de_novo|unknown"}],
  "sex": "male|female|unknown", "age": "<age at report/examination as stated, or null>",
  "age_at_onset": "<age at onset of symptoms as stated (e.g. 'birth', '6 months', 'early childhood'), or null>",
  "observed_phenotypes": [{"text": "<verbatim>", "hpo_label": "<HPO label>"}],
  "ancestry": "<as stated or null>", "ancestry_country": "<country or null>",
  "care_country": "<country or null>", "care_basis": "explicit|our_institution|none",
  "migrant": "yes|no|unknown", "consanguinity": "reported|denied|not_mentioned"}]}
Output JSON only.

Title: {title}
Author affiliation countries (first author first): {affil}
ABSTRACT:
{text}
"""


def post(a, key, prompt, max_tokens):
    body = {"model": a.model, "temperature": 0.0, "max_tokens": max_tokens, "chat_template_kwargs": {"enable_thinking": False},
            "messages": [{"role": "system", "content": SYSTEM}, {"role": "user", "content": prompt}]}
    req = urllib.request.Request(f"{a.endpoint}/v1/chat/completions", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json", "Authorization": f"Bearer {key}"})
    err = None
    for i in range(4):
        try:
            with urllib.request.urlopen(req, timeout=600) as r:
                return json.loads(r.read())["choices"][0]["message"]["content"], None
        except Exception as e:  # noqa: BLE001
            err = str(e); time.sleep(5 * (i + 1))
    return None, err


def run_pool(items, fn, out_path, jobs, label):
    lock = threading.Lock(); out = open(out_path, "a"); q = queue.Queue(); n = {"n": 0}
    for it in items: q.put(it)
    def worker():
        while True:
            try: it = q.get_nowait()
            except queue.Empty: return
            recs = fn(it)
            with lock:
                for r in recs: out.write(json.dumps(r) + "\n")
                out.flush(); n["n"] += 1
                if n["n"] % 100 == 0: print(f"{label}: {n['n']}/{len(items)}", flush=True)
    ts = [threading.Thread(target=worker, daemon=True) for _ in range(jobs)]
    [t.start() for t in ts]; [t.join() for t in ts]
    print(f"{label}: finished", flush=True)


def done_keys(path, field):
    s = set()
    if os.path.exists(path):
        for l in open(path):
            try:
                d = json.loads(l)
                if d.get("ok"): s.add(d[field])
            except Exception:  # noqa: BLE001
                pass
    return s


def stage_setting(a):
    key = open(os.path.expanduser(a.key_file)).read().strip()
    mend = set()
    for l in open(a.patients):
        d = json.loads(l)
        if d.get("ok") and (d.get("result") or {}).get("mendelian_germline"): mend.add(d["uid"])
    affil = {}
    for r in csv.DictReader(open(a.affil), delimiter="\t"):
        cs = [r["first_country"]] + [x.split(":")[0] for x in (r["all_countries"] or "").split(";") if x]
        affil[r["pmid"]] = ", ".join(dict.fromkeys(c for c in cs if c)) or "unknown"
    cands = [c for c in json.load(open(a.cands)) if c["uid"] in mend]
    done = done_keys(a.out, "uid"); todo = [c for c in cands if c["uid"] not in done]
    print(f"setting: {len(cands)} Mendelian, {len(todo)} to do", flush=True)
    def fn(c):
        p = SETTING_PROMPT.replace("{title}", c["title"]).replace("{affil}", affil.get(c["pmid"], "unknown")).replace("{text}", c["text"][:30000])
        txt, err = post(a, key, p, 600)
        parsed = None
        if txt:
            m = re.search(r"\{.*\}", txt, flags=re.S)
            try: parsed = json.loads(m.group(0)) if m else None
            except Exception: parsed = None  # noqa: BLE001
        return [{"uid": c["uid"], "pmid": c["pmid"], "ok": parsed is not None, "result": parsed, "error": err}]
    run_pool(todo, fn, a.out, a.jobs, "setting")


def stage_termclass(a):
    key = open(os.path.expanduser(a.key_file)).read().strip()
    terms = [r for r in csv.DictReader(open(a.terms), delimiter="\t")]
    done = done_keys(a.out, "id"); todo = [t for t in terms if t["id"] not in done]
    batches = [todo[i:i + 25] for i in range(0, len(todo), 25)]
    print(f"termclass: {len(terms)} terms, {len(todo)} to do in {len(batches)} batches", flush=True)
    def fn(b):
        lst = "\n".join(f"{t['id']} | {t['name']} | {t['definition'][:300]}" for t in b)
        txt, err = post(a, key, TERMCLASS_PROMPT.replace("{terms}", lst), 3000)
        got = {}
        if txt:
            m = re.search(r"\[.*\]", txt, flags=re.S)
            try:
                for x in json.loads(m.group(0)) if m else []:
                    if isinstance(x, dict) and x.get("id"): got[x["id"]] = x
            except Exception:  # noqa: BLE001
                pass
        return [{"id": t["id"], "ok": t["id"] in got, "class": got.get(t["id"], {}).get("class"),
                 "needs_test": got.get(t["id"], {}).get("needs_test")} for t in b]
    run_pool(batches, fn, a.out, a.jobs, "termclass")


def stage_select(a):
    key = open(os.path.expanduser(a.key_file)).read().strip()
    rows = [json.loads(l) for l in open(a.cands)]
    done = done_keys(a.out, "key")
    auto = [r for r in rows if r.get("text_exact") and r["key"] not in done]
    with open(a.out, "a") as f:
        for r in auto: f.write(json.dumps({"key": r["key"], "ok": True, "hpo": r["text_exact"], "route": "text_exact"}) + "\n")
    todo = [r for r in rows if not r.get("text_exact") and r["key"] not in done]
    batches = [todo[i:i + a.batch] for i in range(0, len(todo), a.batch)]
    print(f"select: {len(rows)} phrases, {len(auto)} exact text matches, {len(todo)} to LLM in {len(batches)} batches", flush=True)
    def fn(b):
        items = []
        for i, r in enumerate(b):
            c = "\n".join(f"     {x['id']} {x['name']} (parent: {x['parent']}){' - ' + x['def'] if x['def'] else ''}" for x in r["cands"])
            items.append(f"{i+1}. finding: \"{r['text']}\"  proposed label: \"{r['label']}\"\n   candidates:\n{c}")
        txt, err = post(a, key, SELECT_PROMPT.replace("{items}", "\n\n".join(items)), 200 + 40 * len(b))
        got = {}
        if txt:
            m = re.search(r"\[.*\]", txt, flags=re.S)
            try:
                for x in json.loads(m.group(0)) if m else []:
                    if isinstance(x, dict) and "n" in x: got[int(x["n"])] = x.get("hpo")
            except Exception:  # noqa: BLE001
                pass
        out = []
        for i, r in enumerate(b):
            h = got.get(i + 1, "MISSING")
            ok = h != "MISSING"
            if ok and h is not None and h not in {x["id"] for x in r["cands"]}: h = None   # only listed candidates
            out.append({"key": r["key"], "ok": ok, "hpo": h if ok else None, "route": "llm"})
        return out
    run_pool(batches, fn, a.out, a.jobs, "select")


def stage_abstract(a):
    key = open(os.path.expanduser(a.key_file)).read().strip()
    cands = [json.loads(l) for l in open(a.cands)]
    done = done_keys(a.out, "pmid"); todo = [c for c in cands if c["pmid"] not in done]
    print(f"abstract: {len(cands)} candidates, {len(todo)} to do", flush=True)
    def fn(c):
        p = ABSTRACT_PROMPT.replace("{title}", c["title"]).replace("{affil}", c.get("affil") or "unknown").replace("{text}", c["abstract"][:12000])
        txt, err = post(a, key, p, 2500)
        parsed = None
        if txt:
            m = re.search(r"\{.*\}", txt, flags=re.S)
            try: parsed = json.loads(m.group(0)) if m else None
            except Exception: parsed = None  # noqa: BLE001
        return [{"pmid": c["pmid"], "ok": parsed is not None, "result": parsed, "error": err}]
    run_pool(todo, fn, a.out, a.jobs, "abstract")


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="stage", required=True)
    for name in ("setting", "termclass", "select", "abstract"):
        p = sub.add_parser(name)
        p.add_argument("--out", required=True); p.add_argument("--key-file", required=True)
        p.add_argument("--endpoint", default="http://unimatrix01.kaust.edu.sa:8000"); p.add_argument("--model", default="qwen3.8-27b")
        p.add_argument("--jobs", type=int, default=12)
        if name == "setting":
            p.add_argument("--cands", required=True); p.add_argument("--patients", required=True); p.add_argument("--affil", required=True)
        elif name == "abstract":
            p.add_argument("--cands", required=True)
        elif name == "termclass":
            p.add_argument("--terms", required=True)
        else:
            p.add_argument("--cands", required=True); p.add_argument("--batch", type=int, default=8)
    a = ap.parse_args(); {"setting": stage_setting, "termclass": stage_termclass, "select": stage_select, "abstract": stage_abstract}[a.stage](a)


if __name__ == "__main__":
    main()
