#!/usr/bin/env python3
"""Per-patient extraction from PMC-Patients (Zhao et al. 2023) with the BORG Qwen endpoint.

Stage 'patient': for each candidate patient summary (genetic_candidates.json, made by filtering
PMC-Patients-V2 for HPOA gene symbols plus variant vocabulary), decide whether it describes a
germline Mendelian diagnosis, and extract causal gene(s), observed and excluded phenotypes as
verbatim phrases with a proposed HPO label, sex, age, and any origin/consanguinity stated in the
summary. Article-level origin is extracted separately with store_origin_extract.py (fetch + llm).

Usage:
  python3 pmcp_extract.py patient --cands genetic_candidates.json --out patients.jsonl \
      --key-file ~/.borg_llm_key [--limit N] [--jobs 12]
Resumable; the API key is read from a file and never printed.
"""
import argparse, collections, json, os, re, threading, queue, time, urllib.request

SYSTEM = ("You are a careful clinical geneticist extracting structured data from published case "
          "descriptions. Use only what the text states. Answer with JSON only.")

PROMPT = """Below is the description of ONE patient from a published case report.

Task 1. Decide whether this patient has a GERMLINE monogenic (Mendelian) disorder whose causal gene is identified in the text (by sequencing or another genetic test). Somatic tumour mutations, pharmacogenetic findings, risk alleles, chromosomal aneuploidies and large copy-number syndromes do not count; a single-gene deletion or intragenic CNV does count.

Task 2. If yes, extract the patient's phenotype: every distinct clinical abnormality of THIS patient (signs, symptoms, examination, laboratory, imaging, electrophysiology, histology findings). Do not list the disease name itself, genetic results, treatments, or normal findings. Findings that the text explicitly says are ABSENT go into excluded_phenotypes. For each finding give the verbatim phrase and the best matching Human Phenotype Ontology term label (the exact HPO label if you know it; otherwise a short standard clinical term).

Return one JSON object:
{
  "mendelian_germline": true or false,
  "reason": "<one short sentence>",
  "causal_genes": [{"symbol": "<HGNC symbol>", "zygosity": "homozygous|compound_heterozygous|heterozygous|hemizygous|de_novo|unknown", "evidence": "<quote, max 20 words>"}],
  "diagnosis": "<disease name as stated, or null>",
  "sex": "male|female|unknown",
  "age": "<age at report as stated, or null>",
  "observed_phenotypes": [{"text": "<verbatim phrase>", "hpo_label": "<HPO label>"}],
  "excluded_phenotypes": [{"text": "<verbatim phrase>", "hpo_label": "<HPO label>"}],
  "origin": {"country": "<country or null>", "ethnicity_or_ancestry": "<as stated or null>", "evidence": "<quote or null>"},
  "consanguinity": "reported|denied|not_mentioned"
}
If mendelian_germline is false, return empty lists for the phenotype fields. Output JSON only, no markdown.

PATIENT DESCRIPTION (article PMID {pmid}, title: {title}):
{text}
"""


def call(a, key, item):
    body = {"model": a.model, "temperature": 0.0, "max_tokens": 3000,
            "chat_template_kwargs": {"enable_thinking": False},
            "messages": [{"role": "system", "content": SYSTEM},
                         {"role": "user", "content": PROMPT.replace("{pmid}", item["pmid"]).replace("{title}", item["title"]).replace("{text}", item["text"][:30000])}]}
    req = urllib.request.Request(f"{a.endpoint}/v1/chat/completions", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json", "Authorization": f"Bearer {key}"})
    err = None
    for i in range(4):
        try:
            with urllib.request.urlopen(req, timeout=600) as r:
                d = json.loads(r.read())
            c = d["choices"][0]["message"]["content"]
            m = re.search(r"\{.*\}", c, flags=re.S)
            try:
                parsed = json.loads(m.group(0)) if m else None
            except Exception:  # noqa: BLE001
                parsed = None
            u = d.get("usage", {})
            return {"uid": item["uid"], "pmid": item["pmid"], "ok": parsed is not None, "result": parsed,
                    "raw": None if parsed else c[:3000], "prompt_tokens": u.get("prompt_tokens"), "completion_tokens": u.get("completion_tokens")}
        except Exception as e:  # noqa: BLE001
            err = str(e); time.sleep(5 * (i + 1))
    return {"uid": item["uid"], "pmid": item["pmid"], "ok": False, "result": None, "error": err}


def stage_patient(a):
    key = open(os.path.expanduser(a.key_file)).read().strip()
    cands = json.load(open(a.cands))
    if a.limit:
        import random
        random.Random(a.seed).shuffle(cands); cands = cands[:a.limit]
    done = set()
    if os.path.exists(a.out):
        for l in open(a.out):
            try:
                d = json.loads(l)
                if d.get("ok"): done.add(d["uid"])
            except Exception:  # noqa: BLE001
                pass
    todo = [c for c in cands if c["uid"] not in done]
    print(f"patient: {len(done)} done, {len(todo)} to do", flush=True)
    q = queue.Queue(); [q.put(c) for c in todo]
    lock = threading.Lock(); out = open(a.out, "a"); n = {"n": 0}
    def worker():
        while True:
            try: c = q.get_nowait()
            except queue.Empty: return
            r = call(a, key, c)
            with lock:
                out.write(json.dumps(r) + "\n"); out.flush(); n["n"] += 1
                if n["n"] % 100 == 0: print(f"patient: {n['n']}/{len(todo)}", flush=True)
    ts = [threading.Thread(target=worker, daemon=True) for _ in range(a.jobs)]
    [t.start() for t in ts]; [t.join() for t in ts]
    print("patient: finished", flush=True)


# ------------------------------------------------------------------ HPO mapping of unmatched labels
def _norm(t):
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9 ]", " ", (t or "").lower())).strip()


def _toks(t):
    return {w[:-1] if len(w) > 4 and w.endswith("s") else w for w in _norm(t).split() if len(w) > 1}


def load_obo_full(path):
    terms, cur = [], None
    for line in open(path):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = {"syn": [], "allsyn": []}; terms.append(cur)
        elif cur is not None:
            if line.startswith("id: "): cur["id"] = line[4:]
            elif line.startswith("name: "): cur["name"] = line[6:]
            elif line.startswith("synonym: "):
                cur["allsyn"].append(line.split('"')[1])
                if " EXACT" in line: cur["syn"].append(line.split('"')[1])
            elif line.startswith("is_obsolete: true"): cur["obs"] = True
            elif line == "": cur = None
    return [t for t in terms if t.get("id", "").startswith("HP:") and not t.get("obs") and "name" in t]


MAP_PROMPT = """For each numbered clinical finding, choose the Human Phenotype Ontology term that means the same thing or is the closest more general term, from that finding's candidate list. Never choose a term that is more specific than the finding (for example, "seizures" must not become a particular seizure type, and "limited range of motion" must not become a particular joint); prefer the closest more general candidate instead. Answer "none" if no candidate fits (for example if the finding is a disease name, a test result without abnormality, or a treatment).
Return JSON only: {"answers": [{"n": <number>, "id": "HP:xxxxxxx" or "none"}]}

"""


def stage_map(a):
    import math
    key = open(os.path.expanduser(a.key_file)).read().strip()
    terms = load_obo_full(a.obo)
    exact = {}
    for t in terms:
        for x in [t["name"]] + t["syn"]: exact.setdefault(_norm(x), t["id"])
    name = {t["id"]: t["name"] for t in terms}
    # inverted index over names + all synonyms
    docs = []
    for t in terms:
        for x in [t["name"]] + t["allsyn"]: docs.append((t["id"], _toks(x)))
    df = collections.Counter(w for _, ts in docs for w in ts)
    inv = collections.defaultdict(list)
    for k, (_, ts) in enumerate(docs):
        for w in ts: inv[w].append(k)
    N = len(docs); idf = {w: math.log(N / c) for w, c in df.items()}
    def cands(q, n=15):
        qt = _toks(q); sc = collections.Counter()
        for w in qt:
            for k in inv.get(w, ()): sc[k] += idf[w]
        best = {}
        for k, v in sc.most_common(400):
            hid = docs[k][0]; v = v / (1 + 0.15 * len(docs[k][1] - qt))
            if v > best.get(hid, 0): best[hid] = v
        return [h for h, _ in sorted(best.items(), key=lambda x: -x[1])[:n]]
    # unique unmatched (label, text) items
    items = {}
    for l in open(a.patients):
        d = json.loads(l)
        if not d.get("ok") or not (d["result"] or {}).get("mendelian_germline"): continue
        for p in (d["result"].get("observed_phenotypes") or []) + (d["result"].get("excluded_phenotypes") or []):
            if isinstance(p, str): p = {"text": p, "hpo_label": p}
            if not isinstance(p, dict): continue
            lab = p.get("hpo_label") or ""
            if _norm(lab) in exact or _norm(p.get("text")) in exact or not lab: continue
            items.setdefault(_norm(lab), (lab, p.get("text") or ""))
    done = {}
    if os.path.exists(a.out):
        for l in open(a.out):
            r = json.loads(l); done[r["key"]] = r
    todo = [(k, v) for k, v in items.items() if k not in done]
    print(f"map: {len(items)} unique unmatched labels, {len(todo)} to do", flush=True)
    batches = [todo[i:i + 20] for i in range(0, len(todo), 20)]
    q = queue.Queue(); [q.put(b) for b in batches]
    lock = threading.Lock(); out = open(a.out, "a"); cnt = {"n": 0}
    def worker():
        while True:
            try: b = q.get_nowait()
            except queue.Empty: return
            cl = [cands(lab + " " + txt) for _, (lab, txt) in b]
            txt = MAP_PROMPT + "\n".join(f"{i+1}. finding: \"{lab}\" (text: \"{t[:120]}\")\n   candidates: " + "; ".join(f"{h} {name[h]}" for h in c)
                                       for i, ((_, (lab, t)), c) in enumerate(zip(b, cl)))
            body = {"model": a.model, "temperature": 0.0, "max_tokens": 1500, "chat_template_kwargs": {"enable_thinking": False},
                    "messages": [{"role": "user", "content": txt}]}
            req = urllib.request.Request(f"{a.endpoint}/v1/chat/completions", data=json.dumps(body).encode(),
                                         headers={"Content-Type": "application/json", "Authorization": f"Bearer {key}"})
            ans = {}
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(req, timeout=600) as r:
                        c = json.loads(r.read())["choices"][0]["message"]["content"]
                    m = re.search(r"\{.*\}", c, flags=re.S)
                    for x in json.loads(m.group(0))["answers"]: ans[int(x["n"])] = x.get("id")
                    break
                except Exception:  # noqa: BLE001
                    time.sleep(5 * (attempt + 1))
            with lock:
                for i, ((k, (lab, t)), c) in enumerate(zip(b, cl)):
                    h = ans.get(i + 1)
                    h = h if (h in c) else None   # only accept a listed candidate
                    out.write(json.dumps({"key": k, "label": lab, "hpo": h, "n_cands": len(c)}) + "\n")
                out.flush(); cnt["n"] += 1
                if cnt["n"] % 50 == 0: print(f"map: {cnt['n']}/{len(batches)} batches", flush=True)
    ts = [threading.Thread(target=worker, daemon=True) for _ in range(a.jobs)]
    [t.start() for t in ts]; [t.join() for t in ts]
    print("map: finished", flush=True)


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="stage", required=True)
    p = sub.add_parser("patient")
    p.add_argument("--cands", required=True); p.add_argument("--out", required=True); p.add_argument("--key-file", required=True)
    p.add_argument("--endpoint", default=os.environ.get("LLM_ENDPOINT", "http://localhost:8000")); p.add_argument("--model", default="qwen3.8-27b")
    p.add_argument("--jobs", type=int, default=12); p.add_argument("--limit", type=int, default=0); p.add_argument("--seed", type=int, default=1)
    m = sub.add_parser("map")
    m.add_argument("--patients", required=True); m.add_argument("--obo", required=True); m.add_argument("--out", required=True)
    m.add_argument("--key-file", required=True); m.add_argument("--endpoint", default=os.environ.get("LLM_ENDPOINT", "http://localhost:8000"))
    m.add_argument("--model", default="qwen3.8-27b"); m.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args(); {"patient": stage_patient, "map": stage_map}[a.stage](a)


if __name__ == "__main__":
    main()
