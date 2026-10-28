#!/usr/bin/env python3
"""Extract patient geographic/ethnic origin from the source articles of the
GA4GH Phenopacket Store using Europe PMC full text and the BORG Qwen endpoint.

Stages (both resumable, output is JSONL keyed by PMID):
  fetch : Europe PMC fullTextXML (open-access subset) or title+abstract fallback
          -> <workdir>/text/<pmid>.txt and <workdir>/fetch.jsonl
  llm   : one chat completion per article -> <workdir>/origin.jsonl

Usage:
  python3 store_origin_extract.py fetch --pmids store_pmids_union.txt --workdir W
  python3 store_origin_extract.py llm   --workdir W --key-file ~/.borg_llm_key \
          --endpoint http://localhost:8000 --model qwen3.8-27b --jobs 12
The API key is read from a file and never printed.
"""
import argparse, json, os, re, sys, time, threading, queue, urllib.request, urllib.parse
import xml.etree.ElementTree as ET

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MAX_CHARS = 90000  # ~22k tokens; Qwen context is 131k


def http_get(url, timeout=120, retries=4, headers=None):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers or {"User-Agent": "pavs-origin/0.1"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"GET failed {url}: {last}")


def xml_text(elem):
    return " ".join(t.strip() for t in elem.itertext() if t and t.strip())


def fetch_one(pmid):
    """Return dict(pmid, source, text, pmcid)."""
    # 1. search record (gives PMCID + isOpenAccess)
    q = urllib.parse.quote(f"EXT_ID:{pmid} AND SRC:MED")
    rec = json.loads(http_get(f"{EPMC}/search?query={q}&format=json&resultType=core"))
    hits = rec.get("resultList", {}).get("result", [])
    pmcid, oa, title, abstract = None, False, "", ""
    if hits:
        h = hits[0]
        pmcid = h.get("pmcid")
        oa = h.get("isOpenAccess") == "Y"
        title = h.get("title", "")
        abstract = h.get("abstractText", "")
    if pmcid:
        try:
            xml = http_get(f"{EPMC}/{pmcid}/fullTextXML")
            root = ET.fromstring(xml)
            body = root.find(".//body")
            parts = [title]
            if abstract:
                parts.append("ABSTRACT: " + re.sub(r"<[^>]+>", "", abstract))
            if body is not None:
                # keep section titles to help the model find Patients/Methods
                for sec in body.iter("sec"):
                    t = sec.find("title")
                    if t is not None:
                        parts.append("\n## " + xml_text(t))
                    for p in sec.findall("p"):
                        parts.append(xml_text(p))
                # tables often hold "Origin"/"Ethnicity" columns
                for tw in body.iter("table-wrap"):
                    parts.append("\nTABLE: " + xml_text(tw)[:6000])
            txt = "\n".join(parts)
            if len(txt) > 2000:
                return {"pmid": pmid, "pmcid": pmcid, "source": "fulltext", "text": txt[:MAX_CHARS]}
        except Exception as e:  # noqa: BLE001
            sys.stderr.write(f"{pmid}: fulltext failed ({e}); falling back\n")
    if not (title or abstract):
        # PubMed fallback for title/abstract
        xml = http_get(f"{EUTILS}/efetch.fcgi?db=pubmed&retmode=xml&id={pmid}")
        root = ET.fromstring(xml)
        title = root.findtext(".//ArticleTitle") or ""
        abstract = " ".join(t.text or "" for t in root.findall(".//AbstractText"))
    txt = f"{title}\nABSTRACT: {re.sub(r'<[^>]+>', '', abstract)}"
    return {"pmid": pmid, "pmcid": pmcid, "source": "abstract", "text": txt[:MAX_CHARS]}


def stage_fetch(a):
    os.makedirs(os.path.join(a.workdir, "text"), exist_ok=True)
    done = set()
    out_path = os.path.join(a.workdir, "fetch.jsonl")
    if os.path.exists(out_path):
        for line in open(out_path):
            try:
                done.add(json.loads(line)["pmid"])
            except Exception:  # noqa: BLE001
                pass
    pmids = [p for p in open(a.pmids).read().split() if p not in done]
    print(f"fetch: {len(done)} done, {len(pmids)} to do", flush=True)
    with open(out_path, "a") as out:
        for i, pmid in enumerate(pmids):
            try:
                r = fetch_one(pmid)
            except Exception as e:  # noqa: BLE001
                r = {"pmid": pmid, "pmcid": None, "source": "error", "text": "", "error": str(e)}
            with open(os.path.join(a.workdir, "text", f"{pmid}.txt"), "w") as f:
                f.write(r["text"])
            meta = {k: v for k, v in r.items() if k != "text"}
            meta["chars"] = len(r["text"])
            out.write(json.dumps(meta) + "\n")
            out.flush()
            if i % 50 == 0:
                print(f"fetch: {i}/{len(pmids)} {pmid} {meta['source']} {meta['chars']}", flush=True)
            time.sleep(0.34)


def _clean_html(b):
    t = b.decode("utf-8", "ignore")
    t = re.sub(r"(?is)<(script|style|nav|header|footer)[^>]*>.*?</\1>", " ", t)
    t = re.sub(r"(?s)<[^>]+>", " ", t)
    import html as _h
    return re.sub(r"\s+", " ", _h.unescape(t))


def fetch_oa_one(pmid, info, tmpdir):
    """Try OpenAlex open-access locations: PDFs first, then the OA landing page.
    PMC/Europe PMC hosts are skipped (non-OA-subset bulk retrieval is not permitted;
    the OA subset was already read in the 'fetch' stage)."""
    import subprocess
    ua = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) pavs-origin research fetch"}
    skip = ("ncbi.nlm.nih.gov", "europepmc.org")
    urls = [u for u in info.get("pdfs", []) if u and not any(h in u for h in skip)]
    if info.get("oa_url") and not any(h in info["oa_url"] for h in skip):
        urls.append(info["oa_url"])
    for u in dict.fromkeys(urls):
        try:
            b = http_get(u, timeout=90, retries=2, headers=ua)
        except Exception:  # noqa: BLE001
            continue
        if b[:5] == b"%PDF-":
            pdf = os.path.join(tmpdir, f"{pmid}.pdf")
            open(pdf, "wb").write(b)
            r = subprocess.run(["pdftotext", "-layout", pdf, "-"], capture_output=True, timeout=120)
            os.remove(pdf)
            txt = r.stdout.decode("utf-8", "ignore")
            kind = "oa_pdf"
        else:
            txt = _clean_html(b)
            kind = "oa_html"
        if len(txt) > 6000:
            return {"pmid": pmid, "pmcid": None, "source": kind, "url": u, "text": txt[:MAX_CHARS]}
    return None


def stage_fetch_oa(a):
    os.makedirs(os.path.join(a.workdir, "text"), exist_ok=True)
    tmp = os.path.join(a.workdir, "tmp"); os.makedirs(tmp, exist_ok=True)
    info = json.load(open(a.openalex))
    out_path = os.path.join(a.workdir, "fetch.jsonl")
    done = set()
    if os.path.exists(out_path):
        done = {json.loads(l)["pmid"] for l in open(out_path)}
    todo = [p for p in info if p not in done]
    print(f"fetch-oa: {len(todo)} to do", flush=True)
    n_ok = 0
    with open(out_path, "a") as out:
        for i, pmid in enumerate(todo):
            try:
                r = fetch_oa_one(pmid, info[pmid], tmp)
            except Exception as e:  # noqa: BLE001
                r = None
            if r is None:
                out.write(json.dumps({"pmid": pmid, "source": "none", "chars": 0}) + "\n"); out.flush()
                continue
            n_ok += 1
            open(os.path.join(a.workdir, "text", f"{pmid}.txt"), "w").write(r["text"])
            meta = {k: v for k, v in r.items() if k != "text"}; meta["chars"] = len(r["text"])
            out.write(json.dumps(meta) + "\n"); out.flush()
            if i % 25 == 0:
                print(f"fetch-oa: {i}/{len(todo)} ok={n_ok}", flush=True)
            time.sleep(0.5)
    print(f"fetch-oa: finished ok={n_ok}", flush=True)


SYSTEM = (
    "You are a careful biomedical data extractor. You read a case report or cohort "
    "paper and report where the described patients come from. Answer with JSON only."
)

PROMPT = """Read the article text below and extract the geographic and ethnic origin of the PATIENTS described (not the authors). Use only what the text states; do not infer origin from author affiliations unless the text says the patients were recruited or seen at those centres, in which case report that as 'recruitment_country' with evidence.

Return a single JSON object with these fields:
{
  "n_patients_reported": integer or null,
  "patient_origins": [
    {"country": "<country name or null>", "ethnicity_or_ancestry": "<as stated, or null>", "n_patients": integer or null, "evidence": "<short verbatim quote>"}
  ],
  "recruitment_countries": [ {"country": "<name>", "evidence": "<short quote>"} ],
  "consanguinity": "reported" | "denied" | "not_mentioned",
  "consanguinity_evidence": "<short quote or null>",
  "origin_explicit": true if the text states patient origin/ethnicity explicitly, false if only recruitment site or nothing,
  "confidence": "high" | "medium" | "low",
  "notes": "<one sentence>"
}
Rules: list every distinct origin mentioned for the patients; use null when unknown; keep quotes under 25 words; output JSON only, no prose, no markdown fences.

ARTICLE (PMID {pmid}):
{text}
"""


def llm_call(a, key, pmid, text):
    body = {
        "model": a.model,
        "messages": [
            {"role": "system", "content": SYSTEM},
            {"role": "user", "content": PROMPT.replace("{pmid}", pmid).replace("{text}", text)},
        ],
        "temperature": 0.0,
        "max_tokens": 1500,
        "chat_template_kwargs": {"enable_thinking": False},
    }
    req = urllib.request.Request(
        f"{a.endpoint}/v1/chat/completions",
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", "Authorization": f"Bearer {key}"},
    )
    for i in range(4):
        try:
            with urllib.request.urlopen(req, timeout=600) as r:
                d = json.loads(r.read())
            content = d["choices"][0]["message"]["content"]
            usage = d.get("usage", {})
            m = re.search(r"\{.*\}", content, flags=re.S)
            parsed = json.loads(m.group(0)) if m else None
            return {"pmid": pmid, "ok": parsed is not None, "result": parsed, "raw": None if parsed else content[:2000],
                    "prompt_tokens": usage.get("prompt_tokens"), "completion_tokens": usage.get("completion_tokens")}
        except Exception as e:  # noqa: BLE001
            err = str(e)
            time.sleep(5 * (i + 1))
    return {"pmid": pmid, "ok": False, "result": None, "raw": None, "error": err}


def stage_llm(a):
    key = open(os.path.expanduser(a.key_file)).read().strip()
    out_path = os.path.join(a.workdir, "origin.jsonl")
    done = set()
    if os.path.exists(out_path):
        for line in open(out_path):
            try:
                d = json.loads(line)
                if d.get("ok"):
                    done.add(d["pmid"])
            except Exception:  # noqa: BLE001
                pass
    pmids = []
    for line in open(os.path.join(a.workdir, "fetch.jsonl")):
        d = json.loads(line)
        if d["pmid"] not in done and d.get("chars", 0) > 50:
            pmids.append(d["pmid"])
    pmids = sorted(set(pmids))
    print(f"llm: {len(done)} done, {len(pmids)} to do, jobs={a.jobs}", flush=True)
    q = queue.Queue()
    for p in pmids:
        q.put(p)
    lock = threading.Lock()
    out = open(out_path, "a")
    counter = {"n": 0}

    def worker():
        while True:
            try:
                pmid = q.get_nowait()
            except queue.Empty:
                return
            text = open(os.path.join(a.workdir, "text", f"{pmid}.txt")).read()
            r = llm_call(a, key, pmid, text)
            with lock:
                out.write(json.dumps(r) + "\n")
                out.flush()
                counter["n"] += 1
                if counter["n"] % 25 == 0:
                    print(f"llm: {counter['n']}/{len(pmids)}", flush=True)

    threads = [threading.Thread(target=worker, daemon=True) for _ in range(a.jobs)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    out.close()
    print("llm: finished", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)
    f = sub.add_parser("fetch")
    f.add_argument("--pmids", required=True)
    f.add_argument("--workdir", required=True)
    fo = sub.add_parser("fetch-oa")
    fo.add_argument("--openalex", required=True)
    fo.add_argument("--workdir", required=True)
    l = sub.add_parser("llm")
    l.add_argument("--workdir", required=True)
    l.add_argument("--key-file", required=True)
    l.add_argument("--endpoint", default=os.environ.get("LLM_ENDPOINT", "http://localhost:8000"))
    l.add_argument("--model", default="qwen3.8-27b")
    l.add_argument("--jobs", type=int, default=12)
    a = ap.parse_args()
    {"fetch": stage_fetch, "fetch-oa": stage_fetch_oa, "llm": stage_llm}[a.stage](a)


if __name__ == "__main__":
    main()
