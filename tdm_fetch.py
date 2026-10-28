#!/usr/bin/env python3
"""Fetch full text of subscribed articles through publisher text-and-data-mining APIs
(KAUST subscription, campus IP) for the patient-origin extraction.

  Elsevier: Article Retrieval API (api.elsevier.com, X-ELS-APIKey), plain-text full text.
  Wiley:    TDM API (api.wiley.com/onlinelibrary/tdm/v1, Wiley-TDM-Client-Token), PDF -> pdftotext.

Keys are read from an env file (default ~/.config/pavs-tdm/tdm.env: ELSEVIER_API_KEY=..., WILEY_TDM_TOKEN=...)
and never printed. Output matches store_origin_extract.py: <workdir>/text/<pmid>.txt and
<workdir>/fetch.jsonl ({pmid, doi, source, chars, status}). Resumable; rate-limited.

  python3 tdm_fetch.py --dois pmid_doi.tsv --workdir W [--limit N] [--delay 1.5]
"""
import argparse, json, os, re, subprocess, tempfile, time, urllib.error, urllib.parse, urllib.request

ELSEVIER = ("10.1016", "10.1053", "10.1067", "10.1054", "10.1006", "10.1078")
WILEY = ("10.1002", "10.1111", "10.1046", "10.1113", "10.1034")


def keys(path):
    k = {}
    for l in open(os.path.expanduser(path)):
        if "=" in l: n, v = l.strip().split("=", 1); k[n] = v
    return k


def get(url, headers, timeout=120):
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.status, r.headers.get("Content-Type", ""), r.read()


def elsevier(doi, key):
    st, ct, body = get(f"https://api.elsevier.com/content/article/doi/{urllib.parse.quote(doi)}?httpAccept=text/plain",
                       {"X-ELS-APIKey": key, "Accept": "text/plain"})
    return body.decode("utf-8", "replace"), "elsevier"


def wiley(doi, token):
    st, ct, body = get(f"https://api.wiley.com/onlinelibrary/tdm/v1/articles/{urllib.parse.quote(doi, safe='')}",
                       {"Wiley-TDM-Client-Token": token})
    if not body.startswith(b"%PDF"): return "", "wiley:not_pdf"
    with tempfile.NamedTemporaryFile(suffix=".pdf") as f:
        f.write(body); f.flush()
        txt = subprocess.run(["pdftotext", "-layout", f.name, "-"], capture_output=True, timeout=120).stdout.decode("utf-8", "replace")
    return txt, "wiley"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dois", required=True, help="TSV pmid<TAB>doi")
    ap.add_argument("--workdir", required=True); ap.add_argument("--keys", default="~/.config/pavs-tdm/tdm.env")
    ap.add_argument("--limit", type=int, default=0); ap.add_argument("--delay", type=float, default=1.5)
    ap.add_argument("--min-chars", type=int, default=6000)
    a = ap.parse_args()
    K = keys(a.keys)
    os.makedirs(os.path.join(a.workdir, "text"), exist_ok=True)
    out_path = os.path.join(a.workdir, "fetch.jsonl")
    done = {json.loads(l)["pmid"] for l in open(out_path)} if os.path.exists(out_path) else set()
    rows = [l.rstrip("\n").split("\t") for l in open(a.dois) if "\t" in l]
    todo = [(p, d) for p, d in rows if d and p not in done and d.startswith(ELSEVIER + WILEY)]
    if a.limit: todo = todo[:a.limit]
    print(f"tdm: {len(todo)} to fetch ({sum(d.startswith(ELSEVIER) for _, d in todo)} Elsevier, {sum(d.startswith(WILEY) for _, d in todo)} Wiley)", flush=True)
    n_ok = 0
    with open(out_path, "a") as out:
        for i, (pmid, doi) in enumerate(todo):
            rec = {"pmid": pmid, "doi": doi, "source": "none", "chars": 0}
            try:
                txt, src = elsevier(doi, K["ELSEVIER_API_KEY"]) if doi.startswith(ELSEVIER) else wiley(doi, K["WILEY_TDM_TOKEN"])
                txt = re.sub(r"[ \t]+", " ", txt)
                rec["status"] = src
                if len(txt) >= a.min_chars:
                    open(os.path.join(a.workdir, "text", f"{pmid}.txt"), "w").write(txt)
                    rec.update(source=src, chars=len(txt)); n_ok += 1
                else:
                    rec["status"] = f"{src}:short({len(txt)})"
            except urllib.error.HTTPError as e:
                rec["status"] = f"http{e.code}"
            except Exception as e:  # noqa: BLE001
                rec["status"] = type(e).__name__
            out.write(json.dumps(rec) + "\n"); out.flush()
            if i % 25 == 0: print(f"tdm: {i}/{len(todo)} ok={n_ok}", flush=True)
            time.sleep(a.delay)
    print(f"tdm: finished, {n_ok} full texts", flush=True)


if __name__ == "__main__":
    main()
