#!/usr/bin/env python3
"""Look up legal full-text availability for PubMed IDs (no downloads).

  OpenAlex (includes Unpaywall's OA data): oa_status, best OA location (pdf/landing), PMCID.
  Europe PMC: inEPMC (full text in Europe PMC), isOpenAccess, author manuscript, PMCID.

  python3 oa_lookup.py --pmids pmids.txt --out oa.jsonl     (resumable; rate-limited)
"""
import argparse, json, os, time, urllib.parse, urllib.request


def get(url):
    for i in range(5):
        try:
            return json.load(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "pavs-oa-lookup/1.0"}), timeout=90))
        except Exception:  # noqa: BLE001
            time.sleep(3 * (i + 1))
    return None


def openalex(batch):
    url = "https://api.openalex.org/works?" + urllib.parse.urlencode({
        "filter": "ids.pmid:" + "|".join(batch), "per-page": 100,
        "select": "ids,doi,open_access,best_oa_location,primary_location,type"})
    d = get(url) or {}
    out = {}
    for w in d.get("results", []):
        pm = (w.get("ids") or {}).get("pmid", "").rsplit("/", 1)[-1]
        b = w.get("best_oa_location") or {}
        out[pm] = {"doi": (w.get("doi") or "").replace("https://doi.org/", ""), "oa_status": (w.get("open_access") or {}).get("oa_status"),
                   "is_oa": (w.get("open_access") or {}).get("is_oa"), "oa_pdf": b.get("pdf_url"), "oa_landing": b.get("landing_page_url"),
                   "oa_host": ((b.get("source") or {}).get("host_organization_name")), "oa_version": b.get("version"),
                   "pmcid": (w.get("ids") or {}).get("pmcid", "").rsplit("/", 1)[-1] or None,
                   "publisher": (((w.get("primary_location") or {}).get("source") or {}).get("host_organization_name"))}
    return out


def europepmc(batch):
    q = "(" + " OR ".join(f"EXT_ID:{p}" for p in batch) + ") AND SRC:MED"
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q, "format": "json", "resultType": "lite", "pageSize": 100})
    d = get(url) or {}
    out = {}
    for r in (d.get("resultList") or {}).get("result", []):
        out[r.get("pmid")] = {"epmc_inEPMC": r.get("inEPMC"), "epmc_isOA": r.get("isOpenAccess"), "epmc_pmcid": r.get("pmcid"),
                              "epmc_hasPDF": r.get("hasPDF"), "epmc_authMan": r.get("authMan")}
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pmids", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    pm = [p.strip() for p in open(a.pmids) if p.strip()]
    done = {json.loads(l)["pmid"] for l in open(a.out)} if os.path.exists(a.out) else set()
    todo = [p for p in pm if p not in done]
    print(f"oa_lookup: {len(todo)} to look up", flush=True)
    with open(a.out, "a") as f:
        for i in range(0, len(todo), 50):
            b = todo[i:i + 50]
            oa = openalex(b); time.sleep(0.15)
            ep = europepmc(b); time.sleep(0.15)
            for p in b:
                f.write(json.dumps({"pmid": p, **oa.get(p, {"openalex": "not_found"}), **ep.get(p, {})}) + "\n")
            f.flush()
            if (i // 50) % 100 == 0: print(f"  {i + len(b)}/{len(todo)}", flush=True)
    print("done")


if __name__ == "__main__":
    main()
