#!/usr/bin/env python3
"""Turn Qwen per-patient extractions from PMC-Patients into a case table for case_pairs.py.

- keeps patients with mendelian_germline = true and exactly one causal gene (after HGNC
  alias/previous-symbol normalisation);
- maps observed phenotype phrases to HPO: exact match of the proposed HPO label, then of the
  verbatim phrase, against HPO labels and EXACT synonyms (normalised); then a strict fuzzy
  match (difflib ratio >= 0.92) on the label; unmatched phrases are dropped and counted;
  excluded phenotypes are never used as observed terms;
- region: origin stated in the patient summary (country, else ethnicity), otherwise the
  article-level label from labels_llm-style TSV (--article-labels); author region from the
  first-author affiliation proxy (--author-labels).
Usage: pmcp_to_cases.py --patients patients.jsonl --article-labels pmcp_labels_llm.tsv \
         --author-labels pmcp_labels_proxy.tsv --hpo-dir H --hgnc hgnc_complete_set.txt --out cases_pmcp.tsv
"""
import argparse, csv, json, re, os, sys, difflib, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from origin_analyses import region_of

ETH = {"MENA": ["arab", "saudi", "emirati", "qatari", "kuwaiti", "omani", "yemeni", "jordanian", "lebanese", "syrian", "iraqi",
                "egyptian", "moroccan", "tunisian", "algerian", "libyan", "sudanese", "iranian", "persian", "turkish", "kurdish",
                "pakistani", "afghan", "palestinian", "bedouin", "middle eastern", "north african", "baloch", "pashtun", "punjabi"],
       "EASIA": ["chinese", "han", "japanese", "korean", "taiwanese", "thai", "vietnamese", "east asian", "malay", "filipino"],
       "EUR": ["caucasian", "european", "white", "german", "french", "italian", "spanish", "dutch", "british", "english", "polish",
               "swedish", "danish", "norwegian", "finnish", "belgian", "austrian", "swiss", "portuguese", "greek", "czech", "irish",
               "ashkenazi", "scottish", "hungarian", "russian"]}


def norm(t):
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9 ]", " ", (t or "").lower())).strip()


def load_hpo(obo):
    lab, cur, obs = {}, None, False
    terms = []
    for line in open(obo):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur, obs = {"syn": []}, False; terms.append(cur)
        elif cur is not None:
            if line.startswith("id: "): cur["id"] = line[4:]
            elif line.startswith("name: "): cur["name"] = line[6:]
            elif line.startswith("synonym: ") and " EXACT" in line: cur["syn"].append(line.split('"')[1])
            elif line.startswith("is_obsolete: true"): cur["obs"] = True
            elif line == "": cur = None
    for t in terms:
        if t.get("obs") or not t.get("id", "").startswith("HP:"): continue
        for s in [t.get("name", "")] + t["syn"]:
            lab.setdefault(norm(s), t["id"])
    return lab


def load_hgnc(path):
    m = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["status"] != "Approved": continue
        sym = r["symbol"]; m[sym.upper()] = sym
        for k in ("prev_symbol", "alias_symbol"):
            for x in (r.get(k) or "").strip('"').split("|"):
                if x: m.setdefault(x.upper(), sym)
    return m


def region_from_origin(o):
    if not o: return None
    r = region_of(o.get("country") or "")
    if r not in ("UNKNOWN",): return r
    e = (o.get("ethnicity_or_ancestry") or "").lower()
    for reg, words in ETH.items():
        if any(re.search(r"\b" + w + r"\b", e) for w in words): return reg
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", required=True); ap.add_argument("--article-labels"); ap.add_argument("--author-labels")
    ap.add_argument("--hpo-dir", required=True); ap.add_argument("--hgnc", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--map", help="JSONL from pmcp_extract.py map (normalised label -> HPO id)")
    ap.add_argument("--map2", help="JSONL from qwen_tasks.py select (hpo_map2 candidates): replaces all other mapping routes")
    a = ap.parse_args()
    lab = load_hpo(os.path.join(a.hpo_dir, "hp.obo")); keys = list(lab)
    qmap = {}
    if a.map:
        for l in open(a.map):
            r = json.loads(l)
            if r.get("hpo"): qmap[r["key"]] = r["hpo"]
    hg = load_hgnc(a.hgnc)
    art = {r["pmid"]: r["region"] for r in csv.DictReader(open(a.article_labels), delimiter="\t")} if a.article_labels else {}
    aut = {r["pmid"]: r["region"] for r in csv.DictReader(open(a.author_labels), delimiter="\t")} if a.author_labels else {}
    st = collections.Counter(); fuzzy_cache = {}
    m2 = {}
    if a.map2:
        import hashlib
        for l in open(a.map2):
            d = json.loads(l)
            if d.get("ok"): m2[d["key"]] = d.get("hpo")
    def map_phrase(p):
        if isinstance(p, str): p = {"text": p, "hpo_label": p}
        if not isinstance(p, dict): return None, "none"
        if a.map2:
            t, lab_ = (p.get("text") or "").strip(), (p.get("hpo_label") or "").strip()
            k = hashlib.md5((t.lower() + "||" + lab_.lower()).encode()).hexdigest()[:16]
            h = m2.get(k)
            return (h, "map2") if h else (None, "none")
        for cand in (p.get("hpo_label"), p.get("text")):
            k = norm(cand)
            if k in lab: return lab[k], "exact"
        k = norm(p.get("hpo_label"))
        if not k: return None, "none"
        if k in qmap: return qmap[k], "qwen"
        if k not in fuzzy_cache:
            m = difflib.get_close_matches(k, keys, n=1, cutoff=0.92)
            fuzzy_cache[k] = lab[m[0]] if m else None
        return fuzzy_cache[k], ("fuzzy" if fuzzy_cache[k] else "none")
    seen = set()
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "pmid", "gene", "region", "region_source", "author_region", "corpus", "n_phrases", "n_mapped", "consanguinity", "terms"])
        for line in open(a.patients):
            d = json.loads(line)
            if not d.get("ok") or d["uid"] in seen: continue
            seen.add(d["uid"]); r = d["result"] or {}
            st["patients"] += 1
            if not r.get("mendelian_germline"): continue
            st["mendelian"] += 1
            genes = {hg.get((g.get("symbol") or "").upper().strip()) for g in r.get("causal_genes") or []} - {None}
            if len(genes) != 1:
                st["gene_not_unique"] += 1; continue
            gene = genes.pop()
            ts = []
            obs = r.get("observed_phenotypes") or []
            for p in obs:
                h, how = map_phrase(p); st["phr_" + how] += 1
                if h and h not in ts: ts.append(h)
            if not ts:
                st["no_terms"] += 1; continue
            reg = region_from_origin(r.get("origin")); src = "patient"
            if not reg:
                reg = art.get(d["pmid"], "UNKNOWN"); src = "article"
            w.writerow([f"PMCP_{d['uid']}", d["pmid"], gene, reg, src, aut.get(d["pmid"], "UNKNOWN"), "pmcp", len(obs), len(ts), r.get("consanguinity", ""), ";".join(ts)])
            st["cases"] += 1
    tot = st["phr_exact"] + st["phr_fuzzy"] + st["phr_none"] + st["phr_qwen"] + st["phr_map2"]
    print(dict(st)); print(f"map2 {st['phr_map2']/max(tot,1):.1%}; phrase mapping: exact {st['phr_exact']/max(tot,1):.1%} qwen {st['phr_qwen']/max(tot,1):.1%} fuzzy {st['phr_fuzzy']/max(tot,1):.1%} unmapped {st['phr_none']/max(tot,1):.1%}")


if __name__ == "__main__":
    main()
