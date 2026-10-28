#!/usr/bin/env python3
"""Turn Qwen per-patient extractions from PubMed abstracts (qwen_tasks.py abstract) into a case table.

- keeps patients with exactly one causal gene after HGNC normalisation;
- phenotype phrases are mapped with the hpo_map2 selection (qwen_tasks.py select), key
  md5(text.lower()+"||"+label.lower())[:16], as for PMC-Patients (pmcp_to_cases.py --map2);
- region: ancestry region if ancestry is stated (country, else ethnicity words), otherwise the
  region of the place of care; author region from the first-author affiliation country;
- setting: ancestry vs place of care, same groups as cases_pmcp_setting.tsv;
- sex, age and age at onset: age strings normalised to years where a single value is given and
  to bands (prenatal, neonate, infant, child, adolescent, adult) where the value or words allow.
Usage: abs_to_cases.py --extract 'results/abs/abs_*.jsonl' --select 'results/abs/sel/sel_*.jsonl,results/pmcp/hpo_map2_select.jsonl' \
         --abstracts abstracts.jsonl --hgnc hgnc_complete_set.txt --out results/cases_abs.tsv
"""
import argparse, csv, glob, hashlib, json, os, re, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from origin_analyses import region_of
from pmcp_to_cases import ETH, load_hgnc

UNIT = {"y": 1.0, "yr": 1.0, "yrs": 1.0, "year": 1.0, "years": 1.0, "mo": 1 / 12, "mos": 1 / 12, "month": 1 / 12,
        "months": 1 / 12, "m": 1 / 12, "wk": 1 / 52.18, "wks": 1 / 52.18, "week": 1 / 52.18, "weeks": 1 / 52.18,
        "w": 1 / 52.18, "d": 1 / 365.25, "day": 1 / 365.25, "days": 1 / 365.25, "h": 1 / 8766, "hour": 1 / 8766,
        "hours": 1 / 8766}
NUMW = {w: i for i, w in enumerate("zero one two three four five six seven eight nine ten eleven twelve thirteen "
                                   "fourteen fifteen sixteen seventeen eighteen nineteen twenty".split())}
WORDS = [("prenatal", r"fetus|fetal|foetal|prenatal|antenatal|in utero|gestation|stillb"),
         ("neonate", r"neonat|newborn|new-born|at birth|^birth$|congenital|days? (?:\d+ )?of life|first week|perinatal"),
         ("infant", r"(?:first|second|third|\d+(?:st|nd|rd|th)?) months? of life|first year"),
         ("infant", r"infan|toddler"),
         ("adolescent", r"adolescen|teen|juvenile"),
         ("child", r"child|pediatric|paediatric|school"),
         ("adult", r"adult|decade|elderly|middle-aged|late-onset|late onset")]
DECADE = {"second": 15, "third": 25, "fourth": 35, "fifth": 45, "sixth": 55, "seventh": 65, "eighth": 75}


def band(y):
    if y is None: return None
    if y < 28 / 365.25: return "neonate"
    if y < 1: return "infant"
    if y < 12: return "child"
    if y < 18: return "adolescent"
    return "adult"


def parse_age(s):
    """Return (years or None, band or None) for a free-text age."""
    if s is None: return None, None
    t = str(s).lower().strip()
    if t in ("", "none", "null", "unknown", "not stated", "not reported"): return None, None
    if re.search(r"gestation|fetus|fetal|in utero|prenatal|antenatal", t): return None, "prenatal"
    for w, v in DECADE.items():
        if w + " decade" in t: return None, band(v)
    t2 = re.sub(r"\b(" + "|".join(NUMW) + r")\b", lambda m: str(NUMW[m.group(1)]), t)
    t2 = re.sub(r"(\d+)\s+(\d)/(\d)", lambda m: str(int(m.group(1)) + int(m.group(2)) / int(m.group(3))), t2)
    # compound "3 years and 9 months", "2 y, 7 mo", "10-year-2-month-old"
    m = re.fullmatch(r"\D*?(\d+(?:\.\d+)?)[\s-]*(years?|yrs?|y)[\s,-]*(?:and\s*)?(\d+(?:\.\d+)?)[\s-]*(months?|mos?|m)\b.*", t2)
    if m: y = float(m.group(1)) + float(m.group(3)) / 12; return y, band(y)
    vals = [(float(n), UNIT[u]) for n, u in re.findall(r"(\d+(?:\.\d+)?)\s*-?\s*(years?|yrs?|y|months?|mos?|m|weeks?|wks?|wk|w|days?|d|hours?|h)\b", t2)]
    nums = re.findall(r"\d+(?:\.\d+)?", t2)
    if vals and len(nums) == len(vals):
        ys = [n * u for n, u in vals]
        if len(ys) == 1 and not re.search(r"before|after|over|under|>|<|mean|median|range|by age|up to|within", t2):
            return ys[0], band(ys[0])
        bs = {band(y) for y in ys}
        return None, (bs.pop() if len(bs) == 1 and not re.search(r"before|after|over|under|>|<|up to", t2) else None)
    if len(nums) == 1 and re.fullmatch(r"\d+(\.\d+)?", t2.strip()):
        y = float(nums[0]); return y, band(y)
    if len(nums) > 1 and re.search(r"\d\s*(?:-|to|and|,)\s*\d", t2):
        unit = re.search(r"(years?|yrs?|months?|weeks?|days?)\b", t2)
        if unit:
            u = UNIT[unit.group(1)]; bs = {band(float(n) * u) for n in nums}
            if len(bs) == 1 and not re.search(r"before|after|over|under|>|<|up to|mean|median", t2): return None, bs.pop()
        return None, None
    for b, rx in WORDS:
        if re.search(rx, t): return None, b
    return None, None


def anc_region(p):
    if not p.get("ancestry"): return "UNKNOWN"
    r = region_of(p.get("ancestry_country") or "")
    if r != "UNKNOWN": return r
    e = (p.get("ancestry") or "").lower()
    for reg, words in ETH.items():
        if any(re.search(r"\b" + w + r"\b", e) for w in words): return reg
    return "UNKNOWN"


def setting(anc, care):
    if anc == "MENA" and care == "MENA": return "MENAanc_MENAcare"
    if anc == "MENA" and care == "EUR": return "MENAanc_EURcare"
    if anc == "EASIA" and care == "EASIA": return "EASIAanc_EASIAcare"
    if anc == "EASIA" and care == "EUR": return "EASIAanc_EURcare"
    if anc == "UNKNOWN" and care == "MENA": return "local_MENAcare"
    if anc in ("EUR", "UNKNOWN") and care == "EUR": return "EURorlocal_EURcare"
    return "NA"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract", required=True); ap.add_argument("--select", required=True)
    ap.add_argument("--abstracts", required=True); ap.add_argument("--hgnc", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    hg = load_hgnc(a.hgnc)
    m2 = {}
    for f in [g for pat in a.select.split(",") for g in glob.glob(pat)]:
        for l in open(f):
            d = json.loads(l)
            if d.get("ok"): m2[d["key"]] = d.get("hpo")
    meta = {}
    for l in open(a.abstracts):
        d = json.loads(l); meta[d["pmid"]] = d.get("first_country") or ""
    ex = {}
    for f in sorted(glob.glob(a.extract)):
        for l in open(f):
            d = json.loads(l)
            if d["pmid"] not in ex or d.get("ok"): ex[d["pmid"]] = d
    st = collections.Counter()
    cols = ["id", "pmid", "gene", "region", "region_source", "author_region", "corpus", "n_phrases", "n_mapped",
            "consanguinity", "setting", "sex", "age", "age_years", "age_band", "onset", "onset_years", "onset_band", "terms"]
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(cols)
        for pmid in sorted(ex, key=int):
            d = ex[pmid]
            if not d.get("ok"): st["failed"] += 1; continue
            for i, p in enumerate((d["result"] or {}).get("patients") or [], 1):
                st["patients"] += 1
                genes = {hg.get((g.get("symbol") or "").upper().strip()) for g in p.get("causal_genes") or [] if isinstance(g, dict)} - {None}
                if len(genes) != 1: st["gene_not_unique"] += 1; continue
                ts, obs = [], p.get("observed_phenotypes") or []
                for q in obs:
                    if isinstance(q, str): q = {"text": q, "hpo_label": q}
                    if not isinstance(q, dict): continue
                    k = hashlib.md5(((q.get("text") or "").strip().lower() + "||" + (q.get("hpo_label") or "").strip().lower()).encode()).hexdigest()[:16]
                    h = m2.get(k); st["phr"] += 1; st["phr_mapped"] += bool(h)
                    if h and h not in ts: ts.append(h)
                if not ts: st["no_terms"] += 1; continue
                anc = anc_region(p); care = region_of(p.get("care_country") or "")
                reg, src = (anc, "patient") if anc != "UNKNOWN" else (care, "care")
                ay, ab = parse_age(p.get("age")); oy, ob = parse_age(p.get("age_at_onset"))
                sex = p.get("sex") if p.get("sex") in ("male", "female") else "unknown"
                w.writerow([f"ABS_{pmid}-{i}", pmid, genes.pop(), reg, src, region_of(meta.get(pmid, "")), "pubmed_abs",
                            len(obs), len(ts), p.get("consanguinity") or "", setting(anc, care), sex,
                            p.get("age") or "", "" if ay is None else f"{ay:.3g}", ab or "",
                            p.get("age_at_onset") or "", "" if oy is None else f"{oy:.3g}", ob or "", ";".join(ts)])
                st["cases"] += 1; st["reg_" + reg] += 1; st["set_" + setting(anc, care)] += 1
                st["sex_" + sex] += 1; st["ageband_" + (ab or "none")] += 1; st["agey"] += ay is not None
    print(json.dumps(dict(sorted(st.items())), indent=0))
    print(f"phrases mapped {st['phr_mapped'] / max(st['phr'], 1):.1%}")


if __name__ == "__main__":
    main()
