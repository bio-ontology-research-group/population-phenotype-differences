#!/usr/bin/env python3
"""Gene universes of the external tools and the paired cases whose causative gene is outside them (tool_restriction).

  lirical   genes (NCBI Gene id) linked by LIRICAL's data/mim2gene_medgen (rows of type 'phenotype' with a GeneID) to an
            OMIM disease annotated in LIRICAL's data/phenotype.hpoa, i.e. every gene LIRICAL can give a disease-based
            post-test probability (all of them are scored in every run).
  exomiser  genes that Exomiser 15.1.0 enumerates in phenotype-only mode: the rows of a TSV_GENE output (--exomiser-tsv),
            which lists every gene of the hg19 2512 Ensembl transcript set with its phenotype score (0 for most);
            matched by NCBI Gene id, symbol otherwise (some rows have no NCBI Gene id).

Case genes are matched by NCBI Gene id (symbol -> id via genes_to_phenotype.txt), symbol otherwise. Writes
OUT_not_in_universe.csv (tool, id, gene, region) and OUT_universe.json (sizes and definitions), and prints counts per
region. Also counts profile terms that are not current classes in each tool's ontology copy (hp.json / hp.obo).
"""
import argparse, collections, csv, glob, json


def obo_ids(path):
    ids, alt, obs = set(), {}, set(); cur = None
    for line in open(path):
        line = line.strip()
        if line == "[Term]": cur = {"alt": []}
        elif cur is not None:
            if line.startswith("id: "): cur["id"] = line[4:]
            elif line.startswith("alt_id: "): cur["alt"].append(line[8:])
            elif line.startswith("is_obsolete: true"): cur["obs"] = True
            elif line == "":
                if "id" in cur:
                    (obs if cur.get("obs") else ids).add(cur["id"])
                    for x in cur["alt"]: alt[x] = cur["id"]
                cur = None
    return ids, alt, obs


def json_ids(path):
    d = json.load(open(path)); ids, alt, obs = set(), {}, set()
    for n in d["graphs"][0]["nodes"]:
        i = n["id"].rsplit("/", 1)[-1].replace("_", ":")
        if not i.startswith("HP:"): continue
        m = n.get("meta", {})
        if m.get("deprecated"): obs.add(i)
        else: ids.add(i)
        for bv in m.get("basicPropertyValues", []):
            if bv.get("pred", "").endswith("hasAlternativeId"): alt[bv["val"]] = i
    return ids, alt, obs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles", required=True); ap.add_argument("--g2p", required=True)
    ap.add_argument("--lirical-data", required=True); ap.add_argument("--exomiser-tsv", required=True)
    ap.add_argument("--exomiser-hpo", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sym2id = {}
    for line in open(a.g2p):
        if line.startswith("ncbi_gene_id"): continue
        c = line.rstrip("\n").split("\t"); sym2id[c[1]] = c[0]
    # LIRICAL
    ent2sym = {}
    for r in csv.DictReader(open(f"{a.lirical_data}/hgnc_complete_set.txt"), delimiter="\t"):
        if r.get("entrez_id"): ent2sym[r["entrez_id"]] = r["symbol"]
    dis = set()
    for line in open(f"{a.lirical_data}/phenotype.hpoa"):
        if line.startswith("#") or line.startswith("database_id"): continue
        dis.add(line.split("\t", 1)[0])
    lir = set()
    for line in open(f"{a.lirical_data}/mim2gene_medgen"):
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if c[2] == "phenotype" and c[1] != "-" and "OMIM:" + c[0] in dis: lir.add(c[1])
    U = {"lirical": (lir, {ent2sym.get(e, "") for e in lir} - {""})}
    # Exomiser
    ex_id, ex_sym = set(), set()
    files = sorted(glob.glob(a.exomiser_tsv)); n_ex_rows = 0
    for f in files:
        for r in csv.DictReader(open(f), delimiter="\t"):
            n_ex_rows += 1
            if r.get("ENTREZ_GENE_ID"): ex_id.add(r["ENTREZ_GENE_ID"])
            ex_sym.add(r["GENE_SYMBOL"])
    U["exomiser"] = (ex_id, ex_sym)
    cases = [json.loads(l) for l in open(a.profiles)]
    rows = []; cnt = collections.defaultdict(collections.Counter)
    for t, (ids, syms) in U.items():
        for c in cases:
            g = c["gene"]; e = sym2id.get(g)
            if not ((e and e in ids) or g in syms):
                rows.append({"tool": t, "id": c["id"], "gene": g, "region": c["region"]}); cnt[t][c["region"]] += 1
    with open(a.out + "_not_in_universe.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["tool", "id", "gene", "region"]); w.writeheader(); w.writerows(rows)
    # ontology coverage of the profile terms
    onto = {"lirical": json_ids(f"{a.lirical_data}/hp.json"), "exomiser": obo_ids(a.exomiser_hpo)}
    term_stats = {}
    for t, (ids, alt, obs) in onto.items():
        allt = collections.Counter(x for c in cases for x in c["terms_all"])
        miss = {x: n for x, n in allt.items() if x not in ids}
        term_stats[t] = {"distinct_terms": len(allt), "not_current_class": len(miss), "of_which_alt_id": sum(1 for x in miss if x in alt),
                         "of_which_obsolete": sum(1 for x in miss if x in obs), "case_term_occurrences_not_current": sum(miss.values()),
                         "cases_with_any_not_current": sum(1 for c in cases if any(x in miss for x in c["terms_all"]))}
    info = {"lirical": {"N_universe": len(lir), "definition": "NCBI Gene ids linked by data/mim2gene_medgen (type=phenotype, GeneID present) to an OMIM disease in data/phenotype.hpoa"},
            "exomiser": {"N_universe": n_ex_rows, "N_symbols": len(ex_sym), "N_rows_with_entrez": len(ex_id),
                         "definition": f"genes (rows) of {len(files)} phenotype-only TSV_GENE output(s): {', '.join(files)}"},
            "not_in_universe": {t: dict(cnt[t]) | {"total": sum(cnt[t].values())} for t in U}, "profile_terms_vs_tool_ontology": term_stats,
            "n_cases": len(cases)}
    json.dump(info, open(a.out + "_universe.json", "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
