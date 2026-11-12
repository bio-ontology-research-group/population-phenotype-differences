#!/usr/bin/env python3
"""Rank the causative gene of paired cases with Exomiser (phenotype-only, hiPhive) and LIRICAL (phenotype-only) under
three HPO profiles (all findings, bedside-only, random subset of the bedside size) exported by export_profiles.py.

  exomiser  Writes one v1 phenopacket (YAML) per case x profile, runs one Exomiser JVM per sub-batch of cases with the
            `batch` command over both analysis configurations (hiphive_default, hiphive_human), output PARQUET, and
            reads the per-gene phenotype score (genePhenotypeScore = EXOMISER_GENE_PHENO_SCORE at full precision).
  lirical   Runs `lirical prioritize` (no VCF, no --assembly = phenotype-only) once per case x profile (LIRICAL has no
            batch mode), JSON output, and gives each gene the highest post-test probability over the diseases it is
            associated with in LIRICAL's own gene-disease file (data/mim2gene_medgen, OMIM phenotype rows with a GeneID).

Rank of the causative gene: midrank among all scored genes by descending score; if the gene is not scored, the midrank
of the unscored remainder of the tool's gene universe, absent=1. Genes matched by NCBI Gene id (case gene symbol ->
id via genes_to_phenotype.txt), symbol otherwise. Raw tool outputs live in the work dir (node-local scratch) and are
deleted after parsing.

  tool_restriction_run.py --tool exomiser --profiles P.jsonl --ids IDS.txt --chunks 150 --chunk k --work $TMPDIR/x --out OUT.csv ...
"""
import argparse, csv, json, os, shutil, subprocess, sys, time
from collections import defaultdict

PROFILES = ("all", "bed", "rand")
FIELDS = ["id", "gene", "region", "corpus", "n_terms", "n_bed", "tool", "config", "profile", "rank", "score", "n_scored", "absent"]


def load_case_gene_ids(g2p):
    sym2id = {}
    for line in open(g2p):
        if line.startswith("ncbi_gene_id") or line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        sym2id[c[1]] = c[0].replace("NCBIGene:", "")
    return sym2id


def midrank(scores, gid, gsym, universe_n):
    """scores: dict key -> score where key = (entrez or None, symbol). Returns rank, score, n_scored, absent."""
    n = len(scores)
    hit = [v for (e, s), v in scores.items() if gid and e == gid]  # NCBI Gene id first
    if not hit: hit = [v for (e, s), v in scores.items() if s == gsym]  # symbol otherwise
    s0 = max(hit) if hit else None
    if s0 is None:
        return n + (universe_n - n + 1) / 2, "", n, 1
    gt = sum(1 for v in scores.values() if v > s0); eq = sum(1 for v in scores.values() if v == s0)
    return gt + (eq + 1) / 2, s0, n, 0


def load_cases(a):
    ids = set(x.strip() for x in open(a.ids) if x.strip()) if a.ids else None
    cases = [json.loads(l) for l in open(a.profiles)]
    if ids is not None: cases = [c for c in cases if c["id"] in ids]
    return cases[a.chunk::a.chunks]


def safe(s):
    return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in s)


def run_exomiser(a, cases, sym2id, writer):
    """One Exomiser JVM (`batch` command) per sub-batch of --exo-batch cases; PARQUET output, which holds the genes with a
    phenotype score > 0 at full precision (genePhenotypeScore); every other gene of the universe (all genes Exomiser
    enumerates; --universe-n = the row count of a TSV_GENE output, which lists them all, with score 0) scores 0 and gets
    the midrank of the zero-score block (absent=1)."""
    import pyarrow.parquet as pq
    configs = {"hiphive_default": a.exo_default, "hiphive_human": a.exo_human}
    nfail = 0; tt = 0.0; nan = 0; ex_written = False
    for b0 in range(0, len(cases), a.exo_batch):
        work = os.path.join(a.work, f"b{b0}"); os.makedirs(work, exist_ok=True)
        lines = []; jobs = []
        for k, c in enumerate(cases[b0:b0 + a.exo_batch], start=b0):
            for p in PROFILES:
                sp = os.path.join(work, f"s{k}_{p}.yml")
                with open(sp, "w") as f:
                    f.write(f"---\nid: s{k}_{p}\nsubject:\n  id: s{k}_{p}\nphenotypicFeatures:\n")
                    for t in c[f"terms_{p}"]: f.write(f"  - type:\n      id: {t}\n")
                    f.write("metaData:\n  phenopacketSchemaVersion: 1.0\n")
                for cfg, yml in configs.items():
                    fn = f"r{k}_{p}_{cfg}"
                    lines.append(f"--analysis {yml} --sample {sp} --output-format PARQUET --output-directory {work}/out --output-filename {fn}")
                    jobs.append((c, p, cfg, os.path.join(work, "out", fn + ".parquet")))
        bf = os.path.join(work, "batch.txt"); open(bf, "w").write("\n".join(lines) + "\n")
        cmd = ["java", f"-Xmx{a.xmx}", "-XX:ActiveProcessorCount=" + str(a.cpus), f"-Djava.io.tmpdir={work}",
               f"-Dspring.config.location={a.exo_props}", "-jar", a.exo_jar, "batch", bf]
        if b0 == 0: print(" ".join(cmd), flush=True); print("example batch line:", lines[0], flush=True)
        t0 = time.time()
        r = subprocess.run(cmd, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        el = time.time() - t0; tt += el; nan += len(jobs)
        with open(os.path.join(a.logdir, f"exomiser_chunk{a.chunk}.log"), "a") as fl:
            fl.write(f"#### sub-batch {b0}: exit {r.returncode} {el:.1f}s\n")
            fl.write("\n".join(l for l in r.stdout.splitlines() if " WARN " in l or " ERROR " in l or "Exception" in l
                               or "Finished analysis" in l or "Started Main" in l) + "\n")
        print(f"sub-batch {b0}: exomiser exit {r.returncode} in {el:.1f}s for {len(jobs)} analyses", flush=True)
        for c, p, cfg, path in jobs:
            if not os.path.exists(path):
                nfail += 1; print("MISSING", c["id"], p, cfg, path, file=sys.stderr); continue
            if a.chunk == 0 and not ex_written:
                shutil.copy(path, os.path.join(a.logdir, "example_exomiser.parquet")); ex_written = True  # keep one raw output
            t = pq.read_table(path, columns=["geneSymbol", "ncbiGeneId", "genePhenotypeScore"]).to_pydict()
            scores = {}
            for sym, e, v in zip(t["geneSymbol"], t["ncbiGeneId"], t["genePhenotypeScore"]):
                if v is None or v <= 0: continue
                key = (e or None, sym)
                scores[key] = max(v, scores.get(key, v))
            rank, score, n, absent = midrank(scores, sym2id.get(c["gene"]), c["gene"], a.universe_n)
            writer.writerow({"id": c["id"], "gene": c["gene"], "region": c["region"], "corpus": c["corpus"], "n_terms": c["n_terms"],
                             "n_bed": c["n_bed"], "tool": "exomiser", "config": cfg, "profile": p, "rank": rank, "score": score if score != "" else 0.0,
                             "n_scored": n, "absent": absent})
        shutil.rmtree(work, ignore_errors=True)
    print(f"TIMING exomiser analyses={nan} wall={tt:.1f} per_analysis={tt / max(nan, 1):.2f} missing={nfail}", flush=True)
    return nfail


def load_lirical_genes(datadir):
    ent2sym = {}
    with open(os.path.join(datadir, "hgnc_complete_set.txt")) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            if r.get("entrez_id"): ent2sym[r["entrez_id"]] = r["symbol"]
    d2g = defaultdict(set)
    for line in open(os.path.join(datadir, "mim2gene_medgen")):
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if c[2] == "phenotype" and c[1] != "-":
            d2g["OMIM:" + c[0]].add((c[1], ent2sym.get(c[1], "")))
    return d2g


def run_lirical(a, cases, sym2id, writer):
    work = a.work; os.makedirs(work, exist_ok=True)
    d2g = load_lirical_genes(a.lirical_data)
    nfail = 0; tt = 0.0; nrun = 0
    for k, c in enumerate(cases):
        for p in PROFILES:
            prefix = f"l{k}_{p}"
            cmd = ["java", f"-Xmx{a.xmx}", "-XX:ActiveProcessorCount=" + str(a.cpus), "-XX:+UseSerialGC", f"-Djava.io.tmpdir={work}",
                   "-jar", a.lirical_jar, "prioritize", "-d", a.lirical_data, "-p", ",".join(c[f"terms_{p}"]),
                   "-o", work, "-x", prefix, "-f", "json", "--sample-id", prefix]
            t0 = time.time()
            r = subprocess.run(cmd, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            el = time.time() - t0; tt += el; nrun += 1
            warn = [l for l in r.stdout.splitlines() if " WARN " in l or " ERROR " in l]
            if warn:
                with open(os.path.join(a.logdir, f"lirical_warn_chunk{a.chunk}.txt"), "a") as fw:
                    fw.write(f"## {c['id']} {p}\n" + "\n".join(warn) + "\n")
            path = os.path.join(work, prefix + ".json")
            if r.returncode != 0 or not os.path.exists(path):
                nfail += 1
                print("FAIL", c["id"], p, r.returncode, r.stdout[-3000:], file=sys.stderr); continue
            if a.chunk == 0 and nrun == 1: shutil.copy(path, os.path.join(a.logdir, "example_lirical.json"))  # keep one raw output
            res = json.load(open(path))["analysisResults"]
            os.remove(path)
            scores = {}
            for x in res:
                for g in d2g.get(x["diseaseId"], ()):
                    v = x["posttestProbability"]
                    if g not in scores or v > scores[g]: scores[g] = v
            if a.universe_n is None: a.universe_n = len(scores)
            rank, score, n, absent = midrank(scores, sym2id.get(c["gene"]), c["gene"], a.universe_n)
            writer.writerow({"id": c["id"], "gene": c["gene"], "region": c["region"], "corpus": c["corpus"], "n_terms": c["n_terms"],
                             "n_bed": c["n_bed"], "tool": "lirical", "config": "lirical_pheno", "profile": p, "rank": rank,
                             "score": score, "n_scored": n, "absent": absent})
            if nrun == 1: print(" ".join(cmd), flush=True)
    print(f"TIMING lirical analyses={nrun} wall={tt:.1f} per_analysis={tt / max(nrun, 1):.2f} failed={nfail}", flush=True)
    return nfail


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tool", choices=["exomiser", "lirical"], required=True)
    ap.add_argument("--profiles", required=True); ap.add_argument("--ids"); ap.add_argument("--g2p", required=True)
    ap.add_argument("--chunks", type=int, default=1); ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--work", required=True); ap.add_argument("--out", required=True); ap.add_argument("--logdir", required=True)
    ap.add_argument("--xmx", default="4g"); ap.add_argument("--cpus", type=int, default=1)
    ap.add_argument("--universe-n", type=int, default=None)
    ap.add_argument("--exo-jar"); ap.add_argument("--exo-props"); ap.add_argument("--exo-default"); ap.add_argument("--exo-human")
    ap.add_argument("--exo-batch", type=int, default=50)
    ap.add_argument("--lirical-jar"); ap.add_argument("--lirical-data")
    a = ap.parse_args()
    if not a.universe_n: a.universe_n = None
    cases = load_cases(a)
    sym2id = load_case_gene_ids(a.g2p)
    print(f"chunk {a.chunk}/{a.chunks}: {len(cases)} cases", flush=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS); w.writeheader()
        nfail = (run_exomiser if a.tool == "exomiser" else run_lirical)(a, cases, sym2id, w)
    shutil.rmtree(a.work, ignore_errors=True)
    print(f"DONE failures={nfail}")


if __name__ == "__main__":
    main()
