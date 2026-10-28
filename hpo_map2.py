#!/usr/bin/env python3
"""HPO normalisation with dense retrieval + LLM choice (RAG-HPO design, Garcia et al. 2025),
ported into our pipeline.

  index     --obo hp.obo [--addons HPO_addons.csv] --model M --out IDX.npz
            entries: label, EXACT/RELATED/NARROW synonyms (BROAD excluded: matching one would map
            a phrase to a term more specific than the phrase), first sentence of the definition,
            and RAG-HPO's add-on phrases (MIT). Embeddings L2-normalised (float16).
  eval-retrieval --idx IDX.npz [--idx ...] --audit-dir D --termdb T
            recall@k of the audited gold term (Sonnet judgement) for phrase text and Qwen label.
  candidates --idx IDX.npz --phrases P.jsonl --obo hp.obo --out C.jsonl [--k 15]
            P.jsonl lines {"key", "text", "label"}; writes, per key, up to ~30 candidates: dense
            top-k for text and for label, exact label/synonym hits, and the parents of the top
            dense hits (so the LLM can back off to a more general term).
"""
import argparse, collections, csv, json, os, re, sys
import numpy as np


def load_obo(path):
    T, cur = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line == "[Term]": cur = {"syn": [], "is_a": [], "def": "", "obs": False}
        elif line.startswith("[Typedef]"): cur = None
        elif cur is not None:
            if line.startswith("id: "): cur["id"] = line[4:]; T[cur["id"]] = cur
            elif line.startswith("name: "): cur["name"] = line[6:]
            elif line.startswith("def: "): cur["def"] = line[5:].split('" [')[0].strip('"')
            elif line.startswith("synonym: "):
                kind = line.split('"')[2].split()[0]
                cur["syn"].append((line.split('"')[1], kind))
            elif line.startswith("is_a: "): cur["is_a"].append(line[6:].split("!")[0].strip())
            elif line.startswith("is_obsolete: true"): cur["obs"] = True
    return {k: v for k, v in T.items() if k.startswith("HP:") and not v["obs"] and "name" in v}


def model(name):
    from sentence_transformers import SentenceTransformer, models
    if "SapBERT" in name:
        tr = models.Transformer(name, max_seq_length=64)
        return SentenceTransformer(modules=[tr, models.Pooling(tr.get_word_embedding_dimension(), pooling_mode="cls")])
    return SentenceTransformer(name)


def embed(m, texts, bs=256):
    return m.encode(texts, batch_size=bs, normalize_embeddings=True, show_progress_bar=False, convert_to_numpy=True).astype(np.float16)


def cmd_index(a):
    T = load_obo(a.obo); ids, txt, kind = [], [], []
    def add(h, s, k):
        s = (s or "").strip()
        if s: ids.append(h); txt.append(s); kind.append(k)
    for h, t in T.items():
        add(h, t["name"], "label")
        for s, k in t["syn"]:
            if k != "BROAD": add(h, s, "syn_" + k.lower())
        if t["def"]: add(h, re.split(r"(?<=[.;])\s", t["def"])[0][:300], "def")
    if a.addons:
        for r in csv.DictReader(open(a.addons)):
            if r["HP_ID"] in T: add(r["HP_ID"], r["info"][:300], "addon")
    print(f"index entries {len(txt)} ({collections.Counter(kind)})", flush=True)
    E = embed(model(a.model), txt)
    np.savez_compressed(a.out, E=E, ids=np.array(ids), txt=np.array(txt), kind=np.array(kind), model=a.model)
    print("saved", a.out, E.shape)


class Index:
    def __init__(self, path):
        d = np.load(path, allow_pickle=False)
        self.E = d["E"].astype(np.float32); self.ids = d["ids"]; self.txt = d["txt"]; self.name = str(d["model"])
        self.m = None
        self.uid, self.inv = np.unique(self.ids, return_inverse=True)

    def query(self, texts, k):
        if self.m is None: self.m = model(self.name)
        Q = embed(self.m, texts).astype(np.float32)
        out = []
        for i in range(0, len(Q), 512):
            S = Q[i:i + 512] @ self.E.T
            best = np.full((S.shape[0], len(self.uid)), -1.0, np.float32)
            for r in range(S.shape[0]):   # best score per HPO id
                np.maximum.at(best[r], self.inv, S[r])
            top = np.argsort(-best, axis=1)[:, :k]
            out += [[(str(self.uid[j]), float(best[r, j])) for j in top[r]] for r in range(len(top))]
        return out


def audit_items(audit_dir):
    src = {}
    for r in ("EASIA", "EUR", "MENA"):
        for c in json.load(open(os.path.join(audit_dir, f"audit_{r}.json"))):
            src[c["uid"]] = (r, c["qwen_observed"])
    items = []
    for k in (1, 2, 3):
        for c in json.load(open(os.path.join(audit_dir, f"batch{k}_judged.json"))):
            reg, ph = src[c["uid"]]; zero = min(p["n"] for p in c["phrases"]) == 0
            for p in c["phrases"]:
                if p["extraction"] != "correct": continue
                s = next(x for x in ph if x["n"] == p["n"] + (1 if zero else 0))
                gold = s["hpo_id"] if p["mapping"] == "exact" else p.get("best_hpo")
                items.append({"uid": c["uid"], "region": reg, "text": s["text"] or "", "label": s["qwen_label"] or "",
                              "old": s["hpo_id"], "old_route": s["route"], "gold": gold, "judged": p["mapping"]})
    return items


def ancestors(termdb):
    anc = collections.defaultdict(set)
    for l in open(os.path.join(termdb, "graph_path.txt")):
        f = l.split("\t")
        if f[1].startswith("HP:") and f[2].startswith("HP:"): anc[f[2]].add(f[1])
    return anc


def grade(pred, gold, anc):
    if pred is None: return "none"
    if pred == gold: return "exact"
    if pred in anc.get(gold, ()): return "broader"
    if gold in anc.get(pred, ()): return "too_specific"
    return "wrong"


def cmd_eval_retrieval(a):
    items = [x for x in audit_items(a.audit_dir) if x["gold"]]
    print(f"audited phrases with a gold term: {len(items)}")
    for p in a.idx:
        I = Index(p); kk = [1, 5, 10, 15, 30]
        rt = I.query([x["text"] for x in items], max(kk)); rl = I.query([x["label"] for x in items], max(kk))
        line = []
        for k in kk:
            ht = np.mean([x["gold"] in {h for h, _ in r[:k]} for x, r in zip(items, rt)])
            hb = np.mean([x["gold"] in {h for h, _ in r1[:k]} | {h for h, _ in r2[:k]} for x, r1, r2 in zip(items, rt, rl)])
            line.append(f"@{k}: text {ht:.2f} text+label {hb:.2f}")
        print(f"{I.name}: " + "; ".join(line))


def cmd_candidates(a):
    T = load_obo(a.obo); I = Index(a.idx)
    norm = lambda s: re.sub(r"\s+", " ", re.sub(r"[^a-z0-9 ]", " ", (s or "").lower())).strip()
    lex = {}
    for h, t in T.items():
        for s in [t["name"]] + [x for x, k in t["syn"] if k == "EXACT"]: lex.setdefault(norm(s), h)
    P = [json.loads(l) for l in open(a.phrases)]
    rt = I.query([p["text"] or p["label"] for p in P], a.k); rl = I.query([p["label"] or p["text"] for p in P], a.k)
    with open(a.out, "w") as f:
        for p, r1, r2 in zip(P, rt, rl):
            score = collections.defaultdict(float)
            for h, s in r1 + r2: score[h] = max(score[h], s)
            for s in (p["text"], p["label"]):
                if norm(s) in lex: score[lex[norm(s)]] = max(score[lex[norm(s)]], 1.0)
            top = sorted(score, key=lambda h: -score[h])[:a.k + 5]
            par = [q for h in top[:5] for q in T[h]["is_a"] if q in T and q not in score][:6]
            cands = [{"id": h, "name": T[h]["name"], "parent": "; ".join(T[q]["name"] for q in T[h]["is_a"] if q in T)[:120],
                      "def": T[h]["def"][:160], "score": round(score.get(h, 0.0), 3)} for h in top + par]
            text_exact = lex.get(norm(p["text"]))
            f.write(json.dumps({"key": p["key"], "text": p["text"], "label": p["label"], "text_exact": text_exact, "cands": cands}) + "\n")
    print("wrote", a.out, len(P))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["index", "eval-retrieval", "candidates"])
    ap.add_argument("--obo"); ap.add_argument("--addons"); ap.add_argument("--model"); ap.add_argument("--out")
    ap.add_argument("--idx", action="append"); ap.add_argument("--audit-dir"); ap.add_argument("--termdb")
    ap.add_argument("--phrases"); ap.add_argument("--k", type=int, default=15)
    a = ap.parse_args()
    if a.cmd == "candidates": a.idx = a.idx[0]
    {"index": cmd_index, "eval-retrieval": cmd_eval_retrieval, "candidates": cmd_candidates}[a.cmd](a)


if __name__ == "__main__":
    main()
