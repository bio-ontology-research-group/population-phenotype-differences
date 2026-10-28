#!/usr/bin/env python3
"""Vectorised re-implementation of the PAVS gene-ranking evaluation
(evaluation/similarity/scripts/phenotype_similarity_analysis.py in the PAVS repo):
Resnik pairwise similarity with extrinsic (gene-frequency) or intrinsic IC,
symmetric Best Match Average over all genes in genes_to_phenotype.txt.

Resnik(p, g) = max IC over common ancestors, so for a gene profile G
  max_g Resnik(p, g) = max IC over anc(p) & Anc(G)
which makes both BMA directions sparse-matrix operations and lets us re-rank
thousands of cases per minute, including with per-case modified references
(leave-publication-out).
"""
import math, re, json, os
import numpy as np
import scipy.sparse as sp


def parse_obo(path):
    terms, alt = {}, {}
    cur = None
    for line in open(path):
        line = line.strip()
        if line == "[Term]":
            cur = {"id": None, "is_a": [], "obsolete": False, "alt": []}
        elif cur is not None:
            if line.startswith("id: "):
                cur["id"] = line[4:]
            elif line.startswith("is_a: "):
                cur["is_a"].append(line[6:].split("!")[0].strip())
            elif line.startswith("alt_id: "):
                cur["alt"].append(line[8:].strip())
            elif line.startswith("is_obsolete: true"):
                cur["obsolete"] = True
            elif line == "":
                if cur["id"]:
                    terms[cur["id"]] = cur
                    for a in cur["alt"]:
                        alt[a] = cur["id"]
                cur = None
    if cur is not None and cur["id"]:
        terms[cur["id"]] = cur
    return terms, alt


class Engine:
    def __init__(self, obo, g2p, ic_type="extrinsic"):
        self.terms, self.alt = parse_obo(obo)
        # the original code keeps obsolete terms as nodes too (they have no is_a); harmless
        self.tid = {t: i for i, t in enumerate(sorted(self.terms))}
        self.names = sorted(self.terms)
        n = len(self.tid)
        # ancestor closure incl. self
        anc = {}
        def get(t):
            if t in anc:
                return anc[t]
            r = {t}
            for p in self.terms[t]["is_a"]:
                if p in self.terms:
                    r |= get(p)
            anc[t] = r
            return r
        for t in self.terms:
            get(t)
        self.anc_sets = anc
        rows, cols = [], []
        for t, s in anc.items():
            i = self.tid[t]
            for a in s:
                rows.append(i); cols.append(self.tid[a])
        self.A = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(n, n))
        # gene profiles
        g2t = {}
        for line in open(g2p):
            if line.startswith("ncbi_gene_id") or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            gene, hpo = c[1], c[2]
            if hpo in self.terms:
                g2t.setdefault(gene, set()).add(hpo)
        self.genes = sorted(g2t)
        self.gid = {g: i for i, g in enumerate(self.genes)}
        self.g2t = {g: sorted(v) for g, v in g2t.items()}
        self.build_gene_matrices(self.g2t)
        self.set_ic(ic_type)

    def build_gene_matrices(self, g2t):
        rows, cols = [], []
        for g, ts in g2t.items():
            i = self.gid[g]
            for t in ts:
                rows.append(i); cols.append(self.tid[t])
        n = len(self.tid)
        GT = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(len(self.genes), n))
        self.GT = GT
        cnt = np.asarray(GT.sum(axis=1)).ravel()
        cnt[cnt == 0] = 1
        self.GTnorm = sp.diags(1.0 / cnt) @ GT  # row-normalised for mean over gene terms
        GA = (GT @ self.A)
        GA.data[:] = 1.0
        self.GA = GA.tocsc()  # genes x terms, 1 if term is ancestor of some gene term
        self.GAcsr = GA.tocsr()
        self.gene_cnt_terms = np.asarray(GT.sum(axis=1)).ravel()

    def set_ic(self, ic_type):
        n = len(self.tid)
        if ic_type == "extrinsic":
            counts = np.asarray(self.GAcsr.sum(axis=0)).ravel()
            N = len(self.genes)
            ic = np.where(counts > 0, -np.log(np.maximum(counts, 1) / N), 0.0)
        else:
            # intrinsic: -ln(desc(t)/total) with desc incl. self; A^T gives descendants
            desc = np.asarray(self.A.sum(axis=0)).ravel()
            ic = -np.log(desc / n)
        self.ic = ic.astype(np.float64)
        self.ic_type = ic_type

    def map_terms(self, hpos):
        out = set()
        for h in hpos:
            h = self.alt.get(h, h)
            if h in self.tid:
                out.add(h)
        return sorted(out)

    def bma_all_genes(self, hpos, GA=None, GTnorm=None):
        """Return BMA Resnik score vector over all genes for a term list."""
        GA = self.GA if GA is None else GA
        GTnorm = self.GTnorm if GTnorm is None else GTnorm
        P = [self.tid[h] for h in hpos]
        if not P:
            return np.zeros(len(self.genes))
        # direction 1: patient terms -> best gene ancestor
        s1 = np.zeros(len(self.genes))
        for p in P:
            a = np.array(sorted(self.tid[x] for x in self.anc_sets[self.names[p]]))
            sub = GA[:, a].toarray() * self.ic[a]
            s1 += sub.max(axis=1)
        s1 /= len(P)
        # direction 2: gene terms -> best patient ancestor
        v = np.zeros(len(self.tid))
        for p in P:
            for x in self.anc_sets[self.names[p]]:
                j = self.tid[x]
                if self.ic[j] > v[j]:
                    v[j] = self.ic[j]
        M = self.A @ sp.diags(v)
        s2t = np.asarray(M.max(axis=1).todense()).ravel()  # per term: max IC of its ancestors within Anc(P)
        s2 = GTnorm @ s2t
        return (s1 + s2) / 2.0

    def rank_gene(self, hpos, gene, **kw):
        if gene not in self.gid:
            return None, None, None
        scores = self.bma_all_genes(hpos, **kw)
        # original: sort by (-score, gene name)
        order = np.lexsort((np.array(self.genes), -scores))
        pos = np.where(order == self.gid[gene])[0][0] + 1
        truth = np.zeros(len(scores)); truth[self.gid[gene]] = 1
        # AUC = P(score_true > score_other) with ties counted half
        s = scores[self.gid[gene]]
        others = np.delete(scores, self.gid[gene])
        auc = ((others < s).sum() + 0.5 * (others == s).sum()) / len(others)
        return int(pos), float(s), float(auc)


def metrics(ranks, aucs=None):
    r = np.array(ranks, dtype=float)
    out = {"n": len(r), "hits1": float((r == 1).mean() * 100), "hits10": float((r <= 10).mean() * 100),
           "mrr": float((1 / r).mean()), "median_rank": float(np.median(r))}
    if aucs is not None:
        out["auc"] = float(np.nanmean(aucs))
    return out
