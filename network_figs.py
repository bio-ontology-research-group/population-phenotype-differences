#!/usr/bin/env python3
"""Network figures for the population term tests.

shift   Population-shift phenotype network (after the human symptoms-disease network style): nodes are HPO terms
        one or two levels below the organ systems that are frequent in the case collection; edges join terms that
        co-occur in patients more than expected (phi over all cases); node colour is the term's MH z for one
        population vs EUR, node shape the finding class (circle bedside/history, square test-based, diamond mixed or
        unclassified grouping term), node size the number of cases. One panel per population on a shared layout.
dag     Significant-term HPO subgraph per population: terms with BH q < Q, is_a edges between them after transitive
        reduction, drawn radially (graphviz twopi) from the phenotypic-abnormality root; colour = direction and |z|.
popnet  Population similarity network: one node per population (size = cases in the similarity set, colour = excess
        z-similarity of within-population pairs over pairs with all other populations, ring = BH q < 0.05); an arrow
        A -> B when A-A pairs are more (red) or less (blue) similar than A-B pairs at BH q < 0.1 (solid q < 0.05),
        width = |excess|. Panels: total and patient-carried.
spend   Country-level scatter: mean within-stratum residual test share vs health expenditure per capita.

  network_figs.py shift --obo hp.obo --cases A.tsv [...] --terms-dir results/v4/terms_all --termclass T --out F
  network_figs.py dag   --obo hp.obo --terms-dir results/v4/terms_all --q 0.05 --out F
  network_figs.py popnet --contrasts results/v4/allpairs/contrasts.tsv --within results/v4/qvalues.tsv
                          --cases A.tsv [...] --out F
  network_figs.py spend --points results/v4/worldbank/healthspend_country_points.csv --out F
"""
import argparse, collections, csv, json, math, os, subprocess
import numpy as np

ROOT = "HP:0000118"
POPS = ["MENA", "SASIA", "LATAM", "AFR", "EASIA"]
TEST = {"laboratory", "imaging", "physiological_test", "pathology"}


def parse_obo(path):
    terms, cur = {}, None
    for l in open(path):
        l = l.rstrip("\n")
        if l == "[Term]": cur = {"is_a": [], "name": "", "obsolete": False}; tid = None; continue
        if l.startswith("[") and l != "[Term]": cur = None; continue
        if cur is None: continue
        if l.startswith("id: "): tid = l[4:]; terms[tid] = cur
        elif l.startswith("name: "): cur["name"] = l[6:]
        elif l.startswith("is_a: "): cur["is_a"].append(l[6:].split(" ")[0])
        elif l.startswith("is_obsolete: true"): cur["obsolete"] = True
        elif l.startswith("alt_id: "): cur.setdefault("alt", []).append(l[8:])
    return terms


def ancestors(terms):
    memo = {}
    def get(t):
        if t in memo: return memo[t]
        s = {t}
        for p in terms.get(t, {}).get("is_a", []): s |= get(p)
        memo[t] = s; return s
    return get


def depth_from_root(terms):
    ch = collections.defaultdict(list)
    for t, d in terms.items():
        for p in d["is_a"]: ch[p].append(t)
    dep, q = {ROOT: 0}, [ROOT]
    while q:
        nq = []
        for t in q:
            for c in ch[t]:
                if c not in dep: dep[c] = dep[t] + 1; nq.append(c)
        q = nq
    return dep


def load_z(terms_dir, pop):
    f = os.path.join(terms_dir, f"main_{pop}_total_terms.csv")
    return {r["hpo"]: r for r in csv.DictReader(open(f))}


def short(n):
    for p in ("Abnormality of the ", "Abnormality of ", "Abnormal "):
        if n.startswith(p): n = n[len(p):]; break
    return n[:1].upper() + n[1:]


def cmd_shift(a):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, networkx as nx
    terms = parse_obo(a.obo); anc = ancestors(terms); dep = depth_from_root(terms)
    alt = {x: t for t, d in terms.items() for x in d.get("alt", [])}
    Z = {p: load_z(a.terms_dir, p) for p in POPS}
    tc = {x["id"]: x["class"] for x in map(json.loads, open(a.termclass))}
    organ = {t: [b for b in terms if ROOT in terms[b]["is_a"]] for t in ()}
    branches = sorted(b for b in terms if ROOT in terms[b]["is_a"] and not terms[b]["obsolete"])
    # case term sets (propagated) for frequency and co-occurrence
    sets = []
    for path in a.cases:
        for r in csv.DictReader(open(path), delimiter="\t"):
            ts = [alt.get(t, t) for t in r["terms"].split(";") if alt.get(t, t) in terms]
            if ts: sets.append(set().union(*(anc(t) for t in ts)))
    N = len(sets); freq = collections.Counter(t for s in sets for t in s)
    cand = [t for t in freq if dep.get(t) in (2, 3) and all(t in Z[p] for p in ("MENA", "EASIA")) and not terms[t]["obsolete"]]
    cand.sort(key=lambda t: -freq[t]); nodes = cand[:a.n_nodes]
    nid = {t: i for i, t in enumerate(nodes)}
    X = np.zeros((N, len(nodes)), bool)
    for i, s in enumerate(sets):
        for t in s & set(nodes): X[i, nid[t]] = True
    Xf = X.astype(np.float64); co = Xf.T @ Xf; f = Xf.mean(0)
    phi = (co / N - np.outer(f, f)) / np.sqrt(np.outer(f * (1 - f), f * (1 - f)) + 1e-12)
    G = nx.Graph()
    for t in nodes:
        br = next((b for b in branches if b in anc(t)), ROOT)
        # finding class of a grouping term: class of the majority of its observed descendants in the cases
        G.add_node(t, freq=freq[t], branch=br)
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            if nodes[i] in anc(nodes[j]) or nodes[j] in anc(nodes[i]): continue
            if phi[i, j] > a.phi: G.add_edge(nodes[i], nodes[j], w=phi[i, j])
    # class per node from the observed terms below it
    below = collections.defaultdict(collections.Counter)
    for t, c in tc.items():
        if t in terms:
            for x in anc(t):
                if x in nid: below[x]["test" if c in TEST else "bed" if c in ("clinical_exam", "symptom_history") else "other"] += 1
    def shape(t):
        c = below[t]; tot = c["test"] + c["bed"]
        if tot == 0: return "D"
        return "s" if c["test"] / tot > 0.6 else "o" if c["bed"] / tot > 0.6 else "D"
    pos = nx.spring_layout(G, weight="w", seed=7, k=1.6 / math.sqrt(len(nodes)), iterations=300)
    if a.print_layout:   # pull the few far outlying nodes inward so the dense core fills the panel (cosmetic)
        P = np.array([pos[t] for t in nodes]); c0 = np.median(P, 0); r = np.linalg.norm(P - c0, axis=1); rmax = np.percentile(r, 95)
        for t, rr, pp in zip(nodes, r, P):
            if rr > rmax: pos[t] = c0 + (pp - c0) * (rmax + 0.15 * (rr - rmax)) / rr
    bcol = {b: plt.cm.tab20(k % 20) for k, b in enumerate(branches)}
    fs = a.font  # base font size; None keeps the original screen layout
    nr, nc = (3, 2) if a.print_layout else (2, 3)
    W, H = (7.0, 8.6) if a.print_layout else (21, 14)
    fig, axes = plt.subplots(nr, nc, figsize=(W, H)); axes = axes.ravel()
    nscale = 0.25 if a.print_layout else 1.0
    for ax, pop in zip(axes, POPS):
        z = np.array([float(Z[pop][t]["z"]) if t in Z[pop] else 0.0 for t in nodes])
        q = np.array([float(Z[pop][t]["q_bh"]) if t in Z[pop] else 1.0 for t in nodes])
        nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.05 if a.print_layout else 0.12,
                               width=[(0.3 + 1.5 * G[u][v]["w"]) if a.print_layout else (0.5 + 4 * G[u][v]["w"]) for u, v in G.edges()])
        for shp in ("o", "s", "D"):
            idx = [i for i, t in enumerate(nodes) if shape(t) == shp]
            if not idx: continue
            nx.draw_networkx_nodes(G, pos, nodelist=[nodes[i] for i in idx], ax=ax, node_shape=shp,
                                   node_size=[nscale * (20 + 6 * math.sqrt(freq[nodes[i]])) for i in idx],
                                   node_color=[np.clip(z[i], -6, 6) for i in idx], cmap=plt.cm.RdBu_r, vmin=-6, vmax=6,
                                   edgecolors=["black" if q[i] < 0.05 else "#bbbbbb" for i in idx],
                                   linewidths=[(0.9 if a.print_layout else 1.4) if q[i] < 0.05 else 0.3 for i in idx])
        sel = [k for k, t in enumerate(nodes) if abs(z[k]) > a.label_z and q[k] < 0.05]
        if a.max_labels:
            span = np.ptp(np.array(list(pos.values())), axis=0); keep = []
            for k in sorted(sel, key=lambda k: -abs(z[k])):
                if all(abs(pos[nodes[k]][1] - pos[nodes[j]][1]) > 0.06 * span[1] or abs(pos[nodes[k]][0] - pos[nodes[j]][0]) > 0.35 * span[0] for j in keep):
                    keep.append(k)
                if len(keep) >= a.max_labels: break
            sel = keep
        lab = {nodes[k]: short(terms[nodes[k]]["name"])[:28] for k in sel}
        if a.print_layout:
            for t, txt in lab.items():
                ax.annotate(txt, pos[t], fontsize=fs - 1.5, ha="center", va="bottom", xytext=(0, 3), textcoords="offset points",
                            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7), zorder=5)
        else:
            nx.draw_networkx_labels(G, pos, labels=lab, ax=ax, font_size=6.5)
        ax.set_title(f"{pop} vs EUR ({sum(q < 0.05)} of {len(nodes)} terms q<0.05)" if a.print_layout else
                     f"{pop} vs EUR (total; {sum(q < 0.05)} of {len(nodes)} terms q<0.05)", fontsize=fs + 1 if a.print_layout else 12); ax.axis("off")
    ax = axes[5]; ax.axis("off")
    sm = plt.cm.ScalarMappable(cmap=plt.cm.RdBu_r, norm=plt.Normalize(-6, 6)); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.5, shrink=0.6, orientation="horizontal" if a.print_layout else "vertical",
                      label="term z (population minus EUR)" if a.print_layout else "term z (population minus EUR, within gene and source)")
    if a.print_layout: cb.ax.tick_params(labelsize=fs); cb.set_label(cb.ax.get_xlabel(), fontsize=fs)
    for shp, txt in (("o", "mostly bedside/history findings"), ("s", "mostly test-based findings"), ("D", "mixed or unclassified")):
        ax.scatter([], [], marker=shp, c="grey", label=txt)
    ax.scatter([], [], marker="o", facecolors="none", edgecolors="black", linewidths=1.4, label="BH q < 0.05")
    if a.print_layout: ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.6), frameon=False, fontsize=fs)
    else: ax.legend(loc="lower center", frameon=False, fontsize=10)
    ax.text(0.5, 1.0 if a.print_layout else 0.95, f"{len(nodes)} HPO terms 2-3 levels below the root,\nedges: co-occurrence phi > {a.phi} in {N:,} patients\n"
            "size: number of patients", ha="center", va="top", fontsize=fs if a.print_layout else 10, transform=ax.transAxes)
    fig.tight_layout(); fig.savefig(a.out + ".pdf"); fig.savefig(a.out + ".png", dpi=160)
    with open(a.out + "_nodes.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["hpo", "name", "branch", "n_cases", "class"] + [f"z_{p}" for p in POPS] + [f"q_{p}" for p in POPS])
        for t in nodes:
            w.writerow([t, terms[t]["name"], terms[G.nodes[t]["branch"]]["name"] if G.nodes[t]["branch"] in terms else "", freq[t], shape(t)] +
                       [Z[p].get(t, {}).get("z", "") for p in POPS] + [Z[p].get(t, {}).get("q_bh", "") for p in POPS])
    print(f"nodes {len(nodes)}, edges {G.number_of_edges()}, patients {N}")


def cmd_dag(a):
    terms = parse_obo(a.obo); anc = ancestors(terms)
    for pop in POPS:
        Z = load_z(a.terms_dir, pop)
        sig = {t for t, r in Z.items() if float(r["q_bh"]) < a.q and t in terms}
        sig.add(ROOT)
        # transitive reduction: edge t -> nearest significant ancestors
        edges = set()
        for t in sig - {ROOT}:
            up = (anc(t) & sig) - {t}
            for u in up:
                if not any(u in anc(v) and v != u for v in up): edges.add((u, t))
        lines = ["digraph G {", 'graph [overlap=false, splines=true, root="HP:0000118", ranksep=1.3, bgcolor="white"];',
                 'node [style=filled, fontname="Helvetica", fontsize=9, penwidth=0.6];', 'edge [color="#999999", arrowsize=0.4];']
        for t in sig:
            if t == ROOT:
                lines.append(f'"{t}" [label="Phenotypic\\nabnormality", shape=box, fillcolor="#eeeeee"];'); continue
            z = float(Z[t]["z"]); a_ = min(abs(z) / 6, 1)
            col = "#%02x%02x%02x" % ((255, int(255 * (1 - a_)), int(255 * (1 - a_))) if z > 0 else (int(255 * (1 - a_)), int(255 * (1 - a_)), 255))
            w = 0.25 + 0.06 * math.sqrt(int(Z[t]["n_cases"]))
            lab = short(terms[t]["name"]).replace('"', "'")
            lab = lab if len(lab) < 26 else lab[:24] + "..."
            lines.append(f'"{t}" [label="{lab}", shape=ellipse, fillcolor="{col}", width={w:.2f}, fixedsize=false];')
        for u, v in edges: lines.append(f'"{u}" -> "{v}";')
        lines.append(f'labelloc="t"; label="{pop} vs EUR: {len(sig) - 1} HPO terms with BH q < {a.q} (red: more frequent in {pop}, blue: less)"; fontsize=18;')
        lines.append("}")
        dotf = f"{a.out}_{pop}.dot"; open(dotf, "w").write("\n".join(lines))
        for fmt in ("pdf", "png"):
            subprocess.run(["twopi", f"-T{fmt}", dotf, "-o", f"{a.out}_{pop}.{fmt}"], check=True)
        print(pop, len(sig) - 1, "terms", len(edges), "edges")


def cmd_spend(a):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R = list(csv.DictReader(open(a.points)))
    col = {"EUR": "#4c72b0", "MENA": "#dd8452", "EASIA": "#55a868", "SASIA": "#c44e52", "LATAM": "#8172b3", "AFR": "#937860"}
    fig, ax = plt.subplots(figsize=(8, 6))
    for r in R:
        n = int(r["n"]); x = 10 ** float(r["log10_health_exp_pc"]); y = 100 * float(r["mean_resid_test_share"])
        ax.scatter(x, y, s=8 + 2.5 * math.sqrt(n), color=col.get(r["region_mode"], "grey"), alpha=0.75, edgecolors="white", linewidths=0.5)
        if n >= 300: ax.annotate(r["country"], (x, y), fontsize=7, xytext=(3, 2), textcoords="offset points")
    ax.set_xscale("log"); ax.axhline(0, color="grey", lw=0.5)
    ax.set_xlabel("Health expenditure per capita in the country of care (current US$, publication year)")
    ax.set_ylabel("Share of test-based findings, within-gene residual (points)")
    for k, c in col.items(): ax.scatter([], [], color=c, label=k)
    ax.legend(frameon=False, fontsize=8, title="most cases from", title_fontsize=8)
    fig.tight_layout(); fig.savefig(a.out + ".pdf"); fig.savefig(a.out + ".png", dpi=180)


def cmd_popnet(a):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    pops = ["EUR", "MENA", "SASIA", "EASIA", "LATAM", "AFR"]
    n = collections.Counter()
    for f in a.cases:
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r["region"] in pops and r.get("terms"): n[r["region"]] += 1
    edges = [r for r in csv.DictReader(open(a.contrasts), delimiter="\t")]
    within = {(r["population"].split("-")[0], r["analysis"]): (float(r["estimate"]), float(r["q"]))
              for r in csv.DictReader(open(a.within), delimiter="\t") if r["family"] == "within_vs_outside"}
    ang = {p: math.pi / 2 - 2 * math.pi * i / len(pops) for i, p in enumerate(pops)}
    pos = {p: (math.cos(t), math.sin(t)) for p, t in ang.items()}
    cmap = plt.get_cmap("RdBu_r"); vmax = 0.4
    pl = a.print_layout; fs = a.font
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.9) if pl else (12, 6))
    for ax, ver in zip(axes, ("total", "patient-carried")):
        for r in edges:
            if r["version"] != ver: continue
            q = float(r["q_BH"]); e = float(r["excess"])
            if q >= 0.1: continue
            A, rest = r["contrast"].split("_vs_")
            src = next(p for p in pops if A == p + p)
            dst = next(p for p in pops if p != src and rest in (src + p, p + src))
            ax.add_patch(FancyArrowPatch(pos[src], pos[dst], connectionstyle="arc3,rad=0.15", arrowstyle="-|>",
                                         mutation_scale=9 if pl else 14, lw=(0.6 + 3.5 * min(abs(e), 0.8)) if pl else (1 + 6 * min(abs(e), 0.8)), shrinkA=11 if pl else 28, shrinkB=11 if pl else 28,
                                         color="#b2182b" if e > 0 else "#2166ac", alpha=0.8,
                                         linestyle="-" if q < 0.05 else "--"))
        for p in pops:
            e, q = within.get((p, ver), (0.0, 1.0))
            ax.scatter(*pos[p], s=(70 + 600 * math.sqrt(n[p] / max(n.values()))) if pl else (200 + 1800 * math.sqrt(n[p] / max(n.values()))), color=cmap(0.5 + e / (2 * vmax)),
                       edgecolors="black" if q < 0.05 else "grey", linewidths=(1.8 if pl else 2.5) if q < 0.05 else 0.8, zorder=3)
            if pl:
                x0, y0 = pos[p]; rr = 1.33 + 0.4 * math.sqrt(n[p] / max(n.values()))
                ax.annotate(f"{p}\n{n[p]:,}", (x0 * rr, y0 * rr), ha="center", va="center", fontsize=fs - 0.5, zorder=4)
            else:
                ax.annotate(f"{p}\n{n[p]:,}", pos[p], ha="center", va="center", fontsize=9, zorder=4,
                            color="white" if abs(e) > 0.25 else "black")
        ax.set_title({"total": "Total (corpus fixed)", "patient-carried": "Patient-carried (author region fixed)"}[ver], fontsize=fs + 0.5 if pl else None)
        lim = 1.9 if pl else 1.45; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.axis("off")
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color="#b2182b", lw=3, label="A-A more similar than A-B"),
                            Line2D([], [], color="#2166ac", lw=3, label="A-A less similar than A-B"),
                            Line2D([], [], color="grey", lw=2, label="q < 0.05"),
                            Line2D([], [], color="grey", lw=2, ls="--", label="0.05 <= q < 0.1")],
                   loc="lower center", ncol=2 if pl else 4, frameon=False, fontsize=fs if pl else 9, title="arrow A -> B", title_fontsize=fs if pl else 9,
               bbox_to_anchor=(0.45, 0.0 if pl else 0.08))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-vmax, vmax))
    cb = fig.colorbar(sm, ax=axes, shrink=0.6, label="within-population minus\noutside similarity (z)" if pl else "within-population minus outside similarity (z)")
    if pl: cb.ax.tick_params(labelsize=fs); cb.set_label(cb.ax.get_ylabel(), fontsize=fs)
    fig.savefig(a.out + ".pdf", bbox_inches="tight"); fig.savefig(a.out + ".png", dpi=180, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["shift", "dag", "spend", "popnet"])
    ap.add_argument("--obo"); ap.add_argument("--cases", action="append"); ap.add_argument("--terms-dir")
    ap.add_argument("--termclass"); ap.add_argument("--points"); ap.add_argument("--contrasts"); ap.add_argument("--within"); ap.add_argument("--out", required=True)
    ap.add_argument("--n-nodes", type=int, default=140); ap.add_argument("--phi", type=float, default=0.08)
    ap.add_argument("--label-z", type=float, default=3.5); ap.add_argument("--q", type=float, default=0.05)
    ap.add_argument("--print-layout", action="store_true", help="shift/popnet: 178 mm wide layout with print-size fonts")
    ap.add_argument("--font", type=float, default=8, help="base font size for --print-layout")
    ap.add_argument("--max-labels", type=int, default=0, help="shift: label at most N terms per panel (largest |z|)")
    a = ap.parse_args()
    {"shift": cmd_shift, "dag": cmd_dag, "spend": cmd_spend, "popnet": cmd_popnet}[a.cmd](a)


if __name__ == "__main__":
    main()
