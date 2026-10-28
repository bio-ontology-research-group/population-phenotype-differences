#!/usr/bin/env python3
"""Figures for the Genetics in Medicine manuscript (population differences in recorded phenotypes).

Writes PDFs (and PNG previews) to results/v4/figures/:
  gim_fig1_testshare   share of test-based findings vs EUR across analysis versions, and
                       log(1+n) test-based and bedside counts (terms_all/all_outcomes.tsv, qvalues.tsv)
  gim_fig2_map         world map of the within-stratum residual test share by country of care
                       (worldbank/healthspend_country_points.csv; Natural Earth 110m shapes, --shapes)
  gim_fig2_setting     (a) ancestry vs place of care (setting_all/*_outcomes.csv, test_share_q.tsv);
                       (b) country-level residual test share vs health expenditure (worldbank/)
  gim_fig3_shift      (a) MH z-scores of all tested HPO classes by finding class and population
                       (terms_all/main_*_total_terms.csv, termclass.jsonl); (b) dot matrix of the grouping classes with
                       the largest shifts (figures/fig_shift_network_nodes.csv from network_figs.py shift)
  gim_figS_decomp      pair-similarity decomposition for MENA-EUR and EUR-EASIA (allpairs/*_decomposition.csv)
  gim_figS_organ       organ-system differences in log(1+n) counts (terms_all/main_*_total_outcomes.csv, qvalues.tsv)
  gim_figS_recovery    causative-gene recovery, population minus EUR (recovery/recovery_*_recovery.csv, qvalues.tsv)

Usage: python3 gim_figures.py [--res results/v4]
Points are the observed Mantel-Haenszel difference minus the permutation mean ("excess");
bars are the central 95% of the permutation distribution, shifted by the same mean.
"""
import argparse
import csv
import math
import os

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

MM = 1 / 25.4
POPS = ["MENA", "SASIA", "LATAM", "AFR", "EASIA"]
POPNAME = {"MENA": "MENA", "SASIA": "South Asia", "LATAM": "Latin America", "AFR": "Sub-Saharan Africa",
           "EASIA": "East Asia", "EUR": "EUR"}
# Okabe-Ito
OI = {"black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
      "blue": "#0072B2", "red": "#D55E00", "purple": "#CC79A7", "grey": "#999999"}
POPCOL = {"EUR": OI["blue"], "MENA": OI["orange"], "EASIA": OI["green"], "SASIA": OI["red"],
          "LATAM": OI["purple"], "AFR": OI["black"]}

plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "font.family": "DejaVu Sans",
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})


def read_tsv(path, delim="\t"):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def panel_label(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right", va="bottom")


def qlookup(res):
    q = {}
    for r in read_tsv(os.path.join(res, "qvalues.tsv")):
        q[(r["family"], r["analysis"], r["outcome"])] = float(r["q"])
    return q


def save(fig, res, name):
    out = os.path.join(res, "figures", name)
    fig.savefig(out + ".pdf", bbox_inches="tight")
    fig.savefig(out + ".png", dpi=200, bbox_inches="tight")
    print("wrote", out + ".pdf")


def excess_row(r, scale=1.0):
    nm = float(r["null_mean"])
    return ((float(r["diff"]) - nm) * scale, (float(r["null_lo"]) - nm) * scale, (float(r["null_hi"]) - nm) * scale)


def fig1(res, q):
    rows = read_tsv(os.path.join(res, "terms_all", "all_outcomes.tsv"))
    by = {(r["run"], r["outcome"]): r for r in rows}
    versions = [("total", "main_{p}_total", "outcomes:total", OI["black"], "o"),
                ("patient-carried", "main_{p}_pc", "outcomes:pc", OI["sky"], "s"),
                ("age-adjusted", "adjage_{p}", "outcomes:adjage", OI["orange"], "D"),
                ("sex-adjusted", "adjsex_{p}", "outcomes:adjsex", OI["purple"], "^")]
    share = "share of findings needing a test"
    fig, (a, b) = plt.subplots(1, 2, figsize=(178 * MM, 72 * MM), gridspec_kw={"width_ratios": [1.25, 1]})
    off = [0.27, 0.09, -0.09, -0.27]
    for i, p in enumerate(POPS):
        y0 = len(POPS) - 1 - i
        for (lab, run, fam, col, mk), dy in zip(versions, off):
            r = by[(run.format(p=p), share)]
            e, lo, hi = excess_row(r, 100)
            qq = q.get((fam, run.format(p=p), share), 1.0)
            a.plot([lo, hi], [y0 + dy] * 2, color="#d0d0d0", lw=2.2, solid_capstyle="butt", zorder=1)
            a.plot(e, y0 + dy, marker=mk, ms=4.2, color=col, mfc=col if qq < 0.05 else "white", mew=0.9, zorder=3)
    a.axvline(0, color="grey", lw=0.6, zorder=0)
    a.set_yticks(range(len(POPS)))
    a.set_yticklabels([POPNAME[p] for p in POPS[::-1]])
    a.set_xlabel("Share of test-based findings, population minus EUR (percentage points)")
    a.set_ylim(-0.6, len(POPS) - 0.4)
    a.legend(handles=[Line2D([], [], color=c, marker=m, ls="", ms=4.2, label=l) for l, _, _, c, m in versions],
             loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=4, frameon=False, handletextpad=0.2,
             columnspacing=0.8)
    panel_label(a, "a")

    outs = [("log(1+ test-based findings)", "test-based", OI["blue"], "s"),
            ("log(1+ bedside/history findings)", "bedside", OI["red"], "o")]
    for i, p in enumerate(POPS):
        y0 = len(POPS) - 1 - i
        for (o, lab, col, mk), dy in zip(outs, (0.13, -0.13)):
            run = f"main_{p}_total"
            e, lo, hi = excess_row(by[(run, o)])
            qq = q.get(("outcomes:total", run, o), 1.0)
            b.plot([lo, hi], [y0 + dy] * 2, color="#d0d0d0", lw=2.2, solid_capstyle="butt", zorder=1)
            b.plot(e, y0 + dy, marker=mk, ms=4.2, color=col, mfc=col if qq < 0.05 else "white", mew=0.9, zorder=3)
    b.axvline(0, color="grey", lw=0.6, zorder=0)
    b.set_yticks(range(len(POPS)))
    b.set_yticklabels([])
    b.set_ylim(-0.6, len(POPS) - 0.4)
    b.set_xlabel("log(1 + number of findings), population minus EUR")
    b.legend(handles=[Line2D([], [], color=c, marker=m, ls="", ms=4.2, label=l) for _, l, c, m in outs],
             loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False, handletextpad=0.2)
    panel_label(b, "b")
    fig.tight_layout(w_pad=1.0)
    save(fig, res, "gim_fig1_testshare")


def fig2(res):
    sdir = os.path.join(res, "setting_all")
    qs = {r["contrast"]: float(r["q_BH"]) for r in read_tsv(os.path.join(sdir, "test_share_q.tsv"))}
    share = "share of findings needing a test"
    fig, (a, b) = plt.subplots(1, 2, figsize=(178 * MM, 78 * MM), gridspec_kw={"width_ratios": [1, 1.15]})
    contrasts = [("{p}anc_{p}care_vs_EUR", "own-region care vs reference", OI["red"], "o"),
                 ("{p}anc_EURcare_vs_EUR", "EUR/N. America care vs reference", OI["blue"], "s"),
                 ("{p}anc_{p}care_vs_{p}anc_EURcare", "own-region vs EUR/N. America care", OI["black"], "D")]
    order = ["MENA", "SASIA", "AFR", "LATAM", "EASIA"]
    for i, p in enumerate(order):
        y0 = len(order) - 1 - i
        for (c, lab, col, mk), dy in zip(contrasts, (0.22, 0.0, -0.22)):
            name = c.format(p=p)
            r = next(x for x in read_tsv(os.path.join(sdir, f"setting_{name}_outcomes.csv"), ",")
                     if x["outcome"] == share)
            e, lo, hi = excess_row(r, 100)
            qq = qs[name]
            a.plot([lo, hi], [y0 + dy] * 2, color="#d0d0d0", lw=2.2, solid_capstyle="butt", zorder=1)
            a.plot(e, y0 + dy, marker=mk, ms=4.2, color=col, mfc=col if qq < 0.05 else "white", mew=0.9, zorder=3)
    a.axvline(0, color="grey", lw=0.6, zorder=0)
    a.set_yticks(range(len(order)))
    a.set_yticklabels([POPNAME[p] + " ancestry" for p in order[::-1]])
    a.set_ylim(-0.6, len(order) - 0.4)
    a.set_xlabel("Share of test-based findings, difference (percentage points)")
    a.legend(handles=[Line2D([], [], color=c, marker=m, ls="", ms=4.2, label=l) for _, l, c, m in contrasts],
             loc="upper center", bbox_to_anchor=(0.45, -0.17), frameon=False, ncol=1, handletextpad=0.2)
    panel_label(a, "b")

    pts = read_tsv(os.path.join(res, "worldbank", "healthspend_country_points.csv"), ",")
    slope = next(r for r in read_tsv(os.path.join(res, "worldbank", "healthspend_slopes.csv"), ",")
                 if r["covariate"] == "health_exp_pc" and r["outcome"] == share)
    gulf = {"SAU", "ARE", "QAT", "KWT", "BHR"}
    labels = {"SAU": (6, -14), "USA": (-14, 10), "JPN": (-30, 9), "CHN": (-4, 12), "IND": (-8, -14),
              "TUR": (-28, -12)}
    W = sum(int(r["n"]) for r in pts)
    xm = sum(int(r["n"]) * float(r["log10_health_exp_pc"]) for r in pts) / W
    ym = sum(int(r["n"]) * 100 * float(r["mean_resid_test_share"]) for r in pts) / W
    for r in sorted(pts, key=lambda r: -int(r["n"])):
        n = int(r["n"])
        x = float(r["log10_health_exp_pc"])
        y = 100 * float(r["mean_resid_test_share"])
        col = POPCOL.get(r["region_mode"], OI["grey"])
        isg = r["iso"] in gulf
        b.scatter(10 ** x, y, s=3 + 0.9 * math.sqrt(n), color=col, alpha=0.8,
                  edgecolors="black" if isg else "white", linewidths=0.9 if isg else 0.4, zorder=3 if isg else 2)
        if r["iso"] in labels:
            dx, dy = labels[r["iso"]]
            b.annotate(r["country"].replace("United States", "USA"), (10 ** x, y), fontsize=5.5,
                       xytext=(dx, dy), textcoords="offset points", zorder=4,
                       arrowprops=dict(arrowstyle="-", lw=0.4, color="grey"))
    s = 100 * float(slope["slope_per_10fold"])
    xs = [1.0, 4.1]
    b.plot([10 ** v for v in xs], [ym + s * (v - xm) for v in xs], color="grey", lw=0.8, ls="--", zorder=1)
    b.axhline(0, color="#d0d0d0", lw=0.6, zorder=0)
    b.set_xscale("log")
    b.set_xlim(10 ** 1.0, 10 ** 4.15)
    b.set_ylim(-25, 25)
    b.set_xlabel("Health expenditure per capita, country of care (US$)")
    b.set_ylabel("Residual share of test-based findings (points)")
    b.legend(handles=[Line2D([], [], color=POPCOL[k], marker="o", ls="", ms=4, label=POPNAME[k])
                      for k in ["EUR", "MENA", "EASIA", "SASIA", "LATAM", "AFR"]]
             + [Line2D([], [], color="white", marker="o", mec="black", ls="", ms=4, label="Gulf state")],
             loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False, handletextpad=0.1,
             columnspacing=0.8, title="Most frequent population of the country's cases", title_fontsize=6.5)
    panel_label(b, "c")
    fig.tight_layout(w_pad=4.0)
    save(fig, res, "gim_fig2_setting")


def fig3(res, termclass):
    import json
    import numpy as np
    from matplotlib.colors import TwoSlopeNorm
    tc = {x["id"]: x["class"] for x in map(json.loads, open(termclass))}
    test = {"laboratory", "imaging", "physiological_test", "pathology"}
    bed = {"clinical_exam", "symptom_history"}
    fig = plt.figure(figsize=(178 * MM, 120 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.0], wspace=1.05, left=0.08, right=0.93, bottom=0.2, top=0.93)
    a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1])
    # (a) all tested classes
    for i, P in enumerate(POPS):
        z = {"b": [], "t": []}
        for r in read_tsv(os.path.join(res, "terms_all", f"main_{P}_total_terms.csv"), ","):
            c = tc.get(r["hpo"])
            if c in bed: z["b"].append(float(r["z"]))
            elif c in test: z["t"].append(float(r["z"]))
        for k, dx, col in (("b", -0.18, OI["red"]), ("t", 0.18, OI["blue"])):
            v = np.clip(z[k], -6, 6)
            parts = a.violinplot([v], positions=[i + dx], widths=0.34, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor(col); pc.set_edgecolor("none"); pc.set_alpha(0.45)
            a.plot([i + dx - 0.12, i + dx + 0.12], [np.median(z[k])] * 2, color=col, lw=1.6)
    a.axhline(0, color=OI["grey"], lw=0.6)
    a.set_xticks(range(len(POPS))); a.set_xticklabels([POPNAME[p] for p in POPS], rotation=35, ha="right")
    a.set_ylabel("HPO class frequency vs EUR (Mantel-Haenszel z)"); a.set_ylim(-6.3, 6.3)
    a.legend(handles=[Line2D([], [], color=OI["red"], lw=4, alpha=0.6, label="bedside or history class"),
                      Line2D([], [], color=OI["blue"], lw=4, alpha=0.6, label="test-based class"),
                      Line2D([], [], color="black", lw=1.6, label="median")], frameon=False, loc="upper center",
             bbox_to_anchor=(0.5, -0.2), ncol=1)
    panel_label(a, "b", x=-0.25, y=1.05)
    # (b) dot matrix of grouping classes with the largest shifts
    def fz(r, P):
        return float(r[f"z_{P}"]) if r[f"z_{P}"] else float("nan")

    def fq(r, P):
        return float(r[f"q_{P}"]) if r[f"q_{P}"] else 1.0
    rows = read_tsv(os.path.join(res, "figures", "fig_shift_network_nodes.csv"), ",")
    sig = [r for r in rows if any(fq(r, P) < 0.05 for P in POPS)]
    sig.sort(key=lambda r: -np.nanmax([abs(fz(r, P)) for P in POPS]))
    sel = sig[:24]
    order = {"o": 0, "s": 1, "D": 2}
    sel.sort(key=lambda r: (order[r["class"]], -np.nanmean([fz(r, P) for P in POPS[:4]])))
    norm = TwoSlopeNorm(vmin=-8, vcenter=0, vmax=8); cmap = plt.get_cmap("RdBu_r")
    for y, r in enumerate(sel):
        for x, P in enumerate(POPS):
            zv = fz(r, P); qv = fq(r, P)
            if zv != zv:
                b.scatter(x, y, s=6, marker="x", color=OI["grey"], lw=0.6); continue
            b.scatter(x, y, s=12 + 9 * min(abs(zv), 8), color=cmap(norm(max(-8, min(8, zv)))),
                      edgecolors="black" if qv < 0.05 else "none", linewidths=0.8, zorder=3)
    b.set_yticks(range(len(sel)))
    def lab(n):
        n = n.replace("Abnormality of the ", "").replace("Abnormality of ", "").replace("Abnormal ", "")[:42]
        return n[0].upper() + n[1:]
    b.set_yticklabels([lab(r["name"]) for r in sel])
    for t, r in zip(b.get_yticklabels(), sel):
        t.set_color({"o": OI["red"], "s": OI["blue"], "D": "#555555"}[r["class"]])
    b.set_xticks(range(len(POPS))); b.set_xticklabels([POPNAME[p] for p in POPS], rotation=35, ha="right")
    b.set_xlim(-0.6, len(POPS) - 0.4); b.set_ylim(len(sel) - 0.5, -0.5)
    for y in range(1, len(sel)):
        if sel[y]["class"] != sel[y - 1]["class"]: b.axhline(y - 0.5, color=OI["grey"], lw=0.5)
    b.spines["left"].set_visible(False); b.tick_params(axis="y", length=0)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=b, shrink=0.45, pad=0.04); cb.set_label("z vs EUR"); cb.outline.set_linewidth(0.4)
    panel_label(b, "c", x=-1.25, y=1.05)
    b.scatter([], [], s=12 + 9 * 4, color="white", edgecolors="black", linewidths=0.8, label="q < 0.05")
    b.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2))
    save(fig, res, "gim_fig3_shift")


def fig2_map(res, shapes, min_n=20):
    """Choropleth of the mean within-stratum residual share of test-based findings by country of care."""
    import shapefile
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import TwoSlopeNorm
    pts = {r["iso"]: r for r in read_tsv(os.path.join(res, "worldbank", "healthspend_country_points.csv"), ",")}
    rd = shapefile.Reader(shapes)
    fields = [f[0] for f in rd.fields if f[0] != "DeletionFlag"]
    norm = TwoSlopeNorm(vmin=-8, vcenter=0, vmax=8); cmap = plt.get_cmap("RdBu_r")
    fig, ax = plt.subplots(figsize=(178 * MM, 88 * MM))
    polys, cols, shown = [], [], 0
    for sr in rd.iterShapeRecords():
        rec = dict(zip(fields, sr.record))
        iso = rec.get("ISO_A3") if rec.get("ISO_A3") not in (None, "-99") else rec.get("ADM0_A3")
        if iso == "ATA": continue
        r = pts.get(iso)
        col = "#e6e6e6"
        if r and int(r["n"]) >= min_n:
            col = cmap(norm(max(-8, min(8, 100 * float(r["mean_resid_test_share"]))))); shown += 1
        elif r:
            col = "#c8c8c8"
        pts_ = sr.shape.points; parts = list(sr.shape.parts) + [len(pts_)]
        for k in range(len(parts) - 1):
            polys.append(pts_[parts[k]:parts[k + 1]]); cols.append(col)
    ax.add_collection(PolyCollection(polys, facecolors=cols, edgecolors="white", linewidths=0.25))
    ax.set_xlim(-170, 180); ax.set_ylim(-58, 84); ax.set_aspect("equal"); ax.axis("off")
    panel_label(ax, "a", x=0.0, y=0.98)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=ax, orientation="horizontal", shrink=0.45, pad=0.02, aspect=35)
    cb.set_label("Share of test-based findings, within-gene residual of the country of care (points)")
    cb.outline.set_linewidth(0.4)
    fig.tight_layout()
    save(fig, res, "gim_fig2_map")
    return shown


def figS_decomp(res):
    """Selected coefficients of the pair-similarity model for all 15 population pairs (patient-carried model)."""
    import glob, re
    terms = [("curation[store|store]", "Store/Store", (-1, 6)), ("curation[pmcp|store]", "PMC-Patients/Store", (-1, 4)),
             ("same_group", "same group", (-1, 3)), ("author_region_same", "same author region", (-1.5, 2)),
             ("population_same", "same population", (-1.5, 1.5))]
    files = sorted(glob.glob(os.path.join(res, "allpairs", "popmodel_fine_*_decomposition.csv")))
    pairs = []
    for f in files:
        m = re.search(r"popmodel_fine_([A-Z]+)_([A-Z]+)_decomposition", f)
        pairs.append((m.group(1), m.group(2), {r["term"]: r for r in read_tsv(f, ",")}))
    order = ["EUR", "MENA", "SASIA", "LATAM", "AFR", "EASIA"]
    pairs.sort(key=lambda x: (order.index(x[0]), order.index(x[1])))
    fig, axes = plt.subplots(1, len(terms), figsize=(178 * MM, 95 * MM), sharey=True)
    for ax, (t, name, xl) in zip(axes, terms):
        col = OI["blue"] if t.startswith("curation") else (OI["red"] if "group" in t or "region" in t else OI["green"])
        for y, (A, B, d) in enumerate(pairs):
            r = d.get(t)
            if r is None: continue
            c, lo, hi = float(r["coef"]), float(r["lo"]), float(r["hi"])
            ax.plot([max(lo, xl[0]), min(hi, xl[1])], [y, y], color=col, lw=1.1, alpha=0.6); ax.plot(c, y, "o", color=col, ms=3.4)
            if hi > xl[1]: ax.plot(xl[1], y, ">", color=col, ms=3)
            if lo < xl[0]: ax.plot(xl[0], y, "<", color=col, ms=3)
        ax.set_xlim(*xl); ax.axvline(0, color="grey", lw=0.6); ax.set_title(name, fontsize=7); ax.set_xlabel("z-similarity")
    axes[0].set_yticks(range(len(pairs)))
    axes[0].set_yticklabels([f"{POPNAME[A]} / {POPNAME[B]}" for A, B, _ in pairs]); axes[0].invert_yaxis()
    fig.tight_layout(w_pad=0.8)
    save(fig, res, "gim_figS_decomp")


def figS_organ(res, q):
    systems = None
    data = {}
    for p in POPS:
        rows = read_tsv(os.path.join(res, "terms_all", f"main_{p}_total_outcomes.csv"), ",")
        for r in rows:
            if r["outcome"].startswith("organ["):
                sysname, kind = r["outcome"][6:].split("] ")
                e = (float(r["diff"]) - float(r["null_mean"])) * 100
                qq = q.get(("organ:total", f"main_{p}_total", r["outcome"]), 1.0)
                data[(sysname, kind, p)] = (e, qq)
    systems = sorted({k[0] for k in data}, key=str.lower)
    kinds = [("log(1+ findings)", "all findings"), ("log(1+ test-based)", "test-based"), ("log(1+ bedside)", "bedside")]
    fig, axes = plt.subplots(1, 3, figsize=(178 * MM, 110 * MM), sharey=True)
    vmax = 12
    cmap = plt.get_cmap("RdBu_r")
    for ax, (k, title), lab in zip(axes, kinds, "abc"):
        for i, s in enumerate(systems):
            for j, p in enumerate(POPS):
                v = data.get((s, k, p))
                if v is None:
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, color="#f2f2f2"))
                    continue
                e, qq = v
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, color=cmap(0.5 + max(-1, min(1, e / vmax)) / 2)))
                if qq < 0.05:
                    ax.plot(j, i, "o", color="black", ms=2.2)
        ax.set_xlim(-0.5, len(POPS) - 0.5)
        ax.set_ylim(len(systems) - 0.5, -0.5)
        ax.set_xticks(range(len(POPS)))
        ax.set_xticklabels(POPS, rotation=45, ha="right")
        ax.set_title(title)
        for sp in ax.spines.values():
            sp.set_visible(False)
        panel_label(ax, lab)
    axes[0].set_yticks(range(len(systems)))
    axes[0].set_yticklabels(systems)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-vmax, vmax))
    fig.colorbar(sm, ax=axes, shrink=0.5, label="Difference from EUR, log(1+n) × 100")
    save(fig, res, "gim_figS_organ")


def figS_recovery(res, q):
    fig, ax = plt.subplots(figsize=(86 * MM, 62 * MM))
    vers = [("total", OI["black"], "o", 0.12), ("patient-carried", OI["sky"], "s", -0.12)]
    for i, p in enumerate(POPS):
        y0 = len(POPS) - 1 - i
        rows = read_tsv(os.path.join(res, "recovery", f"recovery_{p}_recovery.csv"), ",")
        for v, col, mk, dy in vers:
            r = next(x for x in rows if x["null"] == v)
            nm = float(r["null_mean"])
            e = float(r["stat"]) - nm
            qq = q.get(("recovery", v, "mean log2 rank difference", p), None)
            if qq is None:
                qq = next((float(x["q"]) for x in read_tsv(os.path.join(res, "qvalues.tsv"))
                           if x["family"] == "recovery" and x["population"] == p and x["analysis"] == v), 1.0)
            ax.plot([float(r["lo"]) - nm, float(r["hi"]) - nm], [y0 + dy] * 2, color="#d0d0d0", lw=2.2,
                    solid_capstyle="butt")
            ax.plot(e, y0 + dy, marker=mk, ms=4, color=col, mfc=col if qq < 0.05 else "white", mew=0.9)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_yticks(range(len(POPS)))
    ax.set_yticklabels([POPNAME[p] for p in POPS[::-1]])
    ax.set_xlabel("Causative-gene log2 rank, minus EUR")
    ax.legend(handles=[Line2D([], [], color=c, marker=m, ls="", ms=4, label=v) for v, c, m, _ in vers],
              frameon=False, loc="lower right")
    fig.tight_layout()
    save(fig, res, "gim_figS_recovery")


def figS_agesex(res, q):
    """Deficit in test share by age band and sex (subset analyses, terms_all/age_*/sex_* runs)."""
    share = "share of findings needing a test"
    subsets = [("age_infant", "infant"), ("age_child", "child"), ("age_adolescent", "adolescent"), ("age_adult", "adult"),
               ("sex_male", "male"), ("sex_female", "female")]
    cols = [OI["orange"], OI["red"], OI["purple"], OI["blue"], OI["black"], OI["grey"]]
    mks = ["o", "o", "o", "o", "s", "s"]
    fig, ax = plt.subplots(figsize=(178 * MM, 80 * MM))
    offs = np.linspace(0.32, -0.32, len(subsets))
    for i, p in enumerate(POPS):
        y0 = len(POPS) - 1 - i
        for (run, lab), dy, col, mk in zip(subsets, offs, cols, mks):
            f = os.path.join(res, "terms_all", f"{run}_{p}_outcomes.csv")
            r = next(x for x in read_tsv(f, ",") if x["outcome"] == share)
            e, lo, hi = excess_row(r, 100)
            qq = q.get(("strata", f"{run}_{p}", share), 1.0)
            ax.plot([lo, hi], [y0 + dy] * 2, color="#d0d0d0", lw=1.8, solid_capstyle="butt", zorder=1)
            ax.plot(e, y0 + dy, marker=mk, ms=3.8, color=col, mfc=col if qq < 0.05 else "white", mew=0.9, zorder=3)
    ax.axvline(0, color="grey", lw=0.6, zorder=0)
    ax.set_yticks(range(len(POPS)))
    ax.set_yticklabels([POPNAME[p] for p in POPS[::-1]])
    ax.set_ylim(-0.6, len(POPS) - 0.4)
    ax.set_xlabel("Share of test-based findings, population minus EUR (percentage points)")
    ax.legend(handles=[Line2D([], [], color=c, marker=m, ls="", ms=4, label=l) for (_, l), c, m in zip(subsets, cols, mks)],
              loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=6, frameon=False, handletextpad=0.2)
    fig.tight_layout()
    save(fig, res, "gim_figS_agesex")


def figS_coreport(res, thr=0.10, vmax=0.15):
    """Organ-system co-reporting network (test/bedside split): EUR backbone and population excess."""
    E = read_tsv(os.path.join(res, "comorbidity", "organ_split_edges.csv"), ",")
    N = read_tsv(os.path.join(res, "comorbidity", "organ_split_nodes.csv"), ",")
    G = read_tsv(os.path.join(res, "comorbidity", "organ_split_global.csv"), ",")
    pg = {r["pop"]: float(r["p_global"]) for r in G}
    lab = {int(r["node"]): r["label"] for r in N}
    systems = []
    for i in sorted(lab):
        sname = lab[i].rsplit(" (", 1)[0]
        if sname not in systems: systems.append(sname)
    ang = {}
    k = len(systems)
    for si, sname in enumerate(systems):
        a0 = math.pi / 2 - 2 * math.pi * si / k
        for i in lab:
            if lab[i].rsplit(" (", 1)[0] == sname:
                kind = lab[i].rsplit(" (", 1)[1].rstrip(")")
                ang[i] = a0 + (0.06 if kind == "test" else -0.06)
    pos = {i: (math.cos(a), math.sin(a)) for i, a in ang.items()}
    backbone = [r for r in E if r["pop"] == "MENA" and float(r["r_ref_all"]) > thr]
    keys = {(int(r["i"]), int(r["j"])) for r in backbone}
    fig = plt.figure(figsize=(178 * MM, 215 * MM))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.35, 1, 1])
    ax0 = fig.add_subplot(gs[0, :])
    axes = [fig.add_subplot(gs[r, c]) for r in (1, 2) for c in range(3)]
    def draw_nodes(ax, big):
        for i, (x, y) in pos.items():
            kind = lab[i].rsplit(" (", 1)[1].rstrip(")")
            ax.plot(x, y, marker="s" if kind == "test" else "o", ms=5 if big else 3, color=OI["blue"] if kind == "test" else OI["red"],
                    mec="white", mew=0.4, zorder=3)
        ax.set_xlim(-1.75 if big else -1.15, 1.75 if big else 1.15); ax.set_ylim(-1.45 if big else -1.15, 1.45 if big else 1.15)
        ax.set_aspect("equal"); ax.axis("off")
    for r in backbone:
        i, j = int(r["i"]), int(r["j"]); w = float(r["r_ref_all"])
        ax0.plot(*zip(pos[i], pos[j]), color="#555555", lw=0.4 + 5 * (w - thr), alpha=0.6, zorder=1)
    draw_nodes(ax0, True)
    for si, sname in enumerate(systems):
        a0 = math.pi / 2 - 2 * math.pi * si / k
        rad = 1.1 if abs(math.sin(a0)) < 0.7 else (1.09 if si % 2 == 0 else 1.3)
        rad = {"metabolism/homeostasis": 1.1, "blood and blood-forming tissues": 1.32, "cardiovascular system": 1.12,
               "voice": 1.3}.get(sname, rad)
        x, y = rad * math.cos(a0), rad * math.sin(a0)
        ax0.text(x, y, sname, fontsize=6.5, ha="left" if x > 0.05 else "right" if x < -0.05 else "center",
                 va="bottom" if y > 0.05 else "top" if y < -0.05 else "center")
    ax0.set_title(f"EUR co-reporting backbone (within-gene correlation r > {thr})", fontsize=7.5)
    ax0.plot([], [], "o", color=OI["red"], ms=4, label="bedside findings"); ax0.plot([], [], "s", color=OI["blue"], ms=4, label="test-based findings")
    ax0.legend(loc="upper left", frameon=False, fontsize=6.5)
    panel_label(ax0, "a", x=0.02)
    cmap = plt.get_cmap("RdBu_r")
    for ax, p, l in zip(axes, POPS, "bcdef"):
        ex = {(int(r["i"]), int(r["j"])): float(r["excess"]) for r in E if r["pop"] == p}
        for (i, j) in keys:
            e = ex[(i, j)]
            ax.plot(*zip(pos[i], pos[j]), color=cmap(0.5 + max(-1, min(1, e / vmax)) / 2), lw=1.1, zorder=1)
        draw_nodes(ax, False)
        ax.set_title(f"{POPNAME[p]} (global p = {pg[p]:.3g})", fontsize=7)
        panel_label(ax, l, x=0.05)
    axes[-1].axis("off")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-vmax, vmax))
    cb = fig.colorbar(sm, ax=axes[-1], orientation="horizontal", fraction=0.4, shrink=0.9)
    cb.set_label("excess in Fisher-z correlation,\npopulation minus EUR", fontsize=6.5); cb.ax.tick_params(labelsize=6)
    fig.tight_layout()
    save(fig, res, "gim_figS_coreport")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--res", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "v4"))
    ap.add_argument("--shapes", help="Natural Earth ne_110m_admin_0_countries.shp (for gim_fig2_map; skipped if absent)")
    a = ap.parse_args()
    q = qlookup(a.res)
    fig1(a.res, q)
    fig2(a.res)
    if a.shapes: fig2_map(a.res, a.shapes)
    fig3(a.res, os.path.join(os.path.dirname(a.res.rstrip("/")), "termclass.jsonl"))
    figS_decomp(a.res)
    figS_organ(a.res, q)
    figS_recovery(a.res, q)
    figS_agesex(a.res, q)
    figS_coreport(a.res)


if __name__ == "__main__":
    main()
