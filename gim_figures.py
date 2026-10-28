#!/usr/bin/env python3
"""Figures for the Genetics in Medicine manuscript (population differences in recorded phenotypes).

Writes PDFs (and PNG previews) to results/v4/figures/:
  gim_fig1_testshare   share of test-based findings vs EUR across analysis versions, and
                       log(1+n) test-based and bedside counts (terms_all/all_outcomes.tsv, qvalues.tsv)
  gim_fig2_setting     (a) ancestry vs place of care (setting_all/*_outcomes.csv, test_share_q.tsv);
                       (b) country-level residual test share vs health expenditure (worldbank/)
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
    panel_label(a, "a")

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
    panel_label(b, "b")
    fig.tight_layout(w_pad=4.0)
    save(fig, res, "gim_fig2_setting")


def figS_decomp(res):
    terms = [("curation[store|store]", "Store / Store"), ("curation[pmcp|store]", "PMC-Patients / Store"),
             ("curation[pmcp|pmcp]", "PMC-Patients / PMC-Patients"),
             ("curation[pavs_clin|pubmed_abs]", "Saudi clinical / abstract"),
             ("same_group", "same reporting group"), ("author_region_same", "same author region"),
             ("population_same", "same population"), ("population_n_{A}", "number of {A} cases")]
    fig, axes = plt.subplots(1, 2, figsize=(178 * MM, 62 * MM))
    for ax, (A, B), lab in zip(axes, [("MENA", "EUR"), ("EUR", "EASIA")], "ab"):
        d = {r["term"]: r for r in read_tsv(os.path.join(res, "allpairs", f"popmodel_fine_{A}_{B}_decomposition.csv"), ",")}
        for i, (t, name) in enumerate(terms):
            r = d.get(t.format(A=A))
            y = len(terms) - 1 - i
            if r is None:
                continue
            c, lo, hi = float(r["coef"]), float(r["lo"]), float(r["hi"])
            col = OI["blue"] if t.startswith("curation") else (OI["red"] if "group" in t or "region" in t else OI["green"])
            ax.plot([lo, hi], [y, y], color=col, lw=1.2, alpha=0.6)
            ax.plot(c, y, "o", color=col, ms=3.8)
        ax.axvline(0, color="grey", lw=0.6)
        ax.set_yticks(range(len(terms)))
        ax.set_yticklabels([n.format(A=A) if "{A}" in n else n for _, n in terms][::-1])
        ax.set_xlabel("Coefficient (z-similarity)")
        ax.set_title(f"{A}–{B}")
        panel_label(ax, lab)
    axes[1].set_yticklabels([n.format(A="EUR") if "{A}" in n else n for _, n in terms][::-1])
    axes[1].tick_params(labelleft=True)
    fig.tight_layout(w_pad=1.5)
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


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--res", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "v4"))
    a = ap.parse_args()
    q = qlookup(a.res)
    fig1(a.res, q)
    fig2(a.res)
    figS_decomp(a.res)
    figS_organ(a.res, q)
    figS_recovery(a.res, q)


if __name__ == "__main__":
    main()
