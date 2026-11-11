#!/usr/bin/env python3
"""Study overview schematic (main Figure 1).

Draws the flow from case sources through phenotype extraction and finding classification to the three
analyses. Counts are those of Table 1 of the manuscript. Colours follow gim_figures.py: bedside red,
test-based blue (Okabe-Ito).

  python3 overview_figure.py --out results/v5/figures/fig_overview
"""
import argparse, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

MM = 1 / 25.4
BED, TEST, MIX, INK, MUTE, LINE, FILL = "#D55E00", "#0072B2", "#999999", "#222222", "#555555", "#888888", "#F4F4F4"
plt.rcParams.update({"font.size": 7, "font.family": "DejaVu Sans", "pdf.fonttype": 42})


def box(ax, x, y, w, h, title=None, fc=FILL, ec=LINE, lw=0.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=fc, ec=ec, lw=lw))
    if title: ax.text(x + 2, y + h - 2.2, title, fontsize=7, fontweight="bold", va="top", color=INK)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=MUTE, shrinkA=0, shrinkB=0))


def chip(ax, x, y, label, col, w):
    ax.add_patch(FancyBboxPatch((x, y), w, 4.2, boxstyle="round,pad=0,rounding_size=1.0", fc=col, ec="none", alpha=0.16))
    ax.add_patch(FancyBboxPatch((x, y), 1.1, 4.2, boxstyle="square,pad=0", fc=col, ec="none"))
    ax.text(x + 2.0, y + 2.1, label, fontsize=5.8, va="center", color=INK)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args()
    W, H = 178, 92
    fig = plt.figure(figsize=(W * MM, H * MM)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

    def step(x, n, s):
        ax.text(x, H - 3.5, n, fontsize=9, fontweight="bold", va="center", color=INK)
        ax.text(x + 4.5, H - 3.5, s, fontsize=7.5, fontweight="bold", va="center", color=INK)

    # a: sources
    step(1, "a", "Case sources")
    src = [("Phenopacket Store", "9,517 cases", "curated from full text"),
           ("PMC-Patients", "6,792 cases", "full-text patient summaries"),
           ("PubMed abstracts", "45,911 cases", "case-report abstracts"),
           ("Saudi and US cohorts", "3,978 entries", "clinical records as ranking queries")]
    y = H - 9
    for name, n, d in src:
        y -= 15.5; box(ax, 1, y, 40, 13.5)
        ax.text(3, y + 10.6, name, fontsize=7, fontweight="bold", va="center", color=INK)
        ax.text(3, y + 6.6, n, fontsize=6.5, va="center", color=INK)
        ax.text(3, y + 2.9, d, fontsize=5.8, va="center", color=MUTE)
    ax.text(1, 12, "66,198 case entries\n46,739 publications\n4,769 genes", fontsize=6, ha="left", va="center", color=MUTE, linespacing=1.5)

    # b: processing
    step(49, "b", "Phenotype profile of each case")
    box(ax, 49, H - 33, 58, 24.5, "Extraction and normalization")
    for i, t in enumerate(["Language model extracts findings stated in the text",
                           "Findings mapped to HPO classes",
                           "Causative gene, regional group, stated ancestry,",
                           "place of care and reporting group recorded"]):
        ax.text(51, H - 15.2 - 4.0 * i, t, fontsize=6, va="center", color=INK)
    box(ax, 49, 9.5, 58, 46.5, "Finding class of each HPO class")
    ax.text(51, 49.2, "Assigned once, by how the finding is established", fontsize=6, va="center", color=MUTE)
    rows = [(40.2, "Bedside", BED, ["Hypotonia", "Microcephaly"]), (35.0, "", BED, ["Feeding difficulties"]),
            (28.6, "Test-based", TEST, ["Elevated lactate", "Cerebral atrophy"]), (22.2, "Mixed", MIX, ["Seizure"])]
    for y0, lab, col, chips in rows:
        if lab: ax.text(51, y0 + 2.1, lab, fontsize=6.3, fontweight="bold", va="center", color=col if col != MIX else MUTE)
        for k, c in enumerate(chips): chip(ax, 66 + 20.5 * k, y0, c, col, 26 if len(chips) == 1 and lab == "" else 19.5)
    ax.text(87.5, 24.3, "not counted", fontsize=5.8, va="center", color=MUTE)
    ax.text(51, 18.0, "Test share of a case", fontsize=6.3, fontweight="bold", va="center", color=INK)
    ax.text(51, 13.2, r"$s = \dfrac{T}{T + B}$", fontsize=7.5, va="center", color=INK)
    ax.text(72, 14.8, r"here $T = 2$, $B = 3$,", fontsize=6, va="center", color=MUTE)
    ax.text(72, 11.4, r"$s = 0.40$", fontsize=6, va="center", color=MUTE)
    arrow(ax, 78, H - 33.4, 78, 56.6)

    # c: analyses
    step(115, "c", "Analyses")
    an = [("Regional comparison", "Figure 2", ["51,895 cases in six regional groups", "Same causative gene and corpus;", "permutation of reporting groups"]),
          ("Care setting", "Figure 3", ["Stated ancestry against place of care", "Health expenditure of 103 countries", "(41,155 cases)"]),
          ("Gene ranking", "Figure 4", ["25,234 cases ranked over 5,193 genes", "Full profile, bedside findings only,", "and a random subset of equal size"])]
    y = H - 9
    for title, figref, lines in an:
        y -= 25.7; box(ax, 115, y, 62, 23.2, title)
        ax.text(175, y + 21.0, figref, fontsize=6, ha="right", va="top", color=MUTE)
        for i, t in enumerate(lines): ax.text(117, y + 14.6 - 4.2 * i, t, fontsize=6, va="center", color=INK)

    # connectors: literature sources into extraction, profiles into the three analyses
    def line(xs, ys, **k): ax.plot(xs, ys, color=MUTE, lw=0.8, solid_capstyle="round", **k)
    ys = [H - 17.75, H - 33.25, H - 48.75]
    for yy in ys: line([41, 45], [yy, yy])
    line([45, 45], [ys[-1], ys[0]]); arrow(ax, 45, H - 21, 48.6, H - 21)
    line([107, 111], [33, 33]); ya = [H - 23.1, H - 48.8, H - 74.5]; line([111, 111], [ya[-1], ya[0]])
    for yy in ya: arrow(ax, 111, yy, 114.6, yy)
    # clinical cohorts are queries of the ranking analysis only
    line([38, 38, 146, 146], [21, 3, 3, 3.4], ls=(0, (3, 2))); arrow(ax, 146, 3, 146, 5.5)
    ax.text(92, 4.6, "clinical queries", fontsize=5.8, ha="center", va="center", color=MUTE)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out + ".pdf"); fig.savefig(a.out + ".png", dpi=220); print("wrote", a.out + ".pdf")


if __name__ == "__main__":
    main()
