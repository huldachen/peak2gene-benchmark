"""
Figure 2 - Data flow from RNA/ATAC cohorts into the peak-gene linkage methods.

Schematic of the two data sources (RNA cohort and Multiome cohort), the six
method cards (1, 2, 3, 4A, 4B, 5), and the arrows connecting them colour-coded
by integration quality.

Usage
-----
    python scripts/figures/figure2_data_flow.py
"""
from __future__ import annotations

import logging
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

ROOT    = Path(__file__).resolve().parents[2]
FIG_OUT = ROOT / "results/figures/fig2_data_flow"

COLOR_RNA    = "#4A7C59"
COLOR_ATAC   = "#C7522A"
COLOR_OK     = "#2E8B57"
COLOR_COMP   = "#E39E3A"
COLOR_FAIL   = "#B03A2E"
COLOR_PRIOR  = "#6F6FB8"


def _draw_box(ax, x, y, w, h, text, *, edge, face, text_color="black",
              fontsize=9, fontweight="normal", rounding=0.15, linewidth=2.0):
    """Draw a rounded box with centered multi-line text."""
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={rounding}",
        edgecolor=edge, facecolor=face, linewidth=linewidth,
    )
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, text,
            ha="center", va="center", fontsize=fontsize, fontweight=fontweight,
            color=text_color)
    return box


def _draw_arrow(ax, xy_from, xy_to, *, color, label=None, label_pos=0.5,
                linewidth=2.0, linestyle="-", curved=False):
    """Straight or slightly curved arrow between two points."""
    connectionstyle = "arc3,rad=0.22" if curved else "arc3,rad=0.0"
    arr = FancyArrowPatch(
        xy_from, xy_to,
        arrowstyle="-|>,head_length=10,head_width=6",
        linewidth=linewidth, color=color, linestyle=linestyle,
        connectionstyle=connectionstyle, zorder=2, mutation_scale=1.0,
    )
    ax.add_patch(arr)
    if label:
        mx = xy_from[0] + (xy_to[0] - xy_from[0]) * label_pos
        my = xy_from[1] + (xy_to[1] - xy_from[1]) * label_pos
        ax.text(mx, my, label, ha="center", va="center",
                fontsize=7.5, color=color, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white",
                          ec="none", alpha=0.85))


def build_figure():
    mpl.rcParams["font.family"] = "DejaVu Sans"
    mpl.rcParams["pdf.fonttype"] = 42

    fig, ax = plt.subplots(figsize=(16.0, 10.0))
    ax.set_xlim(0, 16); ax.set_ylim(0, 10)
    ax.axis("off")

    fig.suptitle(
        "Figure 2.  Data flow: how RNA and ATAC cohorts feed the 5 peak–gene linkage methods",
        fontsize=13, fontweight="bold", x=0.04, ha="left", y=0.97,
    )

    rna_x, rna_w = 1.5, 4.8
    mo_x,  mo_w  = 9.7, 4.8
    src_y, src_h = 7.3, 1.8

    _draw_box(ax, rna_x, src_y, rna_w, src_h,
              "RNA cohort\n"
              "B001 / B004 / B005\n"
              "11,604 cells × 45,068 genes\n"
              "(no paired ATAC)",
              edge=COLOR_RNA, face="#E8F0EB", fontsize=9, fontweight="bold")

    _draw_box(ax, mo_x, src_y, mo_w, src_h,
              "Multiome cohort\n"
              "B006 / B008–B012\n"
              "102,453 cells × 1.12 M peaks (ATAC)\n"
              "102,453 cells × 36,601 genes (paired RNA)",
              edge=COLOR_ATAC, face="#FBECE5", fontsize=9, fontweight="bold")

    ax.annotate("", xy=(mo_x - 0.05, src_y + src_h / 2),
                xytext=(rna_x + rna_w + 0.05, src_y + src_h / 2),
                arrowprops=dict(arrowstyle="<->", color="#AAA", lw=1.0,
                                linestyle="--"))
    ax.text((rna_x + rna_w + mo_x) / 2, src_y + src_h / 2 + 0.25,
            "different donors, no overlap",
            ha="center", fontsize=8, color="#888", style="italic")

    method_y, method_h = 3.0, 2.8
    method_w = 2.1
    gap = 0.25
    total_row = 6 * method_w + 5 * gap
    start_x = (16 - total_row) / 2
    method_xs = [start_x + i * (method_w + gap) for i in range(6)]

    methods = [
        {"x": method_xs[0], "name": "Method 1",
         "short": "Distance window",
         "algo": "Link every peak within\n±500 kb of any TSS",
         "color": COLOR_COMP, "result": "3.5 M links",
         "question": "Is the peak within\nthe cis-regulatory\nwindow of the gene?"},
        {"x": method_xs[1], "name": "Method 2",
         "short": "Cicero co-access",
         "algo": "Meta-cell KNN\nPearson r, r ≥ 0.25",
         "color": COLOR_COMP, "result": "599 K links",
         "question": "Do these two peaks\nopen together\n(co-accessibility)?"},
        {"x": method_xs[2], "name": "Method 3",
         "short": "ABC (power-law)",
         "algo": "Activity × Contact,\nper-gene norm. ≥ 0.02",
         "color": COLOR_COMP, "result": "111 K links",
         "question": "Does activity × distance-\nweighted contact\nfavour this gene?"},
        {"x": method_xs[3], "name": "Method 4A",
         "short": "Paired Multiome r",
         "algo": "Per-gene Pearson r\non paired meta-cells",
         "color": COLOR_OK, "result": "48 K links  ★",
         "question": "Does peak access covary\nwith gene expression\nin the same cells?"},
        {"x": method_xs[4], "name": "Method 4B",
         "short": "Cross-cohort r",
         "algo": "Gene-activity KNN →\nimputed RNA r",
         "color": COLOR_FAIL, "result": "~143 K links",
         "question": "Does peak access covary\nwith imputed RNA\nfrom a donor-mismatched ref?"},
        {"x": method_xs[5], "name": "Method 5",
         "short": "rE2G supervised",
         "algo": "8-feature logistic reg.\nthreshold 0.179",
         "color": COLOR_PRIOR, "result": "1.6 M links",
         "question": "Would a K562-trained\nsupervised model call\nthis a regulatory pair?"},
    ]

    for m in methods:
        _draw_box(ax, m["x"], method_y, method_w, method_h,
                  f"{m['name']}\n"
                  f"— {m['short']} —\n\n"
                  f"{m['algo']}\n\n"
                  f"{m['result']}",
                  edge=m["color"], face="#FAFAFA", fontsize=8.5,
                  fontweight="bold", rounding=0.12)

    def top_center(m):
        return (m["x"] + method_w / 2, method_y + method_h)

    rna_anchor_xs = np.linspace(rna_x + 0.6, rna_x + rna_w - 0.6, 2)
    mo_anchor_xs = np.linspace(mo_x + 0.5, mo_x + mo_w - 0.5, 5)

    for i in (0, 1, 2):
        src = (mo_anchor_xs[i], src_y)
        tgt = top_center(methods[i])
        _draw_arrow(ax, src, tgt, color=COLOR_ATAC, linewidth=1.4, curved=True)

    src = (mo_anchor_xs[3], src_y)
    tgt = top_center(methods[3])
    _draw_arrow(ax, src, tgt, color=COLOR_OK, linewidth=2.3, curved=False)

    _draw_arrow(ax, (rna_anchor_xs[1], src_y), top_center(methods[4]),
                color=COLOR_FAIL, linewidth=2.0, curved=True)
    _draw_arrow(ax, (mo_anchor_xs[0] + 0.0, src_y), top_center(methods[4]),
                color=COLOR_FAIL, linewidth=1.6, curved=True, linestyle="--")

    _draw_arrow(ax, (mo_anchor_xs[4], src_y), top_center(methods[5]),
                color=COLOR_ATAC, linewidth=1.4, curved=True)

    k562_x = methods[5]["x"] - 0.05
    k562_y = method_y + method_h + 0.20
    _draw_box(ax, k562_x, k562_y, method_w + 0.1, 0.55,
              "K562 LR weights\n(ENCODE-rE2G, pretrained)",
              edge=COLOR_PRIOR, face="#E8E8F8",
              fontsize=7.5, fontweight="normal", rounding=0.1, linewidth=1.3)
    _draw_arrow(ax, (k562_x + (method_w + 0.1) / 2, k562_y),
                top_center(methods[5]),
                color=COLOR_PRIOR, linewidth=1.4, linestyle=":", curved=False)

    input_labels_y = src_y - 0.28
    input_notes = [
        ("ATAC",                       COLOR_ATAC),
        ("ATAC",                       COLOR_ATAC),
        ("ATAC",                       COLOR_ATAC),
        ("RNA + ATAC\n(paired)",        COLOR_OK),
        ("RNA (ref)\n+ ATAC (target)",  COLOR_FAIL),
        ("ATAC\n+ K562 prior",          COLOR_PRIOR),
    ]
    for m, (label, col) in zip(methods, input_notes):
        ax.text(m["x"] + method_w / 2, input_labels_y,
                label, ha="center", va="top",
                fontsize=7.5, color=col, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.18", fc="white",
                          ec=col, lw=0.8, alpha=0.95))

    rib_y = 1.1
    rib_h = 1.55
    for m in methods:
        _draw_box(ax, m["x"], rib_y, method_w, rib_h,
                  m["question"],
                  edge=m["color"], face="#FFFFFF",
                  fontsize=8, fontweight="normal",
                  text_color="#222", rounding=0.1, linewidth=1.2)

    ax.text(start_x - 0.35, rib_y + rib_h / 2,
            "Question\nanswered →",
            ha="right", va="center", fontsize=9, fontweight="bold",
            color="#444")

    legend_y = 0.35
    legend_items = [
        ("★  true paired data (cell-type valid)",  COLOR_OK),
        ("⚠  cross-cohort integration fails at cell-type level (Finding E)", COLOR_FAIL),
    ]
    legend_x = 3.0
    for i, (label, col) in enumerate(legend_items):
        x = legend_x + i * 6.5
        ax.add_patch(mpatches.Rectangle((x, legend_y), 0.25, 0.22,
                                         facecolor=col, edgecolor="none"))
        ax.text(x + 0.35, legend_y + 0.11, label, fontsize=8.5, va="center",
                color="#222")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    out = FIG_OUT.with_suffix(".png")
    fig.savefig(out, dpi=220, bbox_inches="tight")
    log.info("Wrote %s", out)
    plt.close(fig)


def main():
    build_figure()


if __name__ == "__main__":
    main()
