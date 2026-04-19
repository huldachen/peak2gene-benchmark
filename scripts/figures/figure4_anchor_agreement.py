"""
Figure 4 - Cross-cohort anchor-agreement diagnostic.

Horizontal bar chart per compartment showing, per cell type, the mean
fraction of each ATAC cell's top-30 RNA anchors that share its cell-type
label.

Usage
-----
    python scripts/figures/figure4_anchor_agreement.py
"""
from __future__ import annotations

import logging
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mp_patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

ROOT       = Path(__file__).resolve().parents[2]
ANCHOR_DIR = ROOT / "results/linkage/method4b_anchor_diagnostic"
FIG_OUT    = ROOT / "results/figures/fig4_anchor_agreement"

COMPARTMENTS = ["epithelial", "immune", "stromal"]


def _agreement_color(agreement: float) -> str:
    if agreement >= 0.50:
        return "#2E8B57"
    if agreement >= 0.10:
        return "#E39E3A"
    return "#B03A2E"


def build_figure():
    mpl.rcParams["font.family"] = "DejaVu Sans"
    mpl.rcParams["pdf.fonttype"] = 42

    fig = plt.figure(figsize=(14.0, 11.5))
    gs = fig.add_gridspec(3, 1, hspace=0.55,
                          left=0.22, right=0.92, top=0.915, bottom=0.09,
                          height_ratios=[1.4, 1.2, 1.5])

    for panel_i, comp in enumerate(COMPARTMENTS):
        ax = fig.add_subplot(gs[panel_i, 0])

        f = ANCHOR_DIR / f"{comp}_per_cell_type.tsv"
        if not f.exists():
            ax.text(0.5, 0.5, f"(missing: {f})", ha="center", va="center")
            continue
        df = pd.read_csv(f, sep="\t")

        df = df.sort_values("mean_agreement", ascending=True).reset_index(drop=True)
        # n >= 50: anchor stats are unreliable below this
        df = df[df["n"] >= 50].reset_index(drop=True)

        y_pos = np.arange(len(df))
        bar_colors = [_agreement_color(v) for v in df["mean_agreement"]]

        ax.barh(y_pos, df["mean_agreement"], color=bar_colors,
                edgecolor="black", linewidth=0.5, height=0.75, zorder=2)

        for i, row in df.iterrows():
            n = int(row["n"])
            agr = float(row["mean_agreement"])
            label = f"{agr*100:.1f}%  (n = {n:,})"
            ax.text(max(agr, 0.02) + 0.01, i, label,
                    va="center", ha="left", fontsize=8, color="#333")

        ax.axvline(0.25, color="#888", linestyle=":", linewidth=0.9, zorder=1)
        ax.axvline(0.50, color="#444", linestyle="--", linewidth=0.9, zorder=1)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(df["cell_type_atac"].values, fontsize=8.5)
        for tick_label, col in zip(ax.get_yticklabels(), bar_colors):
            tick_label.set_color(col)
            if col == "#2E8B57":
                tick_label.set_fontweight("bold")

        ax.set_xlim(0, 1.12)
        ax.set_xticks(np.arange(0, 1.01, 0.1))
        ax.set_xticklabels([f"{v:.0%}" for v in np.arange(0, 1.01, 0.1)], fontsize=8)
        if panel_i == len(COMPARTMENTS) - 1:
            ax.set_xlabel("Mean anchor-agreement  (fraction of top-30 RNA anchors sharing cell-type label)",
                          fontsize=10, labelpad=6)

        weighted_mean = float((df["mean_agreement"] * df["n"]).sum() / df["n"].sum())
        n_total = int(df["n"].sum())
        ax.set_title(
            f"{chr(97+panel_i)}.  {comp}  —  {len(df)} cell types, "
            f"{n_total:,} cells  ·  overall mean agreement: {weighted_mean*100:.1f}%",
            loc="left", fontweight="bold", fontsize=11, pad=4,
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(
        "Figure 4.  Cross-cohort anchor-agreement per cell type  "
        "(Finding E: integration fails at cell-type resolution)",
        fontsize=13, fontweight="bold", x=0.02, ha="left", y=0.985,
    )

    y_cap = 0.025
    entries = [
        ("#2E8B57", "green  ≥ 50% anchors correctly typed"),
        ("#E39E3A", "amber  10–50% marginal"),
        ("#B03A2E", "red  < 10% (most cell types)"),
    ]
    x_cap = 0.03
    for col, text in entries:
        fig.patches.append(
            mp_patches.Rectangle((x_cap, y_cap - 0.005), 0.014, 0.018,
                                 transform=fig.transFigure,
                                 facecolor=col, edgecolor="black",
                                 linewidth=0.5, zorder=5)
        )
        fig.text(x_cap + 0.017, y_cap + 0.004, text,
                 fontsize=8.5, color="#333", ha="left", va="center")
        x_cap += 0.017 + 0.005 + 0.008 * len(text)
    fig.text(0.965, y_cap + 0.004,
             "dashed @ 50% = ‘winner’ cut · dotted @ 25% = ‘partial’ cut",
             fontsize=8.5, color="#555", ha="right", va="center")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    out = FIG_OUT.with_suffix(".png")
    fig.savefig(out, dpi=220, bbox_inches="tight")
    log.info("Wrote %s", out)
    plt.close(fig)


def main():
    build_figure()


if __name__ == "__main__":
    main()
