"""
Figure 3 - sLDSC partitioned-heritability enrichment forest plot.

Three side-by-side panels (IBD, Height, EA) showing enrichment +/- SE for
every (method x compartment) combination, annotated with significance stars.

Usage
-----
    python scripts/figures/figure3_ldsc_benchmark.py
"""
from __future__ import annotations

import logging
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

ROOT       = Path(__file__).resolve().parents[2]
ENRICH_TSV = ROOT / "results/ldsc/enrichment_table.tsv"
FIG_OUT    = ROOT / "results/figures/fig3_ldsc_benchmark"

CMAP_COMPARTMENT = {
    "epithelial": "#E4A548",
    "immune":     "#5B7EAE",
    "stromal":    "#9F5EA3",
}

METHOD_LABEL = {
    "method1_distance":         "1. Distance window",
    "method2_cicero":           "2. Cicero co-access",
    "method3_abc":              "3. ABC (power-law)",
    "method4a_paired":          "4A. Paired Multiome r",
    "method4b_crosscohort":     "4B. Cross-cohort r",
    "method5_re2g":             "5. rE2G (K562;\n  binary ≈ continuous)",
}
METHOD_ORDER = list(METHOD_LABEL.keys())
COMP_ORDER = ["epithelial", "immune", "stromal"]


def _stars(p: float) -> str:
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def build_figure():
    mpl.rcParams["font.family"] = "DejaVu Sans"
    mpl.rcParams["pdf.fonttype"] = 42

    df = pd.read_csv(ENRICH_TSV, sep="\t")
    df = df[df["method"].isin(METHOD_ORDER)].copy()

    row_meta = []
    for m in METHOD_ORDER:
        for c in COMP_ORDER:
            row_meta.append((m, c))
    rows = pd.DataFrame(row_meta, columns=["method", "compartment"])
    rows["method_label"] = rows["method"].map(METHOD_LABEL)

    traits = ["IBD", "Height", "EA"]
    trait_subtitle = {
        "IBD":    "IBD  (primary trait — gut-relevant)",
        "Height": "Height  (negative control — expect stromal)",
        "EA":     "Education Years  (negative control — expect no signal)",
    }

    fig = plt.figure(figsize=(16.5, 10.0))
    gs = fig.add_gridspec(1, 3, wspace=0.05,
                          left=0.28, right=0.985, top=0.91, bottom=0.13)

    n_rows = len(rows)
    y_pos = np.arange(n_rows)[::-1]

    band_colors = ["#F6F6F6", "#FFFFFF"]

    axes = []
    for panel_i, trait in enumerate(traits):
        ax = fig.add_subplot(gs[0, panel_i])
        axes.append(ax)

        for mi, m in enumerate(METHOD_ORDER):
            mask = (rows["method"] == m).values
            top_y = y_pos[mask].max() + 0.5
            bot_y = y_pos[mask].min() - 0.5
            ax.axhspan(bot_y, top_y, color=band_colors[mi % 2], zorder=0)

        ax.axvline(1.0, color="#555", linestyle="--", linewidth=0.9, zorder=1)

        enrich = np.full(n_rows, np.nan)
        se     = np.full(n_rows, np.nan)
        pvals  = np.full(n_rows, np.nan)
        comps_of_row = rows["compartment"].values

        for i, (_, r) in enumerate(rows.iterrows()):
            sub = df[(df["method"] == r["method"])
                     & (df["compartment"] == r["compartment"])
                     & (df["trait"] == trait)]
            if not sub.empty:
                enrich[i] = float(sub.iloc[0]["enrichment"])
                se[i]     = float(sub.iloc[0]["enrichment_se"])
                pvals[i]  = float(sub.iloc[0]["enrichment_p"])

        for i, (e, s, p, comp) in enumerate(zip(enrich, se, pvals, comps_of_row)):
            if np.isnan(e):
                continue
            col = CMAP_COMPARTMENT[comp]
            ax.errorbar(e, y_pos[i], xerr=s,
                        fmt="o", markersize=6, color=col,
                        ecolor=col, elinewidth=1.5, capsize=3.0, zorder=3,
                        markeredgecolor="black", markeredgewidth=0.5)
            stars = _stars(p)
            ax.text(e + s + 0.5, y_pos[i],
                    f"{e:.1f}{stars}",
                    va="center", ha="left", fontsize=7.5,
                    color="#222", zorder=4)

        ax.set_title(trait_subtitle[trait], fontsize=11, fontweight="bold",
                     loc="left", pad=8)
        ax.set_xlabel("Enrichment (fold)", fontsize=10)
        ax.set_ylim(-1, n_rows)
        xmax = max(float(np.nanmax(enrich) + np.nanmax(se) + 5), 8)
        xmin = min(float(np.nanmin(enrich) - np.nanmax(se) - 2), -5)
        ax.set_xlim(xmin, xmax + 1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if panel_i == 0:
            labels = [r["compartment"] for _, r in rows.iterrows()]
            ax.set_yticks(y_pos)
            ax.set_yticklabels(labels, fontsize=8.5)
            for tick_label, (_, r) in zip(ax.get_yticklabels(), rows.iterrows()):
                tick_label.set_color(CMAP_COMPARTMENT[r["compartment"]])

            for mi, m in enumerate(METHOD_ORDER):
                mask = (rows["method"] == m).values
                mid_y = float(y_pos[mask].mean())
                ax.text(-0.35, mid_y, METHOD_LABEL[m],
                        transform=ax.get_yaxis_transform(),
                        ha="right", va="center", fontsize=9.5,
                        fontweight="bold", color="#333")
        else:
            ax.set_yticks(y_pos); ax.set_yticklabels([])

    fig.suptitle(
        "Figure 3.  sLDSC partitioned-heritability enrichment across "
        "5 peak–gene linkage methods",
        fontsize=13, fontweight="bold", x=0.02, ha="left", y=0.97,
    )

    legend_handles = [
        plt.Line2D([], [], marker="o", color="white",
                   markerfacecolor=CMAP_COMPARTMENT[c],
                   markersize=9, label=c,
                   markeredgecolor="black", markeredgewidth=0.5)
        for c in COMP_ORDER
    ]
    fig.legend(handles=legend_handles, loc="upper right", ncol=3,
               bbox_to_anchor=(0.985, 0.99), frameon=False, fontsize=10,
               title="Compartment", title_fontsize=10,
               columnspacing=1.2, handletextpad=0.4)

    fig.text(0.5, 0.035,
             "Stars: *** p ≤ 0.001  ·  ** p ≤ 0.01  ·  * p ≤ 0.05     "
             "Dashed line at 1× = no enrichment (null)     "
             "Error bars = ± 1 SE (LDSC block-jackknife)",
             fontsize=9, color="#555", ha="center", va="center")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    out = FIG_OUT.with_suffix(".png")
    fig.savefig(out, dpi=220, bbox_inches="tight")
    log.info("Wrote %s", out)
    plt.close(fig)


def main():
    build_figure()


if __name__ == "__main__":
    main()
