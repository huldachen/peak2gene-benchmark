"""
Figure 3 - sLDSC benchmark: enrichment + per-SNP coefficient (tau-star).

Two-row figure:
  Row A (top):  Fold-enrichment +/- SE  (intuitive; but confounded by annotation size)
  Row B (bottom): Coefficient z-score (tau*) - unique heritability contribution
                  beyond baseline-LD v2.2, size-independent

Each row has three panels (IBD, Height, EA).  Dot size encodes prop_snps so the
reader sees annotation coverage at a glance.

Output: results/figures/fig3_ldsc_benchmark.{png}
"""
from __future__ import annotations

import logging
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
ENRICH_TSV = ROOT / "results/ldsc/enrichment_table.tsv"
FIG_OUT = ROOT / "results/figures/fig3_ldsc_benchmark"

CMAP_COMPARTMENT = {
    "epithelial": "#E4A548",
    "immune":     "#5B7EAE",
    "stromal":    "#9F5EA3",
}

METHOD_LABEL = {
    "method1_distance":         "1. Distance window",
    "method2_cicero":           "2a. Co-access (Pearson)",
    "method2b_glasso":          "2b. Co-access (GLASSO)",
    "method3_abc":              "3. ABC (power-law)",
    "method4a_paired":          "4A. Paired Multiome r",
    "method4b_crosscohort":     "4B. Cross-cohort r",
    "method5_re2g":             "5. rE2G (K562;\n  binary ≈ continuous)",
}
METHOD_ORDER = list(METHOD_LABEL.keys())
COMP_ORDER = ["epithelial", "immune", "stromal"]

# Dot-size scaling for prop_snps (annotation coverage).
# Range: smallest annotation (~0.2%) -> size 4, largest (~9%) -> size 12.
SIZE_MIN, SIZE_MAX = 4.0, 12.0
PROP_MIN, PROP_MAX = 0.002, 0.10


def _prop_to_size(prop: float) -> float:
    """Map prop_snps to marker size (linear scale, clipped)."""
    frac = (np.clip(prop, PROP_MIN, PROP_MAX) - PROP_MIN) / (PROP_MAX - PROP_MIN)
    return SIZE_MIN + frac * (SIZE_MAX - SIZE_MIN)


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


def _z_stars(z: float) -> str:
    """Significance stars based on |z| (two-sided normal)."""
    if pd.isna(z):
        return ""
    az = abs(z)
    if az >= 3.29:    # p < 0.001
        return "***"
    if az >= 2.58:    # p < 0.01
        return "**"
    if az >= 1.96:    # p < 0.05
        return "*"
    return ""


def _draw_panel(ax, rows, df, trait, metric, y_pos, n_rows, panel_i,
                show_ylabels=True):
    """Draw one panel (either enrichment or coefficient_z)."""
    band_colors = ["#F6F6F6", "#FFFFFF"]

    # Alternating method-group bands
    for mi, m in enumerate(METHOD_ORDER):
        mask = (rows["method"] == m).values
        top_y = y_pos[mask].max() + 0.5
        bot_y = y_pos[mask].min() - 0.5
        ax.axhspan(bot_y, top_y, color=band_colors[mi % 2], zorder=0)

    # Null reference line
    null_val = 1.0 if metric == "enrichment" else 0.0
    ax.axvline(null_val, color="#555", linestyle="--", linewidth=0.9, zorder=1)

    # Gather values
    vals      = np.full(n_rows, np.nan)
    ses       = np.full(n_rows, np.nan)
    sig_vals  = np.full(n_rows, np.nan)  # p-value or z-score for stars
    prop_snps = np.full(n_rows, np.nan)
    comps_of_row = rows["compartment"].values

    for i, (_, r) in enumerate(rows.iterrows()):
        sub = df[(df["method"] == r["method"])
                 & (df["compartment"] == r["compartment"])
                 & (df["trait"] == trait)]
        if sub.empty:
            continue
        row_data = sub.iloc[0]
        prop_snps[i] = float(row_data["prop_snps"])
        if metric == "enrichment":
            vals[i] = float(row_data["enrichment"])
            ses[i]  = float(row_data["enrichment_se"])
            sig_vals[i] = float(row_data["enrichment_p"])
        else:
            vals[i] = float(row_data["coefficient_z"])
            ses[i]  = float(row_data["coefficient_se"])
            sig_vals[i] = float(row_data["coefficient_z"])

    # Plot points
    for i, (v, s, sig, ps, comp) in enumerate(
        zip(vals, ses, sig_vals, prop_snps, comps_of_row)
    ):
        if np.isnan(v):
            continue
        col = CMAP_COMPARTMENT[comp]
        ms = _prop_to_size(ps) if not np.isnan(ps) else 6.0

        if metric == "enrichment":
            ax.errorbar(v, y_pos[i], xerr=s,
                        fmt="o", markersize=ms, color=col,
                        ecolor=col, elinewidth=1.5, capsize=3.0, zorder=3,
                        markeredgecolor="black", markeredgewidth=0.5)
            stars = _stars(sig)
            offset = s + 0.5
            label = f"{v:.1f}{stars}"
        else:
            # For tau*, no SE error bar - just plot the z-score point
            ax.plot(v, y_pos[i], "o", markersize=ms, color=col,
                    markeredgecolor="black", markeredgewidth=0.5, zorder=3)
            stars = _z_stars(sig)
            # Place label to the right of positive values, left of negative
            if v >= 0:
                offset = 0.15
                ha_text = "left"
            else:
                offset = -0.15
                ha_text = "right"
            label = f"{v:.2f}{stars}"

        if metric == "enrichment":
            ax.text(v + offset, y_pos[i], label,
                    va="center", ha="left", fontsize=7.5,
                    color="#222", zorder=4)
        else:
            ax.text(v + offset, y_pos[i], label,
                    va="center", ha=ha_text, fontsize=7.0,
                    color="#222", zorder=4)

    # Axis formatting
    ax.set_ylim(-1, n_rows)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if metric == "enrichment":
        xmax = max(float(np.nanmax(vals) + np.nanmax(ses) + 5), 8)
        xmin = min(float(np.nanmin(vals) - np.nanmax(ses) - 2), -5)
        ax.set_xlim(xmin, xmax + 1)
        ax.set_xlabel("Enrichment (fold)", fontsize=10)
    else:
        abs_max = max(float(np.nanmax(np.abs(vals[~np.isnan(vals)]))) + 1.5, 5)
        ax.set_xlim(-abs_max, abs_max)
        ax.set_xlabel("Coefficient z-score (tau*)", fontsize=10)

    if show_ylabels and panel_i == 0:
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
        ax.set_yticks(y_pos)
        ax.set_yticklabels([])


def build_figure():
    mpl.rcParams["font.family"] = "DejaVu Sans"
    df = pd.read_csv(ENRICH_TSV, sep="\t")
    df = df[df["method"].isin(METHOD_ORDER)].copy()

    # Build row labels: method x compartment
    row_meta = []
    for m in METHOD_ORDER:
        for c in COMP_ORDER:
            row_meta.append((m, c))
    rows = pd.DataFrame(row_meta, columns=["method", "compartment"])
    rows["method_label"] = rows["method"].map(METHOD_LABEL)

    traits = ["IBD", "Height", "EA"]
    trait_subtitle = {
        "IBD":    "IBD  (primary - gut-relevant)",
        "Height": "Height  (negative ctrl - expect stromal)",
        "EA":     "Education Yrs  (negative ctrl - no signal)",
    }

    n_rows = len(rows)
    y_pos = np.arange(n_rows)[::-1]

    # Two-row layout: top = enrichment, bottom = coefficient_z
    fig = plt.figure(figsize=(16.5, 16.0))
    gs = fig.add_gridspec(2, 3, wspace=0.05, hspace=0.25,
                          left=0.28, right=0.985, top=0.93, bottom=0.07)

    # --- Row A: Fold-enrichment ---
    for panel_i, trait in enumerate(traits):
        ax = fig.add_subplot(gs[0, panel_i])
        _draw_panel(ax, rows, df, trait, "enrichment", y_pos, n_rows, panel_i)
        ax.set_title(trait_subtitle[trait], fontsize=10.5, fontweight="bold",
                     loc="left", pad=8)

    # --- Row B: Coefficient z-score (tau*) ---
    for panel_i, trait in enumerate(traits):
        ax = fig.add_subplot(gs[1, panel_i])
        _draw_panel(ax, rows, df, trait, "coefficient_z", y_pos, n_rows, panel_i)
        ax.set_title(trait_subtitle[trait], fontsize=10.5, fontweight="bold",
                     loc="left", pad=8)

    # Row labels - centred over the three panels (midpoint of left..right)
    row_label_x = 0.28 + (0.985 - 0.28) / 2   # centre of gridspec
    fig.text(row_label_x, 0.945,
             "A.  Fold-enrichment  (prop_h₂ / prop_SNPs)",
             fontsize=12, fontweight="bold", ha="center")
    fig.text(row_label_x, 0.505,
             "B.  Per-SNP coefficient z-score  (tau*  |  unique signal beyond baseline-LD v2.2)",
             fontsize=12, fontweight="bold", ha="center")

    # Suptitle
    fig.suptitle(
        "Figure 3.  sLDSC partitioned heritability: enrichment and unique contribution (tau*)",
        fontsize=13, fontweight="bold", x=0.02, ha="left", y=0.97,
    )

    # Compartment legend + size legend
    legend_handles = [
        plt.Line2D([], [], marker="o", color="white",
                   markerfacecolor=CMAP_COMPARTMENT[c],
                   markersize=9, label=c,
                   markeredgecolor="black", markeredgewidth=0.5)
        for c in COMP_ORDER
    ]
    # Add size legend entries
    for pct_label, pct_val in [("~0.3% SNPs", 0.003), ("~3% SNPs", 0.03),
                                ("~9% SNPs", 0.09)]:
        legend_handles.append(
            plt.Line2D([], [], marker="o", color="white",
                       markerfacecolor="#999",
                       markersize=_prop_to_size(pct_val),
                       label=pct_label,
                       markeredgecolor="black", markeredgewidth=0.5)
        )
    fig.legend(handles=legend_handles, loc="upper right", ncol=6,
               bbox_to_anchor=(0.985, 0.99), frameon=False, fontsize=9.5,
               title="Compartment                          Annotation size",
               title_fontsize=10,
               columnspacing=1.0, handletextpad=0.3)

    # Caption
    fig.text(0.5, 0.02,
             "Row A: *** p <= 0.001 * ** p <= 0.01 * * p <= 0.05  "
             "(enrichment p from LDSC)     "
             "Row B: *** |z| >= 3.29 * ** |z| >= 2.58 * * |z| >= 1.96     "
             "Dot size = annotation coverage (prop_SNPs)     "
             "Dashed line = null",
             fontsize=8.5, color="#555", ha="center", va="center")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    out = FIG_OUT.with_suffix(".png")
    fig.savefig(out, dpi=220, bbox_inches="tight")
    log.info("Wrote %s", out)
    plt.close(fig)


if __name__ == "__main__":
    build_figure()
