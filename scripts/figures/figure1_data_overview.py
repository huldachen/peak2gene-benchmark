"""
Figure 1 - Dataset overview + QC.

Multi-panel figure covering cohort structure, cell composition, RNA UMAP,
and QC metrics for the RNA and Multiome cohorts.

Usage
-----
    python scripts/figures/figure1_data_overview.py
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

ROOT         = Path(__file__).resolve().parents[2]
RNA_H5AD     = ROOT / "data/processed/hickey2023/rna_qc.h5ad"
ATAC_MERGED  = ROOT / "data/processed/hickey2023/atac_merged.h5ad"
MULTIOME_RNA = ROOT / "data/processed/hickey2023/multiome_rna_merged.h5ad"

FIG_OUT = ROOT / "results/figures/fig1_data_overview"

CMAP_COHORT = {"RNA cohort": "#4A7C59", "Multiome cohort": "#C7522A"}
CMAP_COMPARTMENT = {
    "epithelial": "#E4A548",
    "colon_epithelial": "#E4A548",
    "immune": "#5B7EAE",
    "stromal": "#9F5EA3",
}


def _load_data():
    """Load RNA h5ad + pull ATAC obs from the three compartment h5ads."""
    log.info("Loading h5ad files ...")
    rna = sc.read_h5ad(str(RNA_H5AD))

    atac_obs_pieces = []
    for fn in ("atac_colon_epithelial.h5ad", "atac_immune.h5ad", "atac_stromal.h5ad"):
        p = ROOT / "data/processed/hickey2023" / fn
        if not p.exists():
            continue
        df = ad.read_h5ad(str(p), backed="r").obs.copy()
        atac_obs_pieces.append(df)
    atac_obs = pd.concat(atac_obs_pieces, axis=0) if atac_obs_pieces else pd.DataFrame()

    log.info("RNA:  %d cells x %d genes", rna.n_obs, rna.n_vars)
    log.info("ATAC: %d cells from compartment h5ads (cols: %s)",
             len(atac_obs), list(atac_obs.columns))
    return rna, atac_obs


def _panel_A_cohort_schematic(ax):
    """Two labelled boxes showing the two-cohort structure, with an outer panel frame."""
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.axis("off")

    panel_frame = mpatches.FancyBboxPatch(
        (0.15, 0.7), 9.7, 9.0,
        boxstyle="round,pad=0.02,rounding_size=0.15",
        edgecolor="#BBBBBB", facecolor="#FAFAFA", linewidth=1.0, zorder=0,
    )
    ax.add_patch(panel_frame)

    ax.text(0.35, 9.2, "a.  Two-cohort dataset structure",
            ha="left", va="center", fontweight="bold", fontsize=12, color="#222")

    rna_box = mpatches.FancyBboxPatch(
        (0.55, 2.8), 4.1, 5.6,
        boxstyle="round,pad=0.1,rounding_size=0.2",
        edgecolor=CMAP_COHORT["RNA cohort"], facecolor="#E8F0EB",
        linewidth=2.0, zorder=1,
    )
    ax.add_patch(rna_box)
    ax.text(2.6, 7.9, "RNA cohort",
            ha="center", fontsize=11, fontweight="bold",
            color=CMAP_COHORT["RNA cohort"])
    ax.text(2.6, 7.25, "scRNA-seq only", ha="center", fontsize=9, style="italic")
    ax.text(2.6, 6.7, "3 donors  ·  4 colon samples", ha="center", fontsize=9)
    ax.text(2.6, 6.25, "11,604 cells × 45,068 genes", ha="center", fontsize=9)
    ax.text(2.6, 5.4, "• B001 (Sigmoid, Transverse)",
            ha="center", fontsize=8, color="#333")
    ax.text(2.6, 5.0, "• B004 (Ascending)", ha="center", fontsize=8, color="#333")
    ax.text(2.6, 4.6, "• B005 (Ascending)", ha="center", fontsize=8, color="#333")
    ax.text(2.6, 3.7, "modalities: RNA",
            ha="center", fontsize=9, fontweight="bold", color="#555")

    mo_box = mpatches.FancyBboxPatch(
        (5.35, 2.8), 4.1, 5.6,
        boxstyle="round,pad=0.1,rounding_size=0.2",
        edgecolor=CMAP_COHORT["Multiome cohort"], facecolor="#FBECE5",
        linewidth=2.0, zorder=1,
    )
    ax.add_patch(mo_box)
    ax.text(7.4, 7.9, "Multiome cohort",
            ha="center", fontsize=11, fontweight="bold",
            color=CMAP_COHORT["Multiome cohort"])
    ax.text(7.4, 7.25, "10x Multiome (paired RNA + ATAC)",
            ha="center", fontsize=9, style="italic")
    ax.text(7.4, 6.7, "6 donors  ·  42 samples", ha="center", fontsize=9)
    ax.text(7.4, 6.25, "102,453 cells × 1.12 M peaks", ha="center", fontsize=9)
    ax.text(7.4, 5.7, "+ 102,453 cells × 36,601 genes (paired RNA)",
            ha="center", fontsize=8.5, color="#666")
    ax.text(7.4, 4.9, "B006 / B008 / B009 / B010 / B011 / B012",
            ha="center", fontsize=8, color="#333")
    ax.text(7.4, 4.5, "colon (21 samples) + small intestine (20)",
            ha="center", fontsize=8, color="#333")
    ax.text(7.4, 3.7, "modalities: RNA + ATAC (paired)",
            ha="center", fontsize=9, fontweight="bold", color="#555")

    ax.text(5.0, 2.2, "— different donors, no overlap —",
            ha="center", fontsize=8.5, style="italic", color="#888")

    ax.text(5.0, 1.1,
            "Source: Hickey et al. 2023 Nature  ·  doi:10.5061/dryad.8pk0p2ns8 / 0zpc8672f",
            ha="center", fontsize=7, color="#999")


def _panel_B_cell_composition(ax, rna, atac_obs):
    """Stacked bars: cells per compartment × cohort."""
    rna_comp = rna.obs["compartment"].value_counts()
    atac_comp = atac_obs["compartment"].value_counts()

    # ATAC uses "colon_epithelial"; normalise to "epithelial"
    rna_comp_clean = pd.Series({
        "epithelial": rna_comp.get("colon_epithelial", 0) + rna_comp.get("epithelial", 0),
        "immune":     rna_comp.get("immune", 0),
        "stromal":    rna_comp.get("stromal", 0),
    })
    atac_comp_clean = pd.Series({
        "epithelial": atac_comp.get("colon_epithelial", 0) + atac_comp.get("epithelial", 0),
        "immune":     atac_comp.get("immune", 0),
        "stromal":    atac_comp.get("stromal", 0),
    })

    comps = ["epithelial", "immune", "stromal"]
    cohorts = ["RNA cohort", "Multiome cohort"]
    data = np.array([
        [rna_comp_clean[c] for c in comps],
        [atac_comp_clean[c] for c in comps],
    ])

    x = np.arange(len(cohorts))
    bottom = np.zeros(len(cohorts))
    for i, comp in enumerate(comps):
        ax.bar(x, data[:, i], bottom=bottom, label=comp,
               color=CMAP_COMPARTMENT[comp], edgecolor="white", linewidth=0.5)
        for j, val in enumerate(data[:, i]):
            if val > 3000:
                ax.text(x[j], bottom[j] + val / 2, f"{val:,}",
                        ha="center", va="center", fontsize=7.5, color="white",
                        fontweight="bold")
        bottom += data[:, i]

    ax.set_xticks(x); ax.set_xticklabels(cohorts, fontsize=9)
    ax.set_ylabel("Cell count", fontsize=9)
    ax.set_title("b.  Cell composition by cohort × compartment",
                 loc="left", fontweight="bold", fontsize=12)
    for j, tot in enumerate(data.sum(axis=1)):
        ax.text(x[j], tot + tot * 0.02, f"Σ = {tot:,}",
                ha="center", fontsize=8, color="#555")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, data.sum(axis=1).max() * 1.12)
    ax.yaxis.set_major_formatter(
        mpl.ticker.FuncFormatter(lambda v, _: f"{int(v):,}")
    )


def _panel_C_rna_umap(ax, rna):
    """UMAP of RNA cells coloured by compartment."""
    if "X_umap" not in rna.obsm:
        log.warning("RNA UMAP not in rna.obsm - computing on-the-fly")
        sc.pp.neighbors(rna, n_neighbors=15, use_rep="X_pca")
        sc.tl.umap(rna, random_state=42)

    umap = rna.obsm["X_umap"]
    comp = rna.obs["compartment"].astype(str).fillna("unknown")
    comp = comp.replace({"colon_epithelial": "epithelial"})
    comp_order = ["epithelial", "immune", "stromal"]

    ax.set_aspect("equal")
    for c in comp_order:
        m = (comp == c).values
        if not m.any():
            continue
        ax.scatter(umap[m, 0], umap[m, 1], s=1.8,
                   c=CMAP_COMPARTMENT[c], alpha=0.55,
                   rasterized=True, linewidths=0, label=f"{c}  (n={m.sum():,})")

    ax.set_xlabel("UMAP 1", fontsize=9)
    ax.set_ylabel("UMAP 2", fontsize=9)
    n_ct = rna.obs["cell_type"].nunique()
    ax.set_title(f"c.  RNA cohort UMAP  ({rna.n_obs:,} cells · {n_ct} cell types)",
                 loc="left", fontweight="bold", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.legend(loc="lower left", fontsize=9, frameon=False,
              title="Compartment", title_fontsize=9,
              markerscale=2.5)


def _panel_D_qc_violins(gs_row, fig, rna, atac_obs):
    """QC metrics as 4 mini-violins with individual y-axes."""
    rna_ng = rna.obs["n_genes_by_counts"].values
    rna_mt = rna.obs["pct_counts_mt"].values
    atac_np = atac_obs["n_accessible_peaks"].values if "n_accessible_peaks" in atac_obs.columns else np.array([1])
    atac_ta = atac_obs["total_accessibility"].values if "total_accessibility" in atac_obs.columns else np.array([1])

    specs = [
        ("RNA  ·  n_genes",           "genes detected / cell",     rna_ng,                CMAP_COHORT["RNA cohort"],      False),
        ("RNA  ·  % mitochondrial",    "% of total counts",          rna_mt,                CMAP_COHORT["RNA cohort"],      False),
        ("ATAC  ·  n peaks accessible","log₁₀(peaks / cell)",       atac_np.astype(float), CMAP_COHORT["Multiome cohort"], True),
        ("ATAC  ·  total accessibility", "log₁₀(fragments / cell)",   atac_ta.astype(float), CMAP_COHORT["Multiome cohort"], True),
    ]

    sub_gs = gs_row.subgridspec(1, 4, wspace=0.70)
    axs = [fig.add_subplot(sub_gs[0, i]) for i in range(4)]

    axs[0].text(-0.45, 1.15, "d.  QC metrics (per cell)",
                transform=axs[0].transAxes, fontweight="bold", fontsize=12,
                ha="left")

    for i, (title, y_label, d, col, use_log) in enumerate(specs):
        ax = axs[i]
        if d.size == 0:
            ax.text(0.5, 0.5, "(no data)", ha="center", va="center")
            ax.axis("off"); continue
        if use_log:
            d = np.log10(d + 1)
        parts = ax.violinplot([d], positions=[0], showmeans=False,
                              showmedians=True, widths=0.75)
        for pc in parts["bodies"]:
            pc.set_facecolor(col); pc.set_alpha(0.65)
            pc.set_edgecolor("black"); pc.set_linewidth(0.6)
        for partname in ("cbars", "cmins", "cmaxes", "cmedians"):
            parts[partname].set_edgecolor("black")
            parts[partname].set_linewidth(0.8)
        med = float(np.median(d))
        ax.text(0.35, med, f"median\n{med:.2f}",
                ha="left", va="center", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", fc="white",
                          ec="#888", lw=0.5, alpha=0.9))
        ax.set_xticks([])
        ax.set_title(title, fontsize=9, fontweight="normal", pad=4)
        ax.set_ylabel(y_label, fontsize=8.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(-0.8, 0.8)


def build_figure():
    mpl.rcParams["font.family"] = "DejaVu Sans"
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    rna, atac_obs = _load_data()

    fig = plt.figure(figsize=(15.0, 10.0), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.28,
                          left=0.06, right=0.95, top=0.92, bottom=0.07,
                          height_ratios=[1.0, 1.0])

    ax_A = fig.add_subplot(gs[0, 0])
    ax_B = fig.add_subplot(gs[0, 1])
    ax_C = fig.add_subplot(gs[1, 0])

    _panel_A_cohort_schematic(ax_A)
    _panel_B_cell_composition(ax_B, rna, atac_obs)
    _panel_C_rna_umap(ax_C, rna)
    _panel_D_qc_violins(gs[1, 1], fig, rna, atac_obs)

    fig.suptitle(
        "Figure 1.  Hickey 2023 colon single-cell dataset: "
        "two-cohort architecture + QC",
        fontsize=13, fontweight="bold", y=0.975, x=0.06, ha="left",
    )

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    out = FIG_OUT.with_suffix(".png")
    fig.savefig(out, dpi=220, bbox_inches="tight")
    log.info("Wrote %s", out)
    plt.close(fig)


def main():
    build_figure()


if __name__ == "__main__":
    main()
