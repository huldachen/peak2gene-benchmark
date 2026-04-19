"""
Method 4B diagnostic: anchor-transfer cell-type agreement per compartment.

For each ATAC cell, finds its top-K RNA anchors (same gene-activity KNN as
Method 4B) and measures what fraction share the ATAC cell's cell-type label.
Writes a per-cell TSV plus a per-cell-type summary to
``results/linkage/method4b_anchor_diagnostic/{compartment}{,_per_cell_type}.tsv``.

Usage
-----
    python scripts/linkage/method4b_anchor_diagnostic.py --compartment immune
    python scripts/linkage/method4b_anchor_diagnostic.py --all
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parent))
from method4_crosscohort import (
    CANONICAL_CHROMS,
    DEFAULT_ATAC_H5AD,
    _normalise_log1p,
    compute_gene_activity,
)
from method4_paired import extract_tss_frame, _peaks_from_h5ad

log = logging.getLogger(__name__)

ROOT     = Path(__file__).resolve().parents[2]
GTF      = ROOT / "data/raw/reference/gencode.v44.annotation.gtf.gz"
RNA_H5AD = ROOT / "data/processed/hickey2023/rna_merged.h5ad"
OUT_DIR  = ROOT / "results/linkage/method4b_anchor_diagnostic"


def diagnose_compartment(
    compartment: str,
    rna_h5ad: Path = RNA_H5AD,
    atac_h5ad: Optional[Path] = None,
    n_anchors: int = 30,
    upstream_bp: int = 2_000,
    gene_body_bp: int = 100_000,
) -> pd.DataFrame:
    """Run anchor-agreement diagnostic for one compartment."""
    if atac_h5ad is None:
        atac_h5ad = Path(DEFAULT_ATAC_H5AD[compartment])
    atac_h5ad = Path(atac_h5ad)

    t0 = time.time()
    log.info("[%s] Loading RNA  = %s", compartment, rna_h5ad)
    log.info("[%s] Loading ATAC = %s", compartment, atac_h5ad)
    rna  = sc.read_h5ad(str(rna_h5ad))
    atac = sc.read_h5ad(str(atac_h5ad))

    if "cell_type" not in atac.obs.columns:
        raise ValueError("ATAC h5ad missing cell_type column")
    if "cell_type" not in rna.obs.columns:
        raise ValueError("RNA h5ad missing cell_type column")

    peaks_df = _peaks_from_h5ad(atac)
    keep = peaks_df["chrom"].isin(CANONICAL_CHROMS).values
    atac = atac[:, keep].copy()
    peaks_df = peaks_df.loc[keep].reset_index(drop=True)
    tss_df = extract_tss_frame(str(GTF))

    if "gene_ids" not in rna.var.columns:
        tss_df = tss_df.assign(gene_id=tss_df["gene_name"])

    activity, gene_order = compute_gene_activity(
        atac, tss_df, peaks_df,
        upstream_bp=upstream_bp, gene_body_bp=gene_body_bp,
    )

    if "gene_ids" in rna.var.columns:
        rna_gene_ids = rna.var["gene_ids"].astype(str).values
    else:
        rna_gene_ids = rna.var_names.astype(str).values

    rna_gene_idx = {g: i for i, g in enumerate(rna_gene_ids)}
    gene_activity_idx = {g: i for i, g in enumerate(gene_order)}
    common = [g for g in gene_order if g in rna_gene_idx]
    log.info("[%s] bridge feature space: %d common genes", compartment, len(common))

    act_cols = [gene_activity_idx[g] for g in common]
    rna_cols = [rna_gene_idx[g] for g in common]
    atac_feat = activity[:, act_cols]

    rna_raw = rna.X
    if sp.issparse(rna_raw):
        rna_raw = rna_raw.toarray()
    rna_feat = rna_raw[:, rna_cols]

    atac_feat_norm = _normalise_log1p(atac_feat)
    rna_feat_norm  = _normalise_log1p(rna_feat)

    log.info("[%s] Fitting cosine KNN (k=%d) on RNA (%d cells) ...",
             compartment, n_anchors, rna_feat_norm.shape[0])
    nn = NearestNeighbors(n_neighbors=n_anchors, metric="cosine", n_jobs=-1)
    nn.fit(rna_feat_norm)
    log.info("[%s] Querying KNN for %d ATAC cells ...",
             compartment, atac_feat_norm.shape[0])
    _, idx = nn.kneighbors(atac_feat_norm)

    rna_ct = rna.obs["cell_type"].astype(str).values
    atac_ct = atac.obs["cell_type"].astype(str).values

    n_atac = len(atac_ct)
    agreement = np.zeros(n_atac, dtype=np.float32)
    top_anchor_ct = np.empty(n_atac, dtype=object)
    chunk = 5_000
    for s in range(0, n_atac, chunk):
        e = min(s + chunk, n_atac)
        anchor_ct_matrix = rna_ct[idx[s:e]]
        atac_slice = atac_ct[s:e]
        agreement_chunk = np.mean(
            anchor_ct_matrix == atac_slice[:, None], axis=1
        )
        agreement[s:e] = agreement_chunk.astype(np.float32)
        top_anchor_ct[s:e] = anchor_ct_matrix[:, 0]

    result = pd.DataFrame({
        "barcode":         atac.obs_names,
        "cell_type_atac":  atac_ct,
        "top_anchor_ct":   top_anchor_ct,
        "anchor_agreement": agreement,
        "compartment":     compartment,
    })

    q = result["anchor_agreement"].quantile([0.10, 0.25, 0.50, 0.75, 0.90]).values
    log.info(
        "[%s] Anchor-agreement: p10=%.3f p25=%.3f p50=%.3f p75=%.3f p90=%.3f  (mean=%.3f)",
        compartment, q[0], q[1], q[2], q[3], q[4],
        float(result["anchor_agreement"].mean()),
    )
    per_ct = (result.groupby("cell_type_atac")
                     .agg(n=("anchor_agreement", "size"),
                          mean_agreement=("anchor_agreement", "mean"),
                          median_agreement=("anchor_agreement", "median"))
                     .sort_values("n", ascending=False))
    log.info("[%s] Per cell-type agreement (top 10):\n%s",
             compartment, per_ct.head(10).to_string())

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_tsv = OUT_DIR / f"{compartment}.tsv"
    result.to_csv(out_tsv, sep="\t", index=False)
    per_ct_tsv = OUT_DIR / f"{compartment}_per_cell_type.tsv"
    per_ct.to_csv(per_ct_tsv, sep="\t")
    log.info("[%s] Saved:\n  %s\n  %s", compartment, out_tsv, per_ct_tsv)

    log.info("[%s] Done (%.1f s)", compartment, time.time() - t0)
    return result


def main(argv: Optional[list[str]] = None) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  %(levelname)-8s  %(message)s",
                        datefmt="%H:%M:%S")
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--compartment", choices=list(DEFAULT_ATAC_H5AD.keys()))
    p.add_argument("--all", action="store_true")
    args = p.parse_args(argv)

    if args.all:
        for c in DEFAULT_ATAC_H5AD:
            diagnose_compartment(c)
    else:
        if not args.compartment:
            sys.exit("--compartment required unless --all")
        diagnose_compartment(args.compartment)


if __name__ == "__main__":
    main()
