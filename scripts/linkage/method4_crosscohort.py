"""
Method 4B - Cross-cohort peak-gene correlation via anchor transfer.

  Filename note: file is `method4_crosscohort.py` (no 'b') - canonical
  method label `method4b_crosscohort` is used in figures, tables, config,
  and `build_annotations.py:METHOD_SOURCES`. The diagnostic for THIS
  method's failure mode lives in `method4b_anchor_diagnostic.py` (which
  is *not* a linkage method itself - see that file's header).


Unpaired analog of Method 4A: imputes RNA for every ATAC cell by anchor
transfer from the RNA-only cohort in a shared gene-activity feature space,
then runs the same per-gene meta-cell Pearson correlation as 4A within
+/- window_bp of each TSS.

Usage
-----
    python scripts/linkage/method4_crosscohort.py \\
        --rna-h5ad data/processed/hickey2023/rna_qc.h5ad \\
        --atac-h5ad data/processed/hickey2023/atac_immune.h5ad \\
        --compartment immune \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --out results/linkage/method4_crosscohort/immune.tsv

    python scripts/linkage/method4_crosscohort.py --all-compartments \\
        --rna-h5ad data/processed/hickey2023/rna_qc.h5ad \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --out-dir results/linkage/method4_crosscohort/
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
import pybedtools
import scanpy as sc
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

from method4_paired import (
    CANONICAL_CHROMS,
    DEFAULT_ATAC_H5AD,
    aggregate_metacells,
    build_metacell_indicator,
    extract_tss_frame,
    lsi_embed,
    peak_gene_correlation,
    standardise_cols,
    _peaks_from_h5ad,
    _summarise,
)

log = logging.getLogger(__name__)


def compute_gene_activity(
    atac,
    tss_df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    *,
    upstream_bp: int = 2_000,
    downstream_bp: int = 0,
    gene_body_bp: int = 100_000,
) -> tuple[np.ndarray, list[str]]:
    """Build a (n_atac_cells, n_genes) gene-activity matrix.

    For each gene, sums ATAC signal across peaks overlapping
    [TSS - upstream_bp, TSS + gene_body_bp + downstream_bp] on the gene's strand.
    """
    t0 = time.time()
    rows = []
    for _, g in tss_df.iterrows():
        if g.strand == "+":
            start = max(g.tss - upstream_bp, 0)
            end   = g.tss + gene_body_bp + downstream_bp
        else:
            start = max(g.tss - gene_body_bp - downstream_bp, 0)
            end   = g.tss + upstream_bp
        rows.append((g.chrom, int(start), int(end), g.gene_id))
    gene_df = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene_id"])

    gene_bt  = pybedtools.BedTool.from_dataframe(gene_df).sort()
    peaks_bt = pybedtools.BedTool.from_dataframe(
        peaks_df[["chrom", "start", "end", "peak_id"]]
    ).sort()

    log.info("gene-activity: intersecting %d peaks with %d gene-body windows",
             len(peaks_df), len(gene_df))
    hits = gene_bt.intersect(peaks_bt, wa=True, wb=True)

    peak_idx_map = {pid: i for i, pid in enumerate(peaks_df["peak_id"])}
    gene_to_peaks: dict[str, list[int]] = {}
    for h in hits:
        gid    = h[3]
        p_id   = h[7]
        pcol   = peak_idx_map.get(p_id)
        if pcol is None:
            continue
        gene_to_peaks.setdefault(gid, []).append(pcol)

    gene_order = [g for g in tss_df["gene_id"] if g in gene_to_peaks]
    log.info("gene-activity: %d of %d genes have at least one overlapping peak",
             len(gene_order), len(tss_df))

    n_peaks = atac.n_vars
    rows_i, cols_i = [], []
    for gi, g in enumerate(gene_order):
        for p in gene_to_peaks[g]:
            rows_i.append(gi); cols_i.append(p)
    W = sp.csr_matrix(
        (np.ones(len(rows_i), dtype=np.float32), (rows_i, cols_i)),
        shape=(len(gene_order), n_peaks),
    )

    X = atac.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X_bin = (X > 0).astype(np.float32)
    activity = (X_bin @ W.T).toarray()
    log.info("gene-activity: matrix %s built in %.1fs",
             activity.shape, time.time() - t0)
    return activity.astype(np.float32), gene_order


def _normalise_log1p(X: np.ndarray, target_sum: float = 10_000) -> np.ndarray:
    """Per-row: scale to target_sum, then log1p. Idempotent for already-log1p data."""
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    max_abs = float(np.abs(X).max()) if X.size else 0.0
    if max_abs < 50.0:
        log.debug("_normalise_log1p: max |X|=%.2f -> passing through", max_abs)
        return X
    rs = X.sum(axis=1, keepdims=True) + 1e-6
    out = np.log1p(X * (target_sum / rs))
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)


def anchor_knn_transfer(
    atac_features: np.ndarray,
    rna_features:  np.ndarray,
    rna_counts_raw:  np.ndarray,
    *,
    n_anchors: int = 30,
    metric: str = "cosine",
) -> np.ndarray:
    """Impute RNA per ATAC cell as Gaussian-weighted avg of nearest RNA cells."""
    t0 = time.time()
    log.info("anchor KNN: fitting on RNA (%s features, %d RNA cells) ...",
             metric, rna_features.shape[0])
    nn = NearestNeighbors(n_neighbors=n_anchors, metric=metric, n_jobs=-1)
    nn.fit(rna_features)

    log.info("anchor KNN: querying %d ATAC cells", atac_features.shape[0])
    dists, idx = nn.kneighbors(atac_features)

    # Gaussian kernel: sigma = row-wise median distance (robust)
    sigma = np.maximum(np.median(dists, axis=1, keepdims=True), 1e-6)
    weights = np.exp(-(dists ** 2) / (2 * sigma ** 2))
    weights /= weights.sum(axis=1, keepdims=True)

    n_atac, _ = dists.shape
    n_genes_rna = rna_counts_raw.shape[1]
    imputed = np.zeros((n_atac, n_genes_rna), dtype=np.float32)
    chunk = 10_000
    for s in range(0, n_atac, chunk):
        e = min(s + chunk, n_atac)
        neigh_counts = rna_counts_raw[idx[s:e]]
        w = weights[s:e][:, :, None]
        imputed[s:e] = (neigh_counts * w).sum(axis=1)
    log.info("anchor KNN: imputed RNA shape=%s in %.1fs",
             imputed.shape, time.time() - t0)
    return imputed


def run_compartment_crosscohort(
    rna_h5ad: str,
    atac_h5ad: str,
    compartment: str,
    gtf_path: str,
    out_tsv: str,
    *,
    window_bp: int = 500_000,
    # Use the same ArchR-default threshold as Method 4A for honest comparison.
    r_threshold: float = 0.45,
    fdr_threshold: float = 0.05,
    n_metacells: int = 2_000,
    k_knn: int = 50,
    min_meta_expr: int = 5,
    n_anchors: int = 30,
    upstream_bp: int = 2_000,
    gene_body_bp: int = 100_000,
    random_state: int = 42,
) -> pd.DataFrame:
    """Full Method 4B pipeline for one compartment."""
    t0 = time.time()
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)

    log.info("[%s] Loading RNA reference = %s", compartment, rna_h5ad)
    log.info("[%s] Loading ATAC target   = %s", compartment, atac_h5ad)
    rna  = sc.read_h5ad(rna_h5ad)
    atac = sc.read_h5ad(atac_h5ad)

    tss_df = extract_tss_frame(gtf_path)

    # If RNA uses gene-symbol var_names, rekey tss_df so gene_id == gene_name
    use_gene_name_key = "gene_ids" not in rna.var.columns
    if use_gene_name_key:
        log.warning(
            "[%s] rna.var has no `gene_ids` column - using gene_name as bridge key",
            compartment,
        )
        tss_df = tss_df.assign(gene_id=tss_df["gene_name"])
        tss_df = tss_df.drop_duplicates(subset=["gene_id"]).reset_index(drop=True)

    peaks_df = _peaks_from_h5ad(atac)
    keep = peaks_df["chrom"].isin(CANONICAL_CHROMS).values
    atac = atac[:, keep].copy()
    peaks_df = peaks_df.loc[keep].reset_index(drop=True)
    log.info("[%s] ATAC: %d cells x %d peaks (canonical)",
             compartment, atac.n_obs, atac.n_vars)

    activity, gene_order = compute_gene_activity(
        atac, tss_df, peaks_df,
        upstream_bp=upstream_bp,
        gene_body_bp=gene_body_bp,
    )
    gene_activity_idx = {g: i for i, g in enumerate(gene_order)}

    if "gene_ids" in rna.var.columns:
        rna_gene_ids = rna.var["gene_ids"].astype(str).values
    else:
        log.warning("rna.var has no `gene_ids` column - falling back to gene_name")
        rna_gene_ids = rna.var_names.astype(str).values
    rna_gene_idx = {g: i for i, g in enumerate(rna_gene_ids)}

    common = [g for g in gene_order if g in rna_gene_idx]
    log.info("[%s] bridge feature space: %d common genes (of %d activity / %d RNA)",
             compartment, len(common), len(gene_order), len(rna_gene_ids))
    if len(common) < 500:
        log.warning("[%s] fewer than 500 common genes - integration may be weak.",
                    compartment)

    act_cols = [gene_activity_idx[g] for g in common]
    rna_cols = [rna_gene_idx[g] for g in common]

    atac_feat = activity[:, act_cols]
    rna_raw = rna.X
    if sp.issparse(rna_raw):
        rna_raw = rna_raw.toarray()
    rna_feat = rna_raw[:, rna_cols]

    atac_feat_norm = _normalise_log1p(atac_feat)
    rna_feat_norm  = _normalise_log1p(rna_feat)

    rna_full = rna.X
    if sp.issparse(rna_full):
        rna_full = rna_full.toarray()
    imputed_rna = anchor_knn_transfer(
        atac_features=atac_feat_norm,
        rna_features=rna_feat_norm,
        rna_counts_raw=rna_full.astype(np.float32),
        n_anchors=n_anchors,
    )

    X_atac = atac.X
    if not sp.issparse(X_atac):
        X_atac = sp.csr_matrix(X_atac)
    X_bin = (X_atac > 0).astype(np.float32)
    lsi = lsi_embed(X_bin, random_state=random_state)
    ind = build_metacell_indicator(
        lsi, n_metacells=n_metacells, k_knn=k_knn, random_state=random_state,
    )

    atac_meta = aggregate_metacells(X_bin, ind)
    atac_meta_std = standardise_cols(atac_meta)
    del atac_meta

    rna_meta_raw = aggregate_metacells(sp.csr_matrix(imputed_rna), ind)
    rna_meta_std = standardise_cols(np.log1p(rna_meta_raw))

    gene_idx_map = rna_gene_idx
    links = peak_gene_correlation(
        atac_meta_std=atac_meta_std,
        rna_meta_std=rna_meta_std,
        rna_meta_raw=rna_meta_raw,
        peaks_df=peaks_df,
        gene_idx_map=gene_idx_map,
        tss_df=tss_df,
        compartment=compartment,
        window_bp=window_bp,
        r_threshold=r_threshold,
        fdr_threshold=fdr_threshold,
        min_meta_expr=min_meta_expr,
    )
    links["method"] = "method4_crosscohort"

    links.to_csv(out_tsv, sep="\t", index=False)
    log.info("[%s] %d peak-gene links -> %s  (total %.1f s)",
             compartment, len(links), out_tsv, time.time() - t0)
    _summarise(links, compartment=compartment)
    return links


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--rna-h5ad", required=True, metavar="H5AD",
                   help="RNA reference h5ad (e.g. rna_qc.h5ad from B001/B004/B005).")
    p.add_argument("--atac-h5ad", metavar="H5AD",
                   help="Per-compartment ATAC h5ad (required unless --all-compartments).")
    p.add_argument("--compartment", default="all",
                   help='Compartment label (default "all").')
    p.add_argument("--out", metavar="TSV",
                   help="Output TSV path (required unless --all-compartments).")

    p.add_argument("--gtf", required=True,
                   help="GENCODE annotation GTF (plain or .gz).")

    p.add_argument("--window-bp", type=int, default=500_000, dest="window_bp",
                   help="Half-window (default 500000 = +/-500 kb).")
    p.add_argument("--r-threshold", type=float, default=0.45, dest="r_threshold",
                   help="|r| threshold (default 0.45).")
    p.add_argument("--fdr-threshold", type=float, default=0.05, dest="fdr_threshold",
                   help="Per-gene FDR threshold (default 0.05).")
    p.add_argument("--n-metacells", type=int, default=2_000, dest="n_metacells")
    p.add_argument("--k-knn", type=int, default=50, dest="k_knn")
    p.add_argument("--min-meta-expr", type=int, default=5, dest="min_meta_expr")
    p.add_argument("--n-anchors", type=int, default=30, dest="n_anchors",
                   help="K nearest RNA cells used as anchors per ATAC cell (default 30).")
    p.add_argument("--upstream-bp", type=int, default=2_000, dest="upstream_bp",
                   help="Gene-body window upstream of TSS (default 2000).")
    p.add_argument("--gene-body-bp", type=int, default=100_000, dest="gene_body_bp",
                   help="Gene-body window downstream of TSS (default 100000).")
    p.add_argument("--seed", type=int, default=42, dest="random_state")

    p.add_argument("--all-compartments", action="store_true", dest="all_compartments",
                   help="Run epithelial + immune + stromal using defaults.")
    p.add_argument("--out-dir", dest="out_dir",
                   help="Output directory with --all-compartments.")
    return p


def main(argv: Optional[list[str]] = None) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )
    args = _build_parser().parse_args(argv)

    if args.all_compartments:
        if not args.out_dir:
            sys.exit("--all-compartments requires --out-dir")
        out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        for comp, h5ad in DEFAULT_ATAC_H5AD.items():
            if not Path(h5ad).exists():
                log.warning("ATAC h5ad missing, skipping %s: %s", comp, h5ad)
                continue
            run_compartment_crosscohort(
                rna_h5ad=args.rna_h5ad,
                atac_h5ad=h5ad,
                compartment=comp,
                gtf_path=args.gtf,
                out_tsv=str(out_dir / f"{comp}.tsv"),
                window_bp=args.window_bp,
                r_threshold=args.r_threshold,
                fdr_threshold=args.fdr_threshold,
                n_metacells=args.n_metacells,
                k_knn=args.k_knn,
                min_meta_expr=args.min_meta_expr,
                n_anchors=args.n_anchors,
                upstream_bp=args.upstream_bp,
                gene_body_bp=args.gene_body_bp,
                random_state=args.random_state,
            )
    else:
        if not args.atac_h5ad or not args.out:
            sys.exit("--atac-h5ad and --out are required without --all-compartments")
        run_compartment_crosscohort(
            rna_h5ad=args.rna_h5ad,
            atac_h5ad=args.atac_h5ad,
            compartment=args.compartment,
            gtf_path=args.gtf,
            out_tsv=args.out,
            window_bp=args.window_bp,
            r_threshold=args.r_threshold,
            fdr_threshold=args.fdr_threshold,
            n_metacells=args.n_metacells,
            k_knn=args.k_knn,
            min_meta_expr=args.min_meta_expr,
            n_anchors=args.n_anchors,
            upstream_bp=args.upstream_bp,
            gene_body_bp=args.gene_body_bp,
            random_state=args.random_state,
        )


if __name__ == "__main__":
    main()
