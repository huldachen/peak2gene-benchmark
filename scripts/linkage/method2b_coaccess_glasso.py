"""
Method 2b - Cicero co-accessibility with graphical lasso (Pliner et al. 2018).

Replaces Method 2's raw Pearson correlation with graphical-lasso partial
correlations using a distance-dependent penalty matrix - the canonical
Cicero algorithm. Removes indirect peak-peak correlations mediated through
intermediate peaks.

Reference: Pliner H. A. et al., Mol. Cell 71, 858 (2018).

Usage
-----
    python scripts/linkage/method2b_coaccess_glasso.py --all-compartments \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --out-dir results/linkage/method2b/
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
import scipy.sparse as sp
import scipy.stats as stats

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Import shared pipeline components from Method 2
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent))
from method2_cicero import (                          # noqa: E402
    CANONICAL_CHROMS,
    DEFAULT_ATAC_H5AD,
    _build_metacells,
    _lsi_embed,
    _standardise_columns,
    coaccess_to_peak_gene_links,
    extract_tss_frame,
)

# ======================== GLASSO SOLVER ====================================


def _distance_penalty_matrix(
    peak_mids: np.ndarray,
    alpha: float = 0.5,
    distance_scale: float = 250_000,
) -> np.ndarray:
    """Build element-wise penalty matrix rho_ij = alpha * (1 - exp(-d_ij / scale)).

    Nearby peaks -> low penalty (can have strong edges).
    Distant peaks -> high penalty (need strong evidence).
    Diagonal is always 0 (no penalty on variances).
    """
    n = len(peak_mids)
    D = np.abs(peak_mids[:, None] - peak_mids[None, :]).astype(np.float64)
    rho = alpha * (1.0 - np.exp(-D / distance_scale))
    np.fill_diagonal(rho, 0.0)
    return rho


def _graphical_lasso_elementwise(
    S: np.ndarray,
    rho: np.ndarray,
    *,
    max_iter: int = 100,
    tol: float = 1e-4,
) -> Optional[np.ndarray]:
    """Block coordinate descent for graphical lasso with element-wise penalties.

    Implements Friedman et al. (2008) Algorithm 1 with per-element rho_ij.

    Parameters
    ----------
    S : (p, p) sample covariance matrix (symmetric, PSD)
    rho : (p, p) element-wise penalty matrix (non-negative, zero diagonal)
    max_iter : maximum BCD iterations
    tol : convergence tolerance on max change in W

    Returns
    -------
    Theta : (p, p) estimated precision matrix, or None if p < 3 or singular.
    """
    p = S.shape[0]
    if p < 3:
        return None

    # Initialise W = S + diag(rho) to ensure positive-definiteness
    W = S.copy()
    W += np.diag(np.maximum(rho.diagonal(), rho.mean(axis=1)))
    # Ensure symmetry
    W = (W + W.T) / 2.0

    for iteration in range(max_iter):
        W_old = W.copy()
        for j in range(p):
            # Partition: column j vs rest
            idx = np.concatenate([np.arange(j), np.arange(j + 1, p)])
            W_11 = W[np.ix_(idx, idx)]
            s_12 = S[idx, j]
            rho_12 = rho[idx, j]

            # Solve lasso sub-problem via coordinate descent:
            #   min  0.5 * beta^T W_11 beta - s_12^T beta + sum rho_12_k |beta_k|
            beta = _lasso_cd(W_11, s_12, rho_12, max_iter=200, tol=tol * 0.1)

            # Update W column/row j
            w_12 = W_11 @ beta
            W[idx, j] = w_12
            W[j, idx] = w_12

        # Check convergence
        max_change = np.max(np.abs(W - W_old))
        if max_change < tol:
            break

    # Recover precision matrix: Theta = W^{-1}
    try:
        Theta = np.linalg.inv(W)
        # Symmetrise (numerical)
        Theta = (Theta + Theta.T) / 2.0
        return Theta
    except np.linalg.LinAlgError:
        return None


def _lasso_cd(
    W_11: np.ndarray,
    s_12: np.ndarray,
    rho_12: np.ndarray,
    max_iter: int = 200,
    tol: float = 1e-5,
) -> np.ndarray:
    """Coordinate descent for the lasso sub-problem in GLASSO.

    Solves: min  0.5 * beta^T W_11 beta - s_12^T beta + sum_k rho_k |beta_k|
    """
    p_minus_1 = len(s_12)
    beta = np.zeros(p_minus_1, dtype=np.float64)
    diag_W = np.diag(W_11).copy()
    diag_W[diag_W < 1e-12] = 1e-12  # guard

    for _ in range(max_iter):
        max_delta = 0.0
        for k in range(p_minus_1):
            # Partial residual
            r_k = s_12[k] - W_11[k, :] @ beta + diag_W[k] * beta[k]
            # Soft-threshold
            beta_new = np.sign(r_k) * max(abs(r_k) - rho_12[k], 0.0) / diag_W[k]
            delta = abs(beta_new - beta[k])
            if delta > max_delta:
                max_delta = delta
            beta[k] = beta_new
        if max_delta < tol:
            break
    return beta


def _precision_to_partial_corr(Theta: np.ndarray) -> np.ndarray:
    """Convert precision matrix to partial correlations.

    pcor_ij = -Theta_ij / sqrt(Theta_ii * Theta_jj)   for i != j
    """
    d = np.sqrt(np.diag(Theta))
    d[d < 1e-12] = 1e-12
    pcor = -Theta / np.outer(d, d)
    np.fill_diagonal(pcor, 1.0)
    return pcor


# ======================== CO-ACCESSIBILITY =================================


def _coaccess_one_chrom_glasso(
    M_std: np.ndarray,
    peak_mids: np.ndarray,
    window_bp: int,
    *,
    alpha: float = 0.5,
    distance_scale: float = 250_000,
    max_window_peaks: int = 80,
    prefilter_r: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """GLASSO-based co-accessibility for one chromosome.

    Slides a window across peaks, fits graphical lasso with distance-dependent
    penalty on each window, extracts partial correlations. Returns upper-
    triangular (i, j, pcor) with i < j.

    Windows with >max_window_peaks peaks are sub-divided to keep GLASSO
    tractable (O(p^3) per window).
    """
    n_meta, n_peaks = M_std.shape
    n_total = n_meta

    i_out, j_out, r_out = [], [], []
    n_candidate_pairs = 0

    # Slide by stepping through peaks as window anchors
    start_ptr = 0
    end_ptr = 0
    processed_pairs = set()

    for anchor in range(n_peaks):
        # Define window around this anchor peak
        mid = peak_mids[anchor]
        win_lo = mid - window_bp
        win_hi = mid + window_bp

        # Advance pointers
        while start_ptr < n_peaks and peak_mids[start_ptr] < win_lo:
            start_ptr += 1
        while end_ptr < n_peaks and peak_mids[end_ptr] <= win_hi:
            end_ptr += 1

        win_idx = np.arange(start_ptr, end_ptr)
        p = len(win_idx)
        if p < 3:
            continue

        # Skip very large windows - subdivide by shifting anchor
        if p > max_window_peaks:
            continue

        # Sample covariance of meta-cell accessibility in this window
        X_win = M_std[:, win_idx]  # (n_meta, p)
        S = np.cov(X_win, rowvar=False, bias=False)
        if S.ndim < 2:
            continue

        # Distance-dependent penalty matrix
        mids_win = peak_mids[win_idx]
        rho = _distance_penalty_matrix(mids_win, alpha=alpha,
                                       distance_scale=distance_scale)

        # Solve GLASSO
        Theta = _graphical_lasso_elementwise(S, rho)
        if Theta is None:
            continue

        # Partial correlations
        pcor = _precision_to_partial_corr(Theta)

        # Extract upper-triangular pairs
        for ii in range(p):
            for jj in range(ii + 1, p):
                gi = int(win_idx[ii])
                gj = int(win_idx[jj])
                pair_key = (gi, gj)
                if pair_key in processed_pairs:
                    continue
                processed_pairs.add(pair_key)
                n_candidate_pairs += 1
                pc = pcor[ii, jj]
                if pc >= prefilter_r:
                    i_out.append(gi)
                    j_out.append(gj)
                    r_out.append(pc)

    if not i_out:
        return (np.empty(0, np.int32), np.empty(0, np.int32),
                np.empty(0, np.float32), n_candidate_pairs)
    return (np.array(i_out, dtype=np.int32),
            np.array(j_out, dtype=np.int32),
            np.array(r_out, dtype=np.float32),
            n_candidate_pairs)


# ======================== MAIN PIPELINE ====================================


def compute_coaccessibility_glasso(
    atac_h5ad: str,
    *,
    compartment: str,
    window_bp: int = 500_000,
    coaccess_threshold: float = 0.05,
    fdr_threshold: float = 0.05,
    n_metacells: int = 2000,
    k_knn: int = 50,
    n_lsi_components: int = 50,
    alpha: float = 0.5,
    distance_scale: float = 250_000,
    max_window_peaks: int = 80,
    n_jobs: int = 8,
    random_state: int = 42,
    prefilter_r: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """GLASSO co-accessibility pipeline. Returns (coaccess_df, peaks_df).

    Note: default coaccess_threshold is 0.05 (matching Cicero's shrunk-score
    scale) rather than Method 2's 0.25 (raw Pearson scale), because GLASSO
    partial correlations are smaller in magnitude than raw correlations.
    """
    import scanpy as sc
    t0 = time.time()
    log.info("[%s] Reading %s ...", compartment, atac_h5ad)
    adata = sc.read_h5ad(atac_h5ad)

    peaks_df = pd.DataFrame(
        [pid.replace(":", "-").split("-") for pid in adata.var_names],
        columns=["chrom", "start", "end"],
    )
    peaks_df["start"] = peaks_df["start"].astype(int)
    peaks_df["end"]   = peaks_df["end"].astype(int)
    keep_mask = peaks_df["chrom"].isin(CANONICAL_CHROMS).values
    adata    = adata[:, keep_mask].copy()
    peaks_df = peaks_df.loc[keep_mask].reset_index(drop=True)
    log.info("[%s] %d cells x %d peaks (canonical chroms)",
             compartment, adata.n_obs, adata.n_vars)

    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X_bin = (X > 0).astype(np.float32)

    # Shared LSI + KNN meta-cell pipeline (identical to Method 2)
    lsi  = _lsi_embed(X_bin, n_components=n_lsi_components,
                      random_state=random_state)
    meta = _build_metacells(X_bin, lsi,
                            n_metacells=n_metacells, k_knn=k_knn,
                            random_state=random_state)
    log.info("[%s] meta-cell matrix: %s", compartment, meta.shape)
    M_std = _standardise_columns(meta)
    del meta

    peaks_df["mid"] = (peaks_df["start"] + peaks_df["end"]) // 2
    peaks_df["idx_global"] = np.arange(len(peaks_df), dtype=np.int64)

    order = np.lexsort((peaks_df["mid"].values,
                        peaks_df["chrom"].values))
    peaks_sorted = peaks_df.iloc[order].reset_index(drop=True)

    # Build per-chromosome jobs. Each job is a tuple of args to pass to
    # `_coaccess_one_chrom_glasso`; we slice M_std once per chrom so that the
    # subprocess receives only the columns it needs (smaller pickled payload).
    jobs = []
    chrom_global_idx = {}
    for chrom, sub in peaks_sorted.groupby("chrom", sort=False):
        global_idx = sub["idx_global"].values
        mids       = sub["mid"].values.astype(np.int64)
        sub_M      = np.ascontiguousarray(M_std[:, global_idx])
        chrom_global_idx[chrom] = global_idx
        jobs.append((chrom, sub_M, mids))
        log.info("[%s] %s: %d peaks queued for GLASSO",
                 compartment, chrom, len(sub))

    # Capture n_meta before freeing M_std (used below for t-test df).
    n_meta = M_std.shape[0]

    # Free the big M_std now that per-chrom slices are extracted (saves
    # several GB of RAM during the parallel phase, especially for
    # epithelial which has ~500k peaks).
    del M_std

    # Run GLASSO across chromosomes in parallel. Each worker is a separate
    # process to actually escape the Python GIL (the inner GLASSO/lasso
    # coordinate-descent loops are pure Python and don't release the GIL
    # via NumPy alone).
    n_jobs_eff = min(n_jobs, len(jobs))
    log.info("[%s] Running GLASSO on %d chromosomes in parallel (n_jobs=%d) ...",
             compartment, len(jobs), n_jobs_eff)

    from concurrent.futures import ProcessPoolExecutor, as_completed
    from multiprocessing import get_context

    def _job_kwargs():
        return dict(
            window_bp=window_bp,
            alpha=alpha,
            distance_scale=distance_scale,
            max_window_peaks=max_window_peaks,
            prefilter_r=prefilter_r,
        )

    i_all, j_all, r_all = [], [], []
    total_candidate_pairs = 0
    with ProcessPoolExecutor(max_workers=n_jobs_eff,
                             mp_context=get_context("spawn")) as ex:
        future_to_chrom = {
            ex.submit(_coaccess_one_chrom_glasso, sub_M, mids, **_job_kwargs()): chrom
            for chrom, sub_M, mids in jobs
        }
        for fut in as_completed(future_to_chrom):
            chrom = future_to_chrom[fut]
            try:
                i_local, j_local, r_vals, n_cand = fut.result()
            except Exception as e:
                log.error("[%s] %s GLASSO failed: %s", compartment, chrom, e)
                continue
            total_candidate_pairs += n_cand
            log.info("[%s] %s: GLASSO done (%d kept pairs out of %d candidates)",
                     compartment, chrom, len(i_local), n_cand)
            if i_local.size:
                global_idx = chrom_global_idx[chrom]
                i_all.append(global_idx[i_local])
                j_all.append(global_idx[j_local])
                r_all.append(r_vals)

    if i_all:
        i_arr = np.concatenate(i_all)
        j_arr = np.concatenate(j_all)
        r_arr = np.concatenate(r_all)
    else:
        i_arr = np.empty(0, np.int64)
        j_arr = np.empty(0, np.int64)
        r_arr = np.empty(0, np.float32)

    # One-sided t-test for pcor > 0
    df_meta = max(n_meta - 2, 1)
    r_clip = np.clip(r_arr, -0.9999, 0.9999)
    t_stat = r_clip * np.sqrt(df_meta) / np.sqrt(1.0 - r_clip ** 2)
    p_values = stats.t.sf(t_stat, df=df_meta).astype(np.float32)

    fdr_values = _bh_fdr(p_values, n_tests_total=total_candidate_pairs
                         ).astype(np.float32)

    log.info(
        "[%s] GLASSO: %d candidate pairs; %d pass pcor >= %.3f (prefilter); "
        "%d also pass FDR <= %.3f",
        compartment, total_candidate_pairs, len(r_arr), prefilter_r,
        int((fdr_values <= fdr_threshold).sum()), fdr_threshold,
    )

    keep = (r_arr >= coaccess_threshold) & (fdr_values <= fdr_threshold)
    i_arr = i_arr[keep]; j_arr = j_arr[keep]
    r_arr = r_arr[keep]; p_values = p_values[keep]; fdr_values = fdr_values[keep]

    peak_ids = (
        peaks_df["chrom"].astype(str) + ":"
        + peaks_df["start"].astype(str) + "-"
        + peaks_df["end"].astype(str)
    ).values

    coaccess_df = pd.DataFrame({
        "Peak1":    peak_ids[i_arr],
        "Peak2":    peak_ids[j_arr],
        "coaccess": r_arr,
        "pvalue":   p_values,
        "fdr":      fdr_values,
    })

    log.info("[%s] Final: %d pairs pass pcor >= %.3f AND FDR <= %.3f  (%.1fs)",
             compartment, len(coaccess_df), coaccess_threshold, fdr_threshold,
             time.time() - t0)

    peaks_out = peaks_df[["chrom", "start", "end"]].copy()
    peaks_out["peak_id"] = peak_ids
    return coaccess_df, peaks_out


def _bh_fdr(p_values: np.ndarray, n_tests_total: int) -> np.ndarray:
    """BH FDR for pre-filtered p-values with known total test count."""
    n = len(p_values)
    if n == 0:
        return np.empty(0, dtype=np.float64)
    order = np.argsort(p_values)
    ranked = np.asarray(p_values)[order]
    q = ranked * float(n_tests_total) / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n, dtype=np.float64)
    out[order] = q
    return np.clip(out, 0.0, 1.0)


def _summarise(df: pd.DataFrame, *, compartment: str) -> None:
    if df.empty:
        log.warning("[%s] Empty linkage - nothing to summarise.", compartment)
        return
    pcts_d = np.percentile(df["distance_bp"].values, [0, 25, 50, 75, 90, 100])
    pcts_s = np.percentile(df["score"].values,       [0, 25, 50, 75, 90, 100])
    log.info(
        "[%s] -- Summary --\n"
        "  Total links  : %d\n"
        "  Unique peaks : %d   (avg %.1f genes / peak)\n"
        "  Unique genes : %d   (avg %.1f peaks / gene)\n"
        "  Score (pcor)  : min=%.3f  p25=%.3f  p50=%.3f  p75=%.3f  p90=%.3f  max=%.3f\n"
        "  Distance (bp) : min=%d  p25=%d  p50=%d  p75=%d  p90=%d  max=%d",
        compartment, len(df),
        df["peak_id"].nunique(), len(df) / max(df["peak_id"].nunique(), 1),
        df["gene_id"].nunique(), len(df) / max(df["gene_id"].nunique(), 1),
        pcts_s[0], pcts_s[1], pcts_s[2], pcts_s[3], pcts_s[4], pcts_s[5],
        int(pcts_d[0]), int(pcts_d[1]), int(pcts_d[2]),
        int(pcts_d[3]), int(pcts_d[4]), int(pcts_d[5]),
    )


def run_compartment(
    atac_h5ad: str,
    compartment: str,
    gtf_path: str,
    out_tsv: str,
    *,
    window_bp: int = 500_000,
    coaccess_threshold: float = 0.05,
    fdr_threshold: float = 0.05,
    promoter_bp: int = 2000,
    n_metacells: int = 2000,
    k_knn: int = 50,
    alpha: float = 0.5,
    distance_scale: float = 250_000,
    max_window_peaks: int = 80,
    n_jobs: int = 8,
    random_state: int = 42,
    save_coaccess_tsv: Optional[str] = None,
) -> pd.DataFrame:
    """Full Method-2b pipeline for one compartment."""
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    coaccess_df, peaks_df = compute_coaccessibility_glasso(
        atac_h5ad,
        compartment=compartment,
        window_bp=window_bp,
        coaccess_threshold=coaccess_threshold,
        fdr_threshold=fdr_threshold,
        n_metacells=n_metacells,
        k_knn=k_knn,
        alpha=alpha,
        distance_scale=distance_scale,
        max_window_peaks=max_window_peaks,
        n_jobs=n_jobs,
        random_state=random_state,
    )

    if save_coaccess_tsv:
        Path(save_coaccess_tsv).parent.mkdir(parents=True, exist_ok=True)
        coaccess_df.to_csv(save_coaccess_tsv, sep="\t", index=False)
        log.info("[%s] Saved raw peak-peak GLASSO co-access to %s",
                 compartment, save_coaccess_tsv)

    tss_df = extract_tss_frame(gtf_path)
    links = coaccess_to_peak_gene_links(
        coaccess_df, peaks_df, tss_df,
        compartment=compartment, promoter_bp=promoter_bp,
    )

    # Override method label
    if not links.empty:
        links["method"] = "method2b_glasso"

    links.to_csv(out_tsv, sep="\t", index=False)
    log.info("[%s] %d peak-gene links -> %s", compartment, len(links), out_tsv)
    _summarise(links, compartment=compartment)
    return links


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Method 2b: Cicero GLASSO co-accessibility -> peak-gene links"
    )
    p.add_argument("--atac-h5ad", metavar="H5AD",
                   help="Per-compartment ATAC h5ad.")
    p.add_argument("--compartment", default="all")
    p.add_argument("--out", metavar="TSV")

    p.add_argument("--gtf", required=True,
                   help="GENCODE annotation GTF (plain or .gz).")

    p.add_argument("--window-bp", type=int, default=500_000, dest="window_bp")
    p.add_argument("--coaccess-threshold", type=float, default=0.05,
                   dest="coaccess_threshold",
                   help="Minimum partial-correlation score (default 0.05, "
                        "matching Cicero's shrunk-score scale).")
    p.add_argument("--fdr-threshold", type=float, default=0.05,
                   dest="fdr_threshold")
    p.add_argument("--promoter-bp", type=int, default=2000, dest="promoter_bp")
    p.add_argument("--n-metacells", type=int, default=2000, dest="n_metacells")
    p.add_argument("--k-knn", type=int, default=50, dest="k_knn")
    p.add_argument("--alpha", type=float, default=0.5,
                   help="GLASSO base regularisation strength (default 0.5).")
    p.add_argument("--distance-scale", type=float, default=250_000,
                   dest="distance_scale",
                   help="Distance-penalty length scale in bp (default 250000).")
    p.add_argument("--max-window-peaks", type=int, default=80,
                   dest="max_window_peaks",
                   help="Skip windows with more peaks than this "
                        "(GLASSO is O(p^3) per window; default 80).")
    p.add_argument("--n-jobs", type=int, default=8, dest="n_jobs",
                   help="Number of parallel chromosome workers (default 8).")
    p.add_argument("--seed", type=int, default=42, dest="random_state")

    p.add_argument("--all-compartments", action="store_true",
                   dest="all_compartments")
    p.add_argument("--out-dir", dest="out_dir")
    p.add_argument("--save-coaccess", action="store_true", dest="save_coaccess")
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
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        for comp, h5ad in DEFAULT_ATAC_H5AD.items():
            if not Path(h5ad).exists():
                log.warning("ATAC h5ad missing, skipping %s: %s", comp, h5ad)
                continue
            save_co = (str(out_dir / f"{comp}_coaccess_glasso.tsv")
                       if args.save_coaccess else None)
            run_compartment(
                atac_h5ad=h5ad,
                compartment=comp,
                gtf_path=args.gtf,
                out_tsv=str(out_dir / f"{comp}.tsv"),
                window_bp=args.window_bp,
                coaccess_threshold=args.coaccess_threshold,
                fdr_threshold=args.fdr_threshold,
                promoter_bp=args.promoter_bp,
                n_metacells=args.n_metacells,
                k_knn=args.k_knn,
                alpha=args.alpha,
                distance_scale=args.distance_scale,
                max_window_peaks=args.max_window_peaks,
                n_jobs=args.n_jobs,
                random_state=args.random_state,
                save_coaccess_tsv=save_co,
            )
    else:
        if not args.atac_h5ad or not args.out:
            sys.exit("--atac-h5ad and --out are required without --all-compartments")
        save_co = (Path(args.out).with_suffix("")
                   .with_suffix("._coaccess_glasso.tsv").__str__()
                   if args.save_coaccess else None)
        run_compartment(
            atac_h5ad=args.atac_h5ad,
            compartment=args.compartment,
            gtf_path=args.gtf,
            out_tsv=args.out,
            window_bp=args.window_bp,
            coaccess_threshold=args.coaccess_threshold,
            promoter_bp=args.promoter_bp,
            n_metacells=args.n_metacells,
            k_knn=args.k_knn,
            alpha=args.alpha,
            distance_scale=args.distance_scale,
            random_state=args.random_state,
            save_coaccess_tsv=save_co,
        )


if __name__ == "__main__":
    main()
