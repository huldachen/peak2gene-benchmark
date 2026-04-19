"""
Method 2 - Cicero-style co-accessibility -> peak-gene linkage.

A Python implementation of Cicero (Pliner et al., Mol. Cell 2018): LSI
embedding of binarised single-cell ATAC, KNN meta-cell aggregation,
per-chromosome sliding-window peak-peak Pearson correlation, and promoter-
anchor projection to derive peak-gene links. Emits raw Pearson correlations
(matching ArchR's ``getCoAccessibility``) with a default threshold of 0.25.

Usage
-----
    python scripts/linkage/method2_cicero.py \\
        --atac-h5ad data/processed/hickey2023/atac_immune.h5ad \\
        --compartment immune \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --out results/linkage/method2/immune.tsv

    python scripts/linkage/method2_cicero.py --all-compartments \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --out-dir results/linkage/method2/
"""
from __future__ import annotations

import argparse
import gzip
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
import scipy.stats as stats
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors

log = logging.getLogger(__name__)

CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

DEFAULT_ATAC_H5AD = {
    "epithelial": "data/processed/hickey2023/atac_colon_epithelial.h5ad",
    "immune":     "data/processed/hickey2023/atac_immune.h5ad",
    "stromal":    "data/processed/hickey2023/atac_stromal.h5ad",
}


def _open_gtf(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def extract_tss_frame(gtf_path: str, *, protein_coding_only: bool = True) -> pd.DataFrame:
    """Parse GENCODE GTF, one TSS row per gene (canonical chroms only)."""
    records = []
    with _open_gtf(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = f[0]
            if chrom not in CANONICAL_CHROMS:
                continue
            start, end, strand = int(f[3]) - 1, int(f[4]), f[6]
            attrs: dict[str, str] = {}
            for tok in f[8].split(";"):
                tok = tok.strip()
                if not tok:
                    continue
                kv = tok.split(" ", 1)
                if len(kv) == 2:
                    attrs[kv[0]] = kv[1].strip().strip('"')
            if protein_coding_only and attrs.get("gene_type") != "protein_coding":
                continue
            gene_id   = attrs.get("gene_id", ".").split(".")[0]
            gene_name = attrs.get("gene_name", gene_id)
            tss = start if strand == "+" else end - 1
            records.append((chrom, tss, gene_id, gene_name, strand))
    df = (
        pd.DataFrame(records, columns=["chrom", "tss", "gene_id", "gene_name", "strand"])
        .drop_duplicates(subset=["gene_id"])
        .sort_values(["chrom", "tss"])
        .reset_index(drop=True)
    )
    log.info("Loaded %d protein-coding TSSs from %s", len(df), Path(gtf_path).name)
    return df


def _lsi_embed(
    X_bin: sp.csr_matrix,
    n_components: int = 50,
    random_state: int = 42,
) -> np.ndarray:
    """TF-IDF + truncated SVD (drop component 1, which tracks depth)."""
    log.info("LSI: computing TF-IDF on binarised matrix %s", X_bin.shape)
    row_sum = np.asarray(X_bin.sum(axis=1)).ravel() + 1e-6
    col_sum = np.asarray(X_bin.sum(axis=0)).ravel() + 1e-6
    tf = sp.diags(1.0 / row_sum) @ X_bin
    idf = np.log(X_bin.shape[0] / col_sum)
    tfidf = tf.multiply(idf).tocsr()

    log.info("LSI: truncated SVD (n_components = %d)", n_components)
    svd = TruncatedSVD(n_components=n_components, random_state=random_state)
    emb = svd.fit_transform(tfidf)
    return emb[:, 1:]


def _build_metacells(
    X_bin: sp.csr_matrix,
    lsi_emb: np.ndarray,
    *,
    n_metacells: int = 2000,
    k_knn: int = 50,
    random_state: int = 42,
) -> np.ndarray:
    """Aggregate KNN-neighbourhoods around random seed cells."""
    rng = np.random.default_rng(random_state)
    n_cells = X_bin.shape[0]
    n_metacells = min(n_metacells, n_cells)

    log.info("KNN: fitting k=%d neighbours in %d-D LSI space", k_knn, lsi_emb.shape[1])
    nn = NearestNeighbors(n_neighbors=k_knn, algorithm="auto", n_jobs=-1)
    nn.fit(lsi_emb)

    seed_idx = rng.choice(n_cells, size=n_metacells, replace=False)
    log.info("KNN: querying %d seed cells for neighbourhood indices ...",
             n_metacells)
    _, neighbors = nn.kneighbors(lsi_emb[seed_idx])

    log.info("Aggregating meta-cells (%d x %d) ...",
             n_metacells, X_bin.shape[1])
    rows = np.repeat(np.arange(n_metacells), k_knn)
    cols = neighbors.ravel()
    data = np.ones_like(rows, dtype=np.float32)
    indicator = sp.csr_matrix(
        (data, (rows, cols)), shape=(n_metacells, n_cells),
    )
    meta_X = indicator @ X_bin.astype(np.float32)
    return np.asarray(meta_X.todense(), dtype=np.float32)


def _standardise_columns(M: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """Centre each column to mean-0 and scale to unit variance."""
    mu = M.mean(axis=0, keepdims=True)
    sd = M.std(axis=0, keepdims=True) + eps
    M_std = (M - mu) / sd
    M_std[:, sd.ravel() < eps * 10] = 0.0
    return M_std.astype(np.float32)


def _coaccess_one_chrom(
    M_std: np.ndarray,
    peak_mids: np.ndarray,
    window_bp: int,
    prefilter_r: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Return (i_idx, j_idx, r, n_candidate_pairs) with i<j upper-triangular."""
    n_meta, n_peaks = M_std.shape
    inv_df = 1.0 / (n_meta - 1)

    i_out, j_out, r_out = [], [], []
    n_candidate_pairs = 0
    hi = 0
    for i in range(n_peaks):
        while hi < n_peaks and peak_mids[hi] - peak_mids[i] <= window_bp:
            hi += 1
        k = hi - (i + 1)
        if k <= 0:
            continue
        n_candidate_pairs += k
        r = M_std[:, i] @ M_std[:, i + 1 : hi] * inv_df
        keep = np.where(r >= prefilter_r)[0]
        if keep.size:
            i_out.append(np.full(keep.size, i, dtype=np.int32))
            j_out.append((i + 1 + keep).astype(np.int32))
            r_out.append(r[keep].astype(np.float32))

    if not i_out:
        return (np.empty(0, np.int32), np.empty(0, np.int32),
                np.empty(0, np.float32), n_candidate_pairs)
    return (np.concatenate(i_out), np.concatenate(j_out),
            np.concatenate(r_out), n_candidate_pairs)


def compute_coaccessibility(
    atac_h5ad: str,
    *,
    compartment: str,
    window_bp: int = 500_000,
    coaccess_threshold: float = 0.25,
    fdr_threshold: float = 0.05,
    n_metacells: int = 2000,
    k_knn: int = 50,
    n_lsi_components: int = 50,
    random_state: int = 42,
    prefilter_r: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (coaccess_df, peaks_df). Filters by r, one-sided t-test, BH FDR."""
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

    i_all, j_all, r_all = [], [], []
    total_candidate_pairs = 0
    for chrom, sub in peaks_sorted.groupby("chrom", sort=False):
        global_idx = sub["idx_global"].values
        mids       = sub["mid"].values.astype(np.int64)
        sub_M      = M_std[:, global_idx]
        log.info("[%s] %s: %d peaks, sliding-window correlation ...",
                 compartment, chrom, len(sub))
        i_local, j_local, r_vals, n_cand = _coaccess_one_chrom(
            sub_M, mids, window_bp, prefilter_r=prefilter_r,
        )
        total_candidate_pairs += n_cand
        if i_local.size:
            i_all.append(global_idx[i_local])
            j_all.append(global_idx[j_local])
            r_all.append(r_vals)

    if i_all:
        i_arr = np.concatenate(i_all)
        j_arr = np.concatenate(j_all)
        r_arr = np.concatenate(r_all)
    else:
        i_arr = np.empty(0, np.int64); j_arr = np.empty(0, np.int64); r_arr = np.empty(0, np.float32)

    # One-sided t-test for r > 0
    n_meta = M_std.shape[0]
    df_meta = max(n_meta - 2, 1)
    r_clip = np.clip(r_arr, -0.9999, 0.9999)
    t_stat = r_clip * np.sqrt(df_meta) / np.sqrt(1.0 - r_clip ** 2)
    p_values = stats.t.sf(t_stat, df=df_meta).astype(np.float32)

    # BH FDR uses the full candidate pair count as denominator so prefilter_r
    # doesn't bias q-values of retained pairs.
    fdr_values = _bh_fdr(p_values, n_tests_total=total_candidate_pairs).astype(np.float32)

    log.info(
        "[%s] Correlated %d candidate pairs; %d pass r >= %.3f (prefilter); "
        "%d also pass FDR <= %.3f",
        compartment, total_candidate_pairs, len(r_arr), prefilter_r,
        int((fdr_values <= fdr_threshold).sum()), fdr_threshold,
    )

    keep = (r_arr >= coaccess_threshold) & (fdr_values <= fdr_threshold)
    i_arr  = i_arr[keep];  j_arr  = j_arr[keep]
    r_arr  = r_arr[keep];  p_values = p_values[keep];  fdr_values = fdr_values[keep]

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

    log.info("[%s] Final: %d pairs pass r >= %.3f AND FDR <= %.3f  (total %.1fs)",
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


def coaccess_to_peak_gene_links(
    coaccess_df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    tss_df: pd.DataFrame,
    *,
    compartment: str,
    promoter_bp: int = 2000,
) -> pd.DataFrame:
    """Derive peak-gene links from peak-peak co-accessibility + promoter anchors."""

    peaks_bt = pybedtools.BedTool.from_dataframe(
        peaks_df[["chrom", "start", "end"]]
    ).sort()
    prom_df = tss_df.copy()
    prom_df["start"] = (prom_df["tss"] - promoter_bp).clip(lower=0)
    prom_df["end"]   = prom_df["tss"] + promoter_bp
    prom_bt = pybedtools.BedTool.from_dataframe(
        prom_df[["chrom", "start", "end", "gene_id", "gene_name", "tss", "strand"]]
    ).sort()

    log.info("[%s] Intersecting %d peaks with %d promoter windows (+/-%d bp) ...",
             compartment, len(peaks_df), len(prom_df), promoter_bp)
    hits = peaks_bt.intersect(prom_bt, wa=True, wb=True)

    prom_records = []
    for h in hits:
        p_chrom, p_start, p_end = h[0], int(h[1]), int(h[2])
        g_id, g_name = h[6], h[7]
        tss_pos = int(h[8]); strand = h[9]
        prom_records.append({
            "peak_id":   f"{p_chrom}:{p_start}-{p_end}",
            "chrom":     p_chrom,
            "peak_start": p_start,
            "peak_end":   p_end,
            "gene_id":   g_id,
            "gene_name": g_name,
            "tss_pos":   tss_pos,
            "strand":    strand,
        })
    prom_map = pd.DataFrame(prom_records)
    if prom_map.empty:
        log.warning("[%s] No promoter peaks - output will be empty.", compartment)
        return pd.DataFrame()

    log.info(
        "[%s] %d (peak, gene) promoter assignments "
        "(%d unique peaks, %d unique genes)",
        compartment, len(prom_map),
        prom_map["peak_id"].nunique(), prom_map["gene_id"].nunique(),
    )

    # Peak1 is the promoter -> link Peak2 to Peak1's gene(s); and vice versa
    a = coaccess_df.merge(prom_map, left_on="Peak1", right_on="peak_id", how="inner")
    a = a.rename(columns={"Peak2": "other_peak"})
    b = coaccess_df.merge(prom_map, left_on="Peak2", right_on="peak_id", how="inner")
    b = b.rename(columns={"Peak1": "other_peak"})

    trans = pd.concat([a, b], ignore_index=True)
    log.info("[%s] %d directional peak->gene edges from co-accessibility pairs",
             compartment, len(trans))

    other = trans["other_peak"].str.extract(
        r"^(?P<dist_chrom>[^:]+):(?P<dist_start>\d+)-(?P<dist_end>\d+)$"
    )
    trans["dist_chrom"] = other["dist_chrom"]
    trans["dist_start"] = other["dist_start"].astype(int)
    trans["dist_end"]   = other["dist_end"].astype(int)

    has_stats = {"pvalue", "fdr"}.issubset(coaccess_df.columns)
    if not has_stats:
        trans["pvalue"] = np.nan
        trans["fdr"]    = np.nan

    distal = pd.DataFrame({
        "peak_id":     trans["other_peak"],
        "gene_id":     trans["gene_id"],
        "gene_name":   trans["gene_name"],
        "chrom":       trans["dist_chrom"],
        "peak_start":  trans["dist_start"],
        "peak_end":    trans["dist_end"],
        "tss_pos":     trans["tss_pos"].astype(int),
        "strand":      trans["strand"],
        "distance_bp": (((trans["dist_start"] + trans["dist_end"]) // 2)
                        - trans["tss_pos"]).abs(),
        "score":       trans["coaccess"].astype(float),
        "pvalue":      trans["pvalue"],
        "fdr":         trans["fdr"],
    })

    # Promoter peaks' self-links
    self_prom = pd.DataFrame({
        "peak_id":     prom_map["peak_id"],
        "gene_id":     prom_map["gene_id"],
        "gene_name":   prom_map["gene_name"],
        "chrom":       prom_map["chrom"],
        "peak_start":  prom_map["peak_start"],
        "peak_end":    prom_map["peak_end"],
        "tss_pos":     prom_map["tss_pos"],
        "strand":      prom_map["strand"],
        "distance_bp": (((prom_map["peak_start"] + prom_map["peak_end"]) // 2)
                        - prom_map["tss_pos"]).abs(),
        "score":       1.0,
        "pvalue":      0.0,
        "fdr":         0.0,
    })

    links = pd.concat([self_prom, distal], ignore_index=True)
    links["compartment"] = compartment
    links["method"]      = "method2_cicero"

    links = (
        links.sort_values("score", ascending=False)
             .drop_duplicates(subset=["peak_id", "gene_id"], keep="first")
             .sort_values(["chrom", "peak_start", "gene_id"])
             .reset_index(drop=True)
    )

    col_order = ["peak_id", "gene_id", "gene_name", "chrom", "peak_start",
                 "peak_end", "tss_pos", "strand", "distance_bp",
                 "compartment", "method", "score", "pvalue", "fdr"]
    return links[col_order]


def run_compartment(
    atac_h5ad: str,
    compartment: str,
    gtf_path: str,
    out_tsv: str,
    *,
    window_bp: int = 500_000,
    coaccess_threshold: float = 0.25,
    fdr_threshold: float = 0.05,
    promoter_bp: int = 2000,
    n_metacells: int = 2000,
    k_knn: int = 50,
    random_state: int = 42,
    save_coaccess_tsv: Optional[str] = None,
) -> pd.DataFrame:
    """Full Method-2 pipeline for one compartment."""
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    coaccess_df, peaks_df = compute_coaccessibility(
        atac_h5ad,
        compartment=compartment,
        window_bp=window_bp,
        coaccess_threshold=coaccess_threshold,
        fdr_threshold=fdr_threshold,
        n_metacells=n_metacells,
        k_knn=k_knn,
        random_state=random_state,
    )

    if save_coaccess_tsv:
        Path(save_coaccess_tsv).parent.mkdir(parents=True, exist_ok=True)
        coaccess_df.to_csv(save_coaccess_tsv, sep="\t", index=False)
        log.info("[%s] Saved raw peak-peak co-access to %s",
                 compartment, save_coaccess_tsv)

    tss_df = extract_tss_frame(gtf_path)
    links = coaccess_to_peak_gene_links(
        coaccess_df, peaks_df, tss_df,
        compartment=compartment, promoter_bp=promoter_bp,
    )
    links.to_csv(out_tsv, sep="\t", index=False)
    log.info("[%s] %d peak-gene links -> %s", compartment, len(links), out_tsv)
    _summarise(links, compartment=compartment)
    return links


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
        "  Score  (raw r): min=%.3f  p25=%.3f  p50=%.3f  p75=%.3f  p90=%.3f  max=%.3f\n"
        "  Distance (bp) : min=%d  p25=%d  p50=%d  p75=%d  p90=%d  max=%d",
        compartment, len(df),
        df["peak_id"].nunique(), len(df) / max(df["peak_id"].nunique(), 1),
        df["gene_id"].nunique(), len(df) / max(df["gene_id"].nunique(), 1),
        pcts_s[0], pcts_s[1], pcts_s[2], pcts_s[3], pcts_s[4], pcts_s[5],
        int(pcts_d[0]), int(pcts_d[1]), int(pcts_d[2]),
        int(pcts_d[3]), int(pcts_d[4]), int(pcts_d[5]),
    )


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--atac-h5ad", metavar="H5AD",
                   help="Per-compartment ATAC h5ad (required unless --all-compartments).")
    p.add_argument("--compartment", default="all",
                   help='Compartment label (default "all").')
    p.add_argument("--out", metavar="TSV",
                   help="Output TSV path (required unless --all-compartments).")

    p.add_argument("--gtf", required=True,
                   help="GENCODE annotation GTF (plain or .gz).")

    p.add_argument("--window-bp", type=int, default=500_000, dest="window_bp",
                   help="Co-accessibility window (default 500000 = +/-500 kb).")
    p.add_argument("--coaccess-threshold", type=float, default=0.25,
                   dest="coaccess_threshold",
                   help="Minimum raw-correlation score (default 0.25).")
    p.add_argument("--fdr-threshold", type=float, default=0.05,
                   dest="fdr_threshold",
                   help="BH FDR threshold (default 0.05).")
    p.add_argument("--promoter-bp", type=int, default=2000, dest="promoter_bp",
                   help="Promoter window for TSS-peak anchoring (default 2000 = +/-2 kb).")
    p.add_argument("--n-metacells", type=int, default=2000, dest="n_metacells",
                   help="Number of KNN meta-cells to aggregate (default 2000).")
    p.add_argument("--k-knn", type=int, default=50, dest="k_knn",
                   help="KNN neighbours per seed cell (default 50).")
    p.add_argument("--seed", type=int, default=42, dest="random_state",
                   help="Random seed (default 42).")

    p.add_argument("--all-compartments", action="store_true",
                   dest="all_compartments",
                   help="Run epithelial, immune, stromal using defaults.")
    p.add_argument("--out-dir", dest="out_dir",
                   help="Output directory (with --all-compartments).")
    p.add_argument("--save-coaccess", action="store_true", dest="save_coaccess",
                   help="Also save raw peak-peak co-access TSV next to output.")
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
            save_co = str(out_dir / f"{comp}_coaccess.tsv") if args.save_coaccess else None
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
                random_state=args.random_state,
                save_coaccess_tsv=save_co,
            )
    else:
        if not args.atac_h5ad or not args.out:
            sys.exit("--atac-h5ad and --out are required without --all-compartments")
        save_co = (Path(args.out).with_suffix("")
                   .with_suffix("._coaccess.tsv").__str__()
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
            random_state=args.random_state,
            save_coaccess_tsv=save_co,
        )


if __name__ == "__main__":
    main()
