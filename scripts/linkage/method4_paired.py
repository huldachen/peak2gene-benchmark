"""
Method 4A - Paired peak-gene Pearson correlation on Multiome cells.

For every peak x gene pair within +/- window_bp of the gene's TSS, computes
the Pearson correlation between peak accessibility and gene expression
across KNN-aggregated meta-cells drawn from paired Multiome nuclei. Returns
links with ``|r| >= r_threshold`` and per-gene FDR <= ``fdr_threshold``.

Usage
-----
    python scripts/linkage/method4_paired.py \\
        --rna-h5ad data/processed/hickey2023/multiome_rna_qc.h5ad \\
        --atac-h5ad data/processed/hickey2023/atac_immune.h5ad \\
        --compartment immune \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --out results/linkage/method4_paired/immune.tsv

    python scripts/linkage/method4_paired.py --all-compartments \\
        --rna-h5ad data/processed/hickey2023/multiome_rna_qc.h5ad \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --out-dir results/linkage/method4_paired/
"""
from __future__ import annotations

import argparse
import gzip
import logging
import re as _re
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

# RNA obs index uses `_` between sample_id and raw barcode; ATAC uses `#`
RNA_SEP, ATAC_SEP = "_", "#"

_GEM_SUFFIX_RE = _re.compile(r"_\d+$")


def _open_gtf(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def extract_tss_frame(gtf_path: str, *, protein_coding_only: bool = True) -> pd.DataFrame:
    """Parse GENCODE GTF, one TSS per gene (canonical chroms)."""
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


def _split_barcodes(obs_names: pd.Index, sep: str) -> pd.DataFrame:
    """Split ``{sample_id}{sep}{barcode}`` into (sample_id, barcode).

    Strips any trailing ``_<digits>`` gem-group suffix so the RNA and ATAC
    barcodes of the same nucleus match.
    """
    out = []
    for n in obs_names:
        if sep not in n:
            out.append((pd.NA, n))
            continue
        s, b = n.split(sep, 1)
        b = _GEM_SUFFIX_RE.sub("", b)
        out.append((s, b))
    df = pd.DataFrame(out, columns=["sample_id", "barcode"])
    df["full"] = obs_names
    return df


def _pair_rna_atac(rna, atac):
    """Subset both AnnDatas to the intersection of (sample_id, barcode)."""
    rna_idx  = _split_barcodes(rna.obs_names,  RNA_SEP)
    atac_idx = _split_barcodes(atac.obs_names, ATAC_SEP)

    rna_idx["key"]  = rna_idx["sample_id"].astype(str)  + "|" + rna_idx["barcode"].astype(str)
    atac_idx["key"] = atac_idx["sample_id"].astype(str) + "|" + atac_idx["barcode"].astype(str)

    shared = pd.Index(rna_idx["key"]).intersection(pd.Index(atac_idx["key"]))
    log.info("Shared paired cells: %d (of %d RNA, %d ATAC)",
             len(shared), len(rna_idx), len(atac_idx))

    if len(shared) == 0:
        raise RuntimeError(
            "No shared paired cells between RNA and ATAC. Barcode separators "
            f"expected: RNA `{RNA_SEP}`, ATAC `{ATAC_SEP}`. Check obs index format."
        )

    rna_keep  = rna[rna.obs_names[rna_idx["key"].isin(shared).values]].copy()
    atac_keep = atac[atac.obs_names[atac_idx["key"].isin(shared).values]].copy()

    # Re-align on the suffix-stripped key so row i matches row i
    rna_keep.obs["_pair_key"] = (
        rna_idx.loc[rna_idx["full"].isin(rna_keep.obs_names.values), "key"].values
    )
    atac_keep.obs["_pair_key"] = (
        atac_idx.loc[atac_idx["full"].isin(atac_keep.obs_names.values), "key"].values
    )

    rna_keep = rna_keep[rna_keep.obs.sort_values("_pair_key").index].copy()
    atac_keep = atac_keep[atac_keep.obs.sort_values("_pair_key").index].copy()

    assert (rna_keep.obs["_pair_key"].values
            == atac_keep.obs["_pair_key"].values).all(), \
        "Row alignment failed after pairing."

    return rna_keep, atac_keep


def lsi_embed(
    X_bin: sp.csr_matrix,
    n_components: int = 50,
    random_state: int = 42,
) -> np.ndarray:
    """TF-IDF + truncated SVD, drop component 1 (depth-correlated)."""
    log.info("LSI: TF-IDF on binarised matrix %s", X_bin.shape)
    row_sum = np.asarray(X_bin.sum(axis=1)).ravel() + 1e-6
    col_sum = np.asarray(X_bin.sum(axis=0)).ravel() + 1e-6
    tf = sp.diags(1.0 / row_sum) @ X_bin
    idf = np.log(X_bin.shape[0] / col_sum)
    tfidf = tf.multiply(idf).tocsr()
    log.info("LSI: truncated SVD (n_components = %d)", n_components)
    svd = TruncatedSVD(n_components=n_components, random_state=random_state)
    emb = svd.fit_transform(tfidf)
    return emb[:, 1:]


def build_metacell_indicator(
    lsi_emb: np.ndarray,
    *,
    n_metacells: int = 2000,
    k_knn: int = 50,
    random_state: int = 42,
) -> sp.csr_matrix:
    """Return (n_metacells, n_cells) sparse indicator matrix."""
    rng = np.random.default_rng(random_state)
    n_cells = lsi_emb.shape[0]
    n_metacells = min(n_metacells, n_cells)
    log.info("KNN: k=%d in %d-D LSI space", k_knn, lsi_emb.shape[1])
    nn = NearestNeighbors(n_neighbors=k_knn, algorithm="auto", n_jobs=-1)
    nn.fit(lsi_emb)
    seeds = rng.choice(n_cells, size=n_metacells, replace=False)
    _, neighbors = nn.kneighbors(lsi_emb[seeds])
    rows = np.repeat(np.arange(n_metacells), k_knn)
    cols = neighbors.ravel()
    data = np.ones_like(rows, dtype=np.float32)
    return sp.csr_matrix((data, (rows, cols)), shape=(n_metacells, n_cells))


def aggregate_metacells(X, indicator: sp.csr_matrix) -> np.ndarray:
    """indicator (n_meta, n_cells) @ X (n_cells, n_feats) -> dense float32."""
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    meta = indicator @ X.astype(np.float32)
    return np.asarray(meta.todense(), dtype=np.float32)


def standardise_cols(M: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    mu = M.mean(axis=0, keepdims=True)
    sd = M.std(axis=0, keepdims=True) + eps
    out = (M - mu) / sd
    out[:, sd.ravel() < eps * 10] = 0.0
    return out.astype(np.float32)


def peak_gene_correlation(
    atac_meta_std: np.ndarray,
    rna_meta_std:  np.ndarray,
    rna_meta_raw:  np.ndarray,
    peaks_df: pd.DataFrame,
    gene_idx_map: dict[str, int],
    tss_df: pd.DataFrame,
    *,
    compartment: str,
    window_bp: int = 500_000,
    r_threshold: float = 0.45,
    fdr_threshold: float = 0.05,
    min_meta_expr: int = 5,
) -> pd.DataFrame:
    """Compute Pearson r and FDR-adjusted p-values for all (gene, peak) pairs."""
    t0 = time.time()
    n_meta = atac_meta_std.shape[0]
    inv_df = 1.0 / (n_meta - 1)

    peaks_bt = pybedtools.BedTool.from_dataframe(
        peaks_df[["chrom", "start", "end"]]
    ).sort()
    slop_df = tss_df.copy()
    slop_df["start"] = (slop_df["tss"] - window_bp).clip(lower=0)
    slop_df["end"]   = slop_df["tss"] + window_bp
    slop_bt = pybedtools.BedTool.from_dataframe(
        slop_df[["chrom", "start", "end", "gene_id", "gene_name", "tss", "strand"]]
    ).sort()

    log.info("[%s] bedtools intersect: %d peaks x %d TSS+/-%dbp windows",
             compartment, len(peaks_df), len(tss_df), window_bp)
    hits = slop_bt.intersect(peaks_bt, wa=True, wb=True)

    rows = []
    peak_idx_by_coord = {
        (r.chrom, int(r.start), int(r.end)): i
        for i, r in peaks_df.iterrows()
    }
    for h in hits:
        p_coord = (h[7], int(h[8]), int(h[9]))
        p_idx = peak_idx_by_coord.get(p_coord)
        if p_idx is None:
            continue
        rows.append((h[3], h[4], int(h[5]), h[6], p_idx,
                     p_coord[0], p_coord[1], p_coord[2]))

    cand = pd.DataFrame(rows, columns=["gene_id", "gene_name", "tss_pos",
                                       "strand", "peak_idx", "chrom",
                                       "peak_start", "peak_end"])
    log.info("[%s] %d candidate peak-gene pairs (%d unique genes)",
             compartment, len(cand), cand["gene_id"].nunique())

    cand = cand[cand["gene_id"].isin(gene_idx_map)].copy()
    gene_cols = cand["gene_id"].map(gene_idx_map).values
    expressed_counts = (rna_meta_raw > 0).sum(axis=0)
    keep_gene = expressed_counts[gene_cols] >= min_meta_expr
    cand = cand.loc[keep_gene].reset_index(drop=True)
    log.info("[%s] %d pairs after dropping genes expressed in < %d meta-cells",
             compartment, len(cand), min_meta_expr)

    if cand.empty:
        return _empty_frame()

    log.info("[%s] computing correlations per gene (vectorised) ...", compartment)
    r_values  = np.zeros(len(cand), dtype=np.float32)
    p_values  = np.ones(len(cand),  dtype=np.float32)

    gene_col = cand["gene_id"].map(gene_idx_map).values.astype(np.int64)
    peak_col = cand["peak_idx"].values.astype(np.int64)

    order = np.argsort(gene_col, kind="stable")
    gene_col_sorted = gene_col[order]
    peak_col_sorted = peak_col[order]

    unique_genes, starts = np.unique(gene_col_sorted, return_index=True)
    ends = np.r_[starts[1:], len(gene_col_sorted)]

    for gi, st, en in zip(unique_genes, starts, ends):
        pk_idx = peak_col_sorted[st:en]
        gene_vec = rna_meta_std[:, gi]
        atac_sub = atac_meta_std[:, pk_idx]
        r_slice = (gene_vec @ atac_sub) * inv_df
        r_clip = np.clip(r_slice, -0.9999, 0.9999)
        t_stat = r_clip * np.sqrt(max(n_meta - 2, 1)) / np.sqrt(1.0 - r_clip ** 2)
        p_slice = 2.0 * stats.t.sf(np.abs(t_stat), df=max(n_meta - 2, 1))
        out_idx = order[st:en]
        r_values[out_idx] = r_slice
        p_values[out_idx] = p_slice

    cand["r"] = r_values
    cand["pvalue"] = p_values

    def _bh(p):
        p = np.asarray(p)
        n = len(p)
        order = np.argsort(p)
        ranked = p[order]
        q = ranked * n / (np.arange(n) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        out = np.empty(n)
        out[order] = q
        return np.clip(out, 0, 1)

    cand["fdr"] = (
        cand.groupby("gene_id")["pvalue"]
            .transform(lambda g: pd.Series(_bh(g.values), index=g.index))
    )

    kept = cand[(cand["r"].abs() >= r_threshold)
                & (cand["fdr"] <= fdr_threshold)].copy()

    kept["peak_id"] = (kept["chrom"].astype(str) + ":"
                       + kept["peak_start"].astype(str) + "-"
                       + kept["peak_end"].astype(str))
    kept["distance_bp"] = (((kept["peak_start"] + kept["peak_end"]) // 2)
                           - kept["tss_pos"]).abs()
    kept["compartment"] = compartment
    kept["score"] = kept["r"]

    col_order = ["peak_id", "gene_id", "gene_name", "chrom", "peak_start",
                 "peak_end", "tss_pos", "strand", "distance_bp",
                 "compartment", "method", "score", "r", "pvalue", "fdr"]
    kept["method"] = kept.get("method", "method4_paired")
    kept = kept[col_order].sort_values(["chrom", "peak_start", "gene_id"]).reset_index(drop=True)

    log.info(
        "[%s] %d peak-gene links pass (|r| >= %.3f AND FDR <= %.3f)  -- %.1f s",
        compartment, len(kept), r_threshold, fdr_threshold, time.time() - t0,
    )
    return kept


def _empty_frame() -> pd.DataFrame:
    cols = ["peak_id", "gene_id", "gene_name", "chrom", "peak_start",
            "peak_end", "tss_pos", "strand", "distance_bp",
            "compartment", "method", "score", "r", "pvalue", "fdr"]
    return pd.DataFrame(columns=cols)


def _peaks_from_h5ad(atac) -> pd.DataFrame:
    """Parse ``chr:start-end`` peak IDs -> DataFrame (chrom, start, end, peak_id)."""
    peaks = pd.DataFrame(
        [pid.replace(":", "-").split("-") for pid in atac.var_names],
        columns=["chrom", "start", "end"],
    )
    peaks["start"] = peaks["start"].astype(int)
    peaks["end"]   = peaks["end"].astype(int)
    peaks["peak_id"] = (peaks["chrom"].astype(str) + ":"
                        + peaks["start"].astype(str) + "-"
                        + peaks["end"].astype(str))
    return peaks


def run_compartment_paired(
    rna_h5ad: str,
    atac_h5ad: str,
    compartment: str,
    gtf_path: str,
    out_tsv: str,
    *,
    window_bp: int = 500_000,
    r_threshold: float = 0.45,
    fdr_threshold: float = 0.05,
    n_metacells: int = 2_000,
    k_knn: int = 50,
    min_meta_expr: int = 5,
    random_state: int = 42,
) -> pd.DataFrame:
    t0 = time.time()
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    log.info("[%s] Loading RNA   = %s", compartment, rna_h5ad)
    log.info("[%s] Loading ATAC  = %s", compartment, atac_h5ad)
    rna  = sc.read_h5ad(rna_h5ad)
    atac = sc.read_h5ad(atac_h5ad)
    rna, atac = _pair_rna_atac(rna, atac)

    peaks_df = _peaks_from_h5ad(atac)
    keep = peaks_df["chrom"].isin(CANONICAL_CHROMS).values
    atac = atac[:, keep].copy()
    peaks_df = peaks_df.loc[keep].reset_index(drop=True)
    log.info("[%s] Paired cells: %d   peaks (canonical): %d",
             compartment, atac.n_obs, atac.n_vars)

    X_atac = atac.X
    if not sp.issparse(X_atac):
        X_atac = sp.csr_matrix(X_atac)
    X_bin = (X_atac > 0).astype(np.float32)
    lsi = lsi_embed(X_bin, random_state=random_state)

    # Shared indicator — paired data, same cells on both sides
    ind = build_metacell_indicator(
        lsi, n_metacells=n_metacells, k_knn=k_knn, random_state=random_state,
    )

    atac_meta = aggregate_metacells(X_bin, ind)
    log.info("[%s] ATAC meta-cell matrix: %s", compartment, atac_meta.shape)
    atac_meta_std = standardise_cols(atac_meta)
    del atac_meta

    X_rna = rna.X
    if not sp.issparse(X_rna):
        X_rna = sp.csr_matrix(X_rna)
    rna_meta_raw = aggregate_metacells(X_rna, ind)
    log.info("[%s] RNA  meta-cell matrix: %s", compartment, rna_meta_raw.shape)
    rna_meta_log = np.log1p(rna_meta_raw)
    rna_meta_std = standardise_cols(rna_meta_log)

    tss_df = extract_tss_frame(gtf_path)
    # Prefer Ensembl gene_id as the join key when available
    var_names_str = rna.var_names.astype(str)
    looks_like_ensg = var_names_str.str.startswith("ENSG").any()
    if looks_like_ensg:
        log.info("rna.var_names look like Ensembl IDs - using gene_id join")
        gene_idx_map = {g: i for i, g in enumerate(var_names_str)}
    elif "gene_id" in rna.var.columns:
        gene_idx_map = {g: i for i, g in enumerate(rna.var["gene_id"].astype(str))}
    elif "gene_ids" in rna.var.columns:
        gene_idx_map = {g: i for i, g in enumerate(rna.var["gene_ids"].astype(str))}
    else:
        log.warning("rna.var has no gene_id column and var_names are symbols - "
                    "falling back to gene_name keys")
        gene_idx_map = {g: i for i, g in enumerate(var_names_str)}
        tss_df = tss_df.assign(gene_id=tss_df["gene_name"])

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
    links["method"] = "method4_paired"

    links.to_csv(out_tsv, sep="\t", index=False)
    log.info("[%s] %d peak-gene links -> %s  (total %.1f s)",
             compartment, len(links), out_tsv, time.time() - t0)
    _summarise(links, compartment=compartment)
    return links


def _summarise(df: pd.DataFrame, *, compartment: str) -> None:
    if df.empty:
        log.warning("[%s] Empty Method 4 output - nothing to summarise.", compartment)
        return
    pcts_d = np.percentile(df["distance_bp"].values, [0, 25, 50, 75, 90, 100])
    pcts_r = np.percentile(df["score"].values,       [0, 25, 50, 75, 90, 100])
    n_pos = int((df["score"] > 0).sum())
    log.info(
        "[%s] -- Summary --\n"
        "  Total links     : %d  (positive r: %d, negative r: %d)\n"
        "  Unique peaks    : %d   (avg %.1f genes / peak)\n"
        "  Unique genes    : %d   (avg %.1f peaks / gene)\n"
        "  r (score)       : min=%.3f  p25=%.3f  p50=%.3f  p75=%.3f  p90=%.3f  max=%.3f\n"
        "  Distance (bp)   : min=%d  p25=%d  p50=%d  p75=%d  p90=%d  max=%d",
        compartment, len(df), n_pos, len(df) - n_pos,
        df["peak_id"].nunique(), len(df) / max(df["peak_id"].nunique(), 1),
        df["gene_id"].nunique(), len(df) / max(df["gene_id"].nunique(), 1),
        pcts_r[0], pcts_r[1], pcts_r[2], pcts_r[3], pcts_r[4], pcts_r[5],
        int(pcts_d[0]), int(pcts_d[1]), int(pcts_d[2]),
        int(pcts_d[3]), int(pcts_d[4]), int(pcts_d[5]),
    )


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--rna-h5ad", required=True, metavar="H5AD",
                   help="Multiome RNA h5ad (all cells; will be subset to paired).")
    p.add_argument("--atac-h5ad", metavar="H5AD",
                   help="Per-compartment Multiome ATAC h5ad (required unless --all-compartments).")
    p.add_argument("--compartment", default="all",
                   help='Compartment label (default "all").')
    p.add_argument("--out", metavar="TSV",
                   help="Output TSV path (required unless --all-compartments).")

    p.add_argument("--gtf", required=True,
                   help="GENCODE annotation GTF (plain or .gz).")

    p.add_argument("--window-bp", type=int, default=500_000, dest="window_bp",
                   help="Half-window (default 500000 = +/-500 kb).")
    p.add_argument("--r-threshold", type=float, default=0.45, dest="r_threshold",
                   help="Absolute Pearson r threshold (default 0.45, ArchR default).")
    p.add_argument("--fdr-threshold", type=float, default=0.05, dest="fdr_threshold",
                   help="Per-gene FDR threshold (default 0.05).")
    p.add_argument("--n-metacells", type=int, default=2_000, dest="n_metacells",
                   help="Number of KNN meta-cells (default 2000).")
    p.add_argument("--k-knn", type=int, default=50, dest="k_knn",
                   help="KNN neighbours per seed (default 50).")
    p.add_argument("--min-meta-expr", type=int, default=5, dest="min_meta_expr",
                   help="Min meta-cells expressing a gene (default 5).")
    p.add_argument("--seed", type=int, default=42, dest="random_state",
                   help="Random seed (default 42).")

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
            run_compartment_paired(
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
                random_state=args.random_state,
            )
    else:
        if not args.atac_h5ad or not args.out:
            sys.exit("--atac-h5ad and --out are required without --all-compartments")
        run_compartment_paired(
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
            random_state=args.random_state,
        )


if __name__ == "__main__":
    main()
