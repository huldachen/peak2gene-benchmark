"""
Method 5 - ENCODE-rE2G supervised peak-gene linkage.

Applies the pretrained ``atac_megamap`` logistic-regression model from
EngreitzLab/ENCODE_rE2G (Gschwind 2025 Nature) to candidate peak-gene pairs.
The model outputs a calibrated probability per pair; the default 0.179
threshold matches the authors' 70%-recall calibration on K562.

Usage
-----
    python scripts/linkage/method5_re2g.py \\
        --atac-h5ad data/processed/hickey2023/atac_immune.h5ad \\
        --compartment immune \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --out results/linkage/method5/immune.tsv

    python scripts/linkage/method5_re2g.py --all-compartments \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --out-dir results/linkage/method5/
"""
from __future__ import annotations

import argparse
import gzip
import logging
import pickle
import sys
import time
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pybedtools
import scanpy as sc
import scipy.sparse as sp

log = logging.getLogger(__name__)

CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

DEFAULT_ATAC_H5AD = {
    "epithelial": "data/processed/hickey2023/atac_colon_epithelial.h5ad",
    "immune":     "data/processed/hickey2023/atac_immune.h5ad",
    "stromal":    "data/processed/hickey2023/atac_stromal.h5ad",
}

DEFAULT_MODEL_DIR = "external/ENCODE_rE2G/models/atac_megamap"
DEFAULT_UBIQUITOUS_TSV = (
    "external/ENCODE_rE2G/resources/external_features/"
    "gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv"
)


def _open_gtf(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def extract_tss_frame(gtf_path: str, *, protein_coding_only: bool = True) -> pd.DataFrame:
    """Parse GENCODE GTF, one TSS per gene (canonical chroms only)."""
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


def _peaks_from_h5ad(adata) -> pd.DataFrame:
    df = pd.DataFrame(
        [pid.replace(":", "-").split("-") for pid in adata.var_names],
        columns=["chrom", "start", "end"],
    )
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    df["peak_id"] = (df["chrom"].astype(str) + ":"
                     + df["start"].astype(str) + "-" + df["end"].astype(str))
    df["mid"] = (df["start"] + df["end"]) // 2
    return df


def _pseudobulk_activity(adata) -> np.ndarray:
    X = adata.X
    if sp.issparse(X):
        return np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
    return np.asarray(X.sum(axis=0)).astype(np.float64)


def _quantile_normalise(x: np.ndarray) -> np.ndarray:
    """Rank-transform to [0, 1]."""
    ranks = np.argsort(np.argsort(x))
    return ranks.astype(np.float64) / max(len(x) - 1, 1)


def load_rE2G_model(model_dir: str):
    """Load the pickled LR model + feature names + threshold from ``model_dir``."""
    mdir = Path(model_dir)
    model_path = mdir / "model.pkl"
    if not model_path.exists():
        raise FileNotFoundError(f"Missing model.pkl in {mdir}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with open(model_path, "rb") as f:
            model = pickle.load(f)

    feature_table = pd.read_csv(mdir / "feature_table.tsv", sep="\t")
    feature_names = list(feature_table["feature"])

    # Threshold encoded in a filename like "threshold_.179" or "threshold_0.179"
    threshold_files = list(mdir.glob("threshold_*"))
    if len(threshold_files) != 1:
        raise FileNotFoundError(f"Expected exactly one threshold_* file in {mdir}")
    name = threshold_files[0].name
    num = name.split("_", 1)[1]
    if num.startswith("."):
        num = "0" + num
    threshold = float(num)

    log.info("Loaded rE2G model: %s, %d features, threshold = %.3f",
             type(model).__name__, len(feature_names), threshold)
    return model, feature_names, threshold


def _power_law_contact(distance_bp, *, gamma=0.87, offset=5000):
    """Nasser 2021 power-law proxy."""
    d = np.abs(np.asarray(distance_bp)).astype(np.float64)
    return (d + offset) ** (-gamma)


def _counts_between(points_sorted: np.ndarray, lo_vals: np.ndarray,
                    hi_vals: np.ndarray) -> np.ndarray:
    """Count points strictly between lo and hi via binary search."""
    left  = np.minimum(lo_vals, hi_vals)
    right = np.maximum(lo_vals, hi_vals)
    right_idx = np.searchsorted(points_sorted, right, side="left")
    left_idx  = np.searchsorted(points_sorted, left,  side="right")
    return np.clip(right_idx - left_idx, 0, None)


def compute_features(
    atac,
    tss_df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    ubiquitous_map: dict[str, int],
    *,
    compartment: str,
    window_bp: int = 500_000,
    promoter_bp: int = 2_000,
    nearby_bp: int = 5_000,
    power_law_gamma: float = 0.87,
    power_law_offset: int = 5_000,
) -> pd.DataFrame:
    """Compute the 8 rE2G features per candidate (enhancer peak, target gene) pair."""
    t0 = time.time()

    activity_raw = _pseudobulk_activity(atac)
    peaks_df = peaks_df.copy()
    peaks_df["activity"] = activity_raw
    peaks_df["activity_q"] = _quantile_normalise(activity_raw)
    log.info("[%s] Pseudobulk activity: min=%.1f p50=%.1f max=%.1f",
             compartment, activity_raw.min(), float(np.median(activity_raw)),
             activity_raw.max())

    slop_df = tss_df.copy()
    slop_df["start"] = (slop_df["tss"] - window_bp).clip(lower=0)
    slop_df["end"]   = slop_df["tss"] + window_bp
    slop_bt = pybedtools.BedTool.from_dataframe(
        slop_df[["chrom", "start", "end", "gene_id", "gene_name", "tss", "strand"]]
    ).sort()
    peaks_bt = pybedtools.BedTool.from_dataframe(
        peaks_df[["chrom", "start", "end", "peak_id"]]
    ).sort()
    log.info("[%s] Intersecting %d peaks with %d TSS +/- %dbp windows",
             compartment, len(peaks_df), len(tss_df), window_bp)
    hits = slop_bt.intersect(peaks_bt, wa=True, wb=True)

    records = []
    for h in hits:
        records.append((h[3], h[4], int(h[5]), h[6], h[7], int(h[8]), int(h[9]), h[10]))
    cand = pd.DataFrame(records, columns=["gene_id", "gene_name", "tss_pos", "strand",
                                          "chrom", "peak_start", "peak_end", "peak_id"])
    cand["peak_mid"]    = (cand["peak_start"] + cand["peak_end"]) // 2
    cand["distance_bp"] = (cand["peak_mid"] - cand["tss_pos"]).abs().astype(np.int64)
    log.info("[%s] %d candidate pairs (avg %.1f per gene)",
             compartment, len(cand),
             len(cand) / max(cand["gene_id"].nunique(), 1))

    cand = cand.merge(
        peaks_df[["peak_id", "activity", "activity_q"]],
        on="peak_id", how="left",
    )
    cand = cand.rename(columns={"activity": "activity_enh",
                                "activity_q": "activity_enh_q"})

    cand["distanceToTSS"] = cand["distance_bp"].astype(np.float64)

    cand["contactFrequency"] = _power_law_contact(
        cand["distance_bp"].values, gamma=power_law_gamma, offset=power_law_offset,
    )

    cand["abc_num"] = cand["activity_enh"] * cand["contactFrequency"]
    per_gene_sum = cand.groupby("gene_id")["abc_num"].transform("sum")
    cand["ABC.Score"] = np.where(per_gene_sum > 0, cand["abc_num"] / per_gene_sum, 0.0)

    prom_df = tss_df.copy()
    prom_df["start"] = (prom_df["tss"] - promoter_bp).clip(lower=0)
    prom_df["end"]   = prom_df["tss"] + promoter_bp
    prom_bt = pybedtools.BedTool.from_dataframe(
        prom_df[["chrom", "start", "end", "gene_id"]]
    ).sort()
    prom_hits = peaks_bt.intersect(prom_bt, wa=True, wb=True)
    prom_records = []
    for h in prom_hits:
        prom_records.append((h[3], h[7]))
    prom_map = pd.DataFrame(prom_records, columns=["peak_id", "gene_id"])
    prom_merge = prom_map.merge(peaks_df[["peak_id", "activity_q"]],
                                on="peak_id", how="left")
    gene_prom_q = (
        prom_merge.groupby("gene_id")["activity_q"].max().rename("normalizedATAC_prom")
    )
    cand = cand.merge(gene_prom_q, on="gene_id", how="left")
    cand["normalizedATAC_prom"] = cand["normalizedATAC_prom"].fillna(0.0)

    log.info("[%s] Computing sumNearbyEnhancers (+/-%d bp) ...", compartment, nearby_bp)
    sum_near = np.zeros(len(peaks_df), dtype=np.float64)
    id_to_pos = {pid: i for i, pid in enumerate(peaks_df["peak_id"])}
    for chrom, grp in peaks_df.sort_values(["chrom", "mid"]).groupby("chrom"):
        mids = grp["mid"].values
        act  = grp["activity"].values
        cum  = np.concatenate([[0.0], np.cumsum(act)])
        lo_idx = np.searchsorted(mids, mids - nearby_bp, side="left")
        hi_idx = np.searchsorted(mids, mids + nearby_bp, side="right")
        s = cum[hi_idx] - cum[lo_idx]
        for i, pid in enumerate(grp["peak_id"].values):
            sum_near[id_to_pos[pid]] = s[i]
    peaks_df_sumnear = peaks_df[["peak_id"]].copy()
    peaks_df_sumnear["sumNearbyEnhancers"] = sum_near
    cand = cand.merge(peaks_df_sumnear, on="peak_id", how="left")

    log.info("[%s] Computing numTSSEnhGene / numCandidateEnhGene ...", compartment)
    n_tss  = np.zeros(len(cand), dtype=np.int32)
    n_cand = np.zeros(len(cand), dtype=np.int32)
    tss_by_chrom = {c: np.sort(g["tss"].values)
                    for c, g in tss_df.groupby("chrom")}
    peaks_by_chrom_mid = {c: np.sort(g["mid"].values)
                          for c, g in peaks_df.groupby("chrom")}
    for chrom, grp in cand.groupby("chrom", sort=False):
        idx = grp.index.values
        lo  = grp["peak_mid"].values
        hi  = grp["tss_pos"].values
        tss_chr = tss_by_chrom.get(chrom, np.array([], dtype=np.int64))
        peak_mid_chr = peaks_by_chrom_mid.get(chrom, np.array([], dtype=np.int64))
        tss_between = _counts_between(tss_chr, lo, hi)
        n_tss[idx]  = np.clip(tss_between - 1, 0, None)  # exclude target TSS
        n_cand[idx] = _counts_between(peak_mid_chr, lo, hi)
    cand["numTSSEnhGene"]       = n_tss
    cand["numCandidateEnhGene"] = n_cand

    cand["ubiquitousExpressedGene"] = (
        cand["gene_name"].map(ubiquitous_map).fillna(0).astype(np.int8)
    )

    log.info(
        "[%s] Features done in %.1fs  candidate pairs: %d  "
        "ubiquitous gene lookups: %d / %d",
        compartment, time.time() - t0, len(cand),
        int(cand["ubiquitousExpressedGene"].sum()),
        cand["gene_id"].nunique(),
    )
    return cand


def apply_rE2G(cand: pd.DataFrame, model, feature_names: list[str],
               threshold: float, epsilon: float = 0.01) -> pd.DataFrame:
    """Run the pretrained rE2G model on the candidate feature matrix.

    Applies rE2G's trained log transform ``log(|X| + epsilon)`` before
    ``predict_proba``.
    """
    X = np.zeros((len(cand), len(feature_names)), dtype=np.float64)
    for i, name in enumerate(feature_names):
        if name in cand.columns:
            X[:, i] = cand[name].astype(np.float64).fillna(0.0).values
        else:
            log.warning("Feature missing in candidate frame, filling 0.0: %s", name)

    X_logged = np.log(np.abs(X) + epsilon)

    log.info("Applying rE2G LR model to %d candidates x %d features  (log-transformed, eps=%g)",
             *X.shape, epsilon)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        probs = model.predict_proba(X_logged)[:, 1].astype(np.float64)

    out = cand.copy()
    out["rE2G_probability"] = probs
    out["rE2G_binary"]      = (probs >= threshold).astype(np.int8)
    return out


def _load_ubiquitous_map(path: str) -> dict[str, int]:
    df = pd.read_csv(path, sep="\t")
    if not {"TargetGene", "is_ubiquitous_uniform"}.issubset(df.columns):
        raise ValueError(f"Unexpected columns in {path}: {df.columns.tolist()}")
    mapped = df["is_ubiquitous_uniform"].astype(str).str.lower().map(
        {"true": 1, "1": 1, "false": 0, "0": 0}
    ).fillna(0).astype(np.int8)
    out = dict(zip(df["TargetGene"].astype(str), mapped))
    log.info("Ubiquitous-gene lookup: %d genes, %d ubiquitous",
             len(out), int(sum(out.values())))
    return out


def run_compartment_re2g(
    atac_h5ad: str,
    compartment: str,
    gtf_path: str,
    model_dir: str,
    ubiquitous_tsv: str,
    out_tsv: str,
    *,
    window_bp: int = 500_000,
    promoter_bp: int = 2_000,
    nearby_bp: int = 5_000,
    keep_all: bool = False,
) -> pd.DataFrame:
    """Full Method-5 pipeline for one compartment."""
    t_start = time.time()
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    model, feature_names, threshold = load_rE2G_model(model_dir)
    ubiquitous_map = _load_ubiquitous_map(ubiquitous_tsv)

    log.info("[%s] Loading %s ...", compartment, atac_h5ad)
    atac = sc.read_h5ad(atac_h5ad)
    peaks_df = _peaks_from_h5ad(atac)
    keep = peaks_df["chrom"].isin(CANONICAL_CHROMS).values
    atac = atac[:, keep].copy()
    peaks_df = peaks_df.loc[keep].reset_index(drop=True)
    log.info("[%s] %d cells x %d peaks (canonical)",
             compartment, atac.n_obs, atac.n_vars)

    tss_df = extract_tss_frame(gtf_path)

    cand = compute_features(
        atac, tss_df, peaks_df, ubiquitous_map,
        compartment=compartment,
        window_bp=window_bp, promoter_bp=promoter_bp, nearby_bp=nearby_bp,
    )

    scored = apply_rE2G(cand, model, feature_names, threshold)

    if keep_all:
        final = scored
    else:
        final = scored[scored["rE2G_binary"] == 1].copy()

    final["method"]      = "method5_re2g"
    final["compartment"] = compartment
    final["score"]       = final["rE2G_probability"]

    col_order = ["peak_id", "gene_id", "gene_name", "chrom", "peak_start",
                 "peak_end", "tss_pos", "strand", "distance_bp",
                 "compartment", "method", "score",
                 "rE2G_probability", "rE2G_binary"]
    final = final[col_order].sort_values(["chrom", "peak_start", "gene_id"]).reset_index(drop=True)
    final.to_csv(out_tsv, sep="\t", index=False)
    log.info("[%s] %d peak-gene links -> %s  (%.1fs total)",
             compartment, len(final), out_tsv, time.time() - t_start)
    _summarise(final, compartment=compartment, threshold=threshold)
    return final


def _summarise(df: pd.DataFrame, *, compartment: str, threshold: float) -> None:
    if df.empty:
        log.warning("[%s] Empty rE2G output - nothing to summarise.", compartment)
        return
    pcts_d = np.percentile(df["distance_bp"].values, [0, 25, 50, 75, 90, 100])
    pcts_s = np.percentile(df["score"].values, [0, 25, 50, 75, 90, 100])
    log.info(
        "[%s] -- Summary (threshold = %.3f) --\n"
        "  Total links  : %d\n"
        "  Unique peaks : %d   (avg %.1f genes / peak)\n"
        "  Unique genes : %d   (avg %.1f peaks / gene)\n"
        "  rE2G prob    : min=%.3f  p25=%.3f  p50=%.3f  p75=%.3f  p90=%.3f  max=%.3f\n"
        "  Distance (bp): min=%d  p25=%d  p50=%d  p75=%d  p90=%d  max=%d",
        compartment, threshold, len(df),
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
    p.add_argument("--compartment", default="all")
    p.add_argument("--out", metavar="TSV",
                   help="Output TSV path (required unless --all-compartments).")

    p.add_argument("--gtf", required=True, help="GENCODE annotation GTF.")
    p.add_argument("--model-dir", default=DEFAULT_MODEL_DIR, dest="model_dir",
                   help="Pretrained rE2G model directory (default: atac_megamap).")
    p.add_argument("--ubiquitous-tsv", default=DEFAULT_UBIQUITOUS_TSV,
                   dest="ubiquitous_tsv",
                   help="Path to rE2G's gene_promoter_class lookup TSV.")

    p.add_argument("--window-bp", type=int, default=500_000, dest="window_bp",
                   help="Candidate-pair window +/- bp (default 500000).")
    p.add_argument("--promoter-bp", type=int, default=2_000, dest="promoter_bp",
                   help="Promoter window for normalizedATAC_prom (default 2000).")
    p.add_argument("--nearby-bp", type=int, default=5_000, dest="nearby_bp",
                   help="+/-bp for sumNearbyEnhancers (default 5000).")

    p.add_argument("--keep-all", action="store_true", dest="keep_all",
                   help="Keep ALL candidate pairs with rE2G probability (no threshold).")

    p.add_argument("--all-compartments", action="store_true", dest="all_compartments",
                   help="Run epithelial + immune + stromal using defaults.")
    p.add_argument("--out-dir", dest="out_dir",
                   help="Output directory when using --all-compartments.")
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
            run_compartment_re2g(
                atac_h5ad=h5ad,
                compartment=comp,
                gtf_path=args.gtf,
                model_dir=args.model_dir,
                ubiquitous_tsv=args.ubiquitous_tsv,
                out_tsv=str(out_dir / f"{comp}.tsv"),
                window_bp=args.window_bp,
                promoter_bp=args.promoter_bp,
                nearby_bp=args.nearby_bp,
                keep_all=args.keep_all,
            )
    else:
        if not args.atac_h5ad or not args.out:
            sys.exit("--atac-h5ad and --out are required without --all-compartments")
        run_compartment_re2g(
            atac_h5ad=args.atac_h5ad,
            compartment=args.compartment,
            gtf_path=args.gtf,
            model_dir=args.model_dir,
            ubiquitous_tsv=args.ubiquitous_tsv,
            out_tsv=args.out,
            window_bp=args.window_bp,
            promoter_bp=args.promoter_bp,
            nearby_bp=args.nearby_bp,
            keep_all=args.keep_all,
        )


if __name__ == "__main__":
    main()
