"""
Method 3 - Activity-by-Contact (ABC) peak-gene linkage.

Implementation of the ABC model (Fulco 2019 Nat Genet) with the Nasser 2021
power-law Hi-C proxy. Activity uses the pseudobulk ATAC fragment count;
contact defaults to ``(|d| + offset)^(-gamma)`` with gamma=0.87 and can
optionally come from a real Hi-C .cool file. Per-gene ABC scores are
normalised to sum to 1.0, and links with score >= ``abc_threshold``
(default 0.02) are kept.

Usage
-----
    python scripts/linkage/method3_abc.py \\
        --atac-h5ad data/processed/hickey2023/atac_immune.h5ad \\
        --compartment immune \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --out results/linkage/method3/immune.tsv

    python scripts/linkage/method3_abc.py --all-compartments \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --out-dir results/linkage/method3/
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


def _power_law_contact(distance_bp: np.ndarray, *, gamma: float = 1.0,
                       offset: int = 5_000) -> np.ndarray:
    """Nasser 2021 power-law Hi-C proxy: C(d) = (|d| + offset)^(-gamma)."""
    d = np.abs(distance_bp).astype(np.float64)
    return (d + offset) ** (-gamma)


def _hic_contact_lookup(
    hic_cool_path: str,
    pairs: pd.DataFrame,
) -> np.ndarray:
    """Fetch balanced Hi-C contact frequency for each (chrom, pos_a, pos_b)."""
    try:
        import cooler
    except ImportError as exc:
        raise RuntimeError(
            "--hic-cool requires the `cooler` package (pip install cooler)."
        ) from exc

    clr = cooler.Cooler(hic_cool_path)
    res = clr.binsize
    balanced_mat = clr.matrix(balance=True)

    out = np.zeros(len(pairs), dtype=np.float64)
    for chrom, group in pairs.groupby("chrom"):
        try:
            M = balanced_mat.fetch(chrom)
        except KeyError:
            log.warning("Chromosome %s not in Hi-C cooler; 0 contact", chrom)
            continue
        bin_a = (group["pos_a"].values // res).astype(int)
        bin_b = (group["pos_b"].values // res).astype(int)
        bin_a = np.clip(bin_a, 0, M.shape[0] - 1)
        bin_b = np.clip(bin_b, 0, M.shape[1] - 1)
        vals  = M[bin_a, bin_b]
        vals  = np.nan_to_num(vals, nan=0.0)
        out[group.index.values] = vals
    return out


def _pseudobulk_activity(adata) -> np.ndarray:
    """Sum accessibility across all cells -> 1-D array of length n_peaks."""
    X = adata.X
    if sp.issparse(X):
        return np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
    return X.sum(axis=0).astype(np.float64)


def compute_abc(
    atac_h5ad: str,
    compartment: str,
    gtf_path: str,
    *,
    window_bp: int = 5_000_000,
    abc_threshold: float = 0.02,
    hic_cool: Optional[str] = None,
    power_law_gamma: float = 0.87,
    power_law_offset: int = 5_000,
) -> pd.DataFrame:
    """Full ABC pipeline for one compartment; returns links DataFrame."""
    t0 = time.time()

    log.info("[%s] Loading %s ...", compartment, atac_h5ad)
    adata = sc.read_h5ad(atac_h5ad)

    peaks_df = pd.DataFrame(
        [pid.replace(":", "-").split("-") for pid in adata.var_names],
        columns=["chrom", "start", "end"],
    )
    peaks_df["start"] = peaks_df["start"].astype(int)
    peaks_df["end"]   = peaks_df["end"].astype(int)
    peaks_df["mid"]   = (peaks_df["start"] + peaks_df["end"]) // 2

    keep = peaks_df["chrom"].isin(CANONICAL_CHROMS).values
    adata = adata[:, keep].copy()
    peaks_df = peaks_df.loc[keep].reset_index(drop=True)
    peaks_df["peak_id"] = (peaks_df["chrom"].astype(str) + ":"
                           + peaks_df["start"].astype(str) + "-"
                           + peaks_df["end"].astype(str))
    log.info("[%s] %d cells x %d peaks (canonical chroms)",
             compartment, adata.n_obs, adata.n_vars)

    activity = _pseudobulk_activity(adata)
    peaks_df["activity"] = activity
    peaks_df["peak_idx"] = np.arange(len(peaks_df), dtype=np.int64)
    log.info("[%s] Activity stats: min=%.1f median=%.1f max=%.1f",
             compartment, activity.min(),
             float(np.median(activity)), activity.max())

    tss_df = extract_tss_frame(gtf_path)

    slop_df = tss_df.copy()
    slop_df["start"] = (slop_df["tss"] - window_bp).clip(lower=0)
    slop_df["end"]   = slop_df["tss"] + window_bp
    slop_bt = pybedtools.BedTool.from_dataframe(
        slop_df[["chrom", "start", "end", "gene_id", "gene_name", "tss", "strand"]]
    ).sort()
    peaks_bt = pybedtools.BedTool.from_dataframe(
        peaks_df[["chrom", "start", "end", "peak_idx"]]
    ).sort()

    log.info("[%s] Intersecting %d peaks with %d TSS +/- %d bp windows ...",
             compartment, len(peaks_df), len(tss_df), window_bp)
    hits = slop_bt.intersect(peaks_bt, wa=True, wb=True)

    records = []
    for h in hits:
        records.append((int(h[10]), h[3], h[4], int(h[5]), h[6]))

    if not records:
        log.warning("[%s] No candidate peak-gene pairs found.", compartment)
        return _empty_frame()

    cand = pd.DataFrame(records, columns=["peak_idx", "gene_id", "gene_name",
                                          "tss_pos", "strand"])
    cand = cand.merge(
        peaks_df[["peak_idx", "peak_id", "chrom", "start", "end", "mid", "activity"]],
        on="peak_idx", how="left",
    )
    cand = cand.rename(columns={"start": "peak_start", "end": "peak_end"})
    cand["distance_bp"] = (cand["mid"] - cand["tss_pos"]).abs().astype(np.int64)
    log.info("[%s] %d (peak, gene) candidate pairs (avg %.1f per gene)",
             compartment, len(cand),
             len(cand) / max(cand["gene_id"].nunique(), 1))

    if hic_cool is not None:
        log.info("[%s] Fetching Hi-C contact from %s ...", compartment, hic_cool)
        pairs = pd.DataFrame({
            "chrom": cand["chrom"].values,
            "pos_a": cand["mid"].values,
            "pos_b": cand["tss_pos"].values,
        })
        pairs.index = cand.index
        cand["contact"] = _hic_contact_lookup(hic_cool, pairs)
    else:
        log.info("[%s] Using power-law Hi-C proxy (gamma=%.2f, offset=%d)",
                 compartment, power_law_gamma, power_law_offset)
        cand["contact"] = _power_law_contact(cand["distance_bp"].values,
                                             gamma=power_law_gamma,
                                             offset=power_law_offset)

    cand["abc_numerator"] = cand["activity"] * cand["contact"]
    per_gene_sum = cand.groupby("gene_id")["abc_numerator"].transform("sum")
    cand["score"] = np.where(per_gene_sum > 0,
                             cand["abc_numerator"] / per_gene_sum,
                             0.0)
    log.info("[%s] Pre-threshold ABC: mean=%.4f median=%.4f max=%.4f",
             compartment,
             float(cand["score"].mean()),
             float(cand["score"].median()),
             float(cand["score"].max()))

    kept = cand.loc[cand["score"] >= abc_threshold].copy()
    kept["compartment"] = compartment
    kept["method"]      = "method3_abc"

    col_order = ["peak_id", "gene_id", "gene_name", "chrom", "peak_start",
                 "peak_end", "tss_pos", "strand", "distance_bp",
                 "compartment", "method", "score"]
    kept = kept[col_order].sort_values(["chrom", "peak_start", "gene_id"]).reset_index(drop=True)

    log.info("[%s] %d peak-gene links pass ABC >= %.3f (of %d candidates) in %.1fs",
             compartment, len(kept), abc_threshold, len(cand), time.time() - t0)
    return kept


def _empty_frame() -> pd.DataFrame:
    cols = ["peak_id", "gene_id", "gene_name", "chrom", "peak_start",
            "peak_end", "tss_pos", "strand", "distance_bp",
            "compartment", "method", "score"]
    return pd.DataFrame(columns=cols)


def _summarise(df: pd.DataFrame, *, compartment: str) -> None:
    if df.empty:
        log.warning("[%s] Empty ABC output - nothing to summarise.", compartment)
        return
    pcts_d = np.percentile(df["distance_bp"].values, [0, 25, 50, 75, 90, 100])
    pcts_s = np.percentile(df["score"].values,       [0, 25, 50, 75, 90, 100])
    log.info(
        "[%s] -- Summary --\n"
        "  Total links  : %d\n"
        "  Unique peaks : %d   (avg %.1f genes / peak)\n"
        "  Unique genes : %d   (avg %.1f peaks / gene)\n"
        "  ABC score    : min=%.3f  p25=%.3f  p50=%.3f  p75=%.3f  p90=%.3f  max=%.3f\n"
        "  Distance (bp): min=%d  p25=%d  p50=%d  p75=%d  p90=%d  max=%d",
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
    window_bp: int = 5_000_000,
    abc_threshold: float = 0.02,
    hic_cool: Optional[str] = None,
    power_law_gamma: float = 0.87,
    power_law_offset: int = 5_000,
) -> pd.DataFrame:
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    df = compute_abc(
        atac_h5ad=atac_h5ad,
        compartment=compartment,
        gtf_path=gtf_path,
        window_bp=window_bp,
        abc_threshold=abc_threshold,
        hic_cool=hic_cool,
        power_law_gamma=power_law_gamma,
        power_law_offset=power_law_offset,
    )
    df.to_csv(out_tsv, sep="\t", index=False)
    log.info("[%s] %d peak-gene links written to %s", compartment, len(df), out_tsv)
    _summarise(df, compartment=compartment)
    return df


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

    p.add_argument("--window-bp", type=int, default=5_000_000, dest="window_bp",
                   help="Half-window in bp (default 5000000 = +/-5 Mb).")
    p.add_argument("--abc-threshold", type=float, default=0.02,
                   dest="abc_threshold",
                   help="Minimum ABC score (default 0.02 per Fulco 2019).")

    p.add_argument("--hic-cool", dest="hic_cool", default=None,
                   help="Optional Hi-C cooler file. If omitted, uses power-law proxy.")
    p.add_argument("--power-law-gamma", type=float, default=0.87,
                   dest="power_law_gamma",
                   help="Exponent for power-law contact (default 0.87 per Nasser 2021).")
    p.add_argument("--power-law-offset", type=int, default=5_000,
                   dest="power_law_offset",
                   help="Pseudocount offset for power-law contact (default 5000).")

    p.add_argument("--all-compartments", action="store_true",
                   dest="all_compartments",
                   help="Run epithelial, immune, stromal using defaults.")
    p.add_argument("--out-dir", dest="out_dir",
                   help="Output directory (with --all-compartments).")
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
            run_compartment(
                atac_h5ad=h5ad,
                compartment=comp,
                gtf_path=args.gtf,
                out_tsv=str(out_dir / f"{comp}.tsv"),
                window_bp=args.window_bp,
                abc_threshold=args.abc_threshold,
                hic_cool=args.hic_cool,
                power_law_gamma=args.power_law_gamma,
                power_law_offset=args.power_law_offset,
            )
    else:
        if not args.atac_h5ad or not args.out:
            sys.exit("--atac-h5ad and --out are required without --all-compartments")
        run_compartment(
            atac_h5ad=args.atac_h5ad,
            compartment=args.compartment,
            gtf_path=args.gtf,
            out_tsv=args.out,
            window_bp=args.window_bp,
            abc_threshold=args.abc_threshold,
            hic_cool=args.hic_cool,
            power_law_gamma=args.power_law_gamma,
            power_law_offset=args.power_law_offset,
        )


if __name__ == "__main__":
    main()
