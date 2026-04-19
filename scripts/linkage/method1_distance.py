"""
Method 1 - Distance-window baseline linkage.

For every ATAC peak, links all protein-coding gene TSSs within +/- window_kb.
No expression or regulatory weighting; this is the null comparator for the
other linkage methods.

Usage
-----
    python scripts/linkage/method1_distance.py \\
        --peaks data/raw/hickey2023/atac/immune_peaks.tab.bed \\
        --gtf   data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --chrom-sizes data/raw/reference/hg38.chrom.sizes \\
        --compartment immune \\
        --out   results/linkage/method1/immune.tsv

    python scripts/linkage/method1_distance.py --all-compartments \\
        --gtf data/raw/reference/gencode.v44.annotation.gtf.gz \\
        --chrom-sizes data/raw/reference/hg38.chrom.sizes \\
        --out-dir results/linkage/method1/
"""
from __future__ import annotations

import argparse
import gzip
import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pybedtools

log = logging.getLogger(__name__)

CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

DEFAULT_PEAKS = {
    "epithelial": "data/raw/hickey2023/atac/colon_epithelial_peaks.tab.bed",
    "immune":     "data/raw/hickey2023/atac/immune_peaks.tab.bed",
    "stromal":    "data/raw/hickey2023/atac/stromal_peaks.tab.bed",
}


def _open_gtf(gtf_path: str):
    """Return a text-mode file handle for a GTF file (plain or gzip)."""
    p = str(gtf_path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p)


def extract_tss_bed(
    gtf_path: str,
    *,
    protein_coding_only: bool = True,
) -> pybedtools.BedTool:
    """Parse a GENCODE GTF and return a BedTool of gene-level TSS intervals.

    Each TSS is a 1-bp interval. + strand genes use transcript-start, - strand
    genes use transcript-end.
    """
    records: list[tuple] = []

    log.info("Parsing GTF: %s  (protein_coding_only=%s)", gtf_path, protein_coding_only)

    with _open_gtf(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            chrom = fields[0]
            if chrom not in CANONICAL_CHROMS:
                continue

            # GTF is 1-based inclusive; convert to 0-based half-open
            start  = int(fields[3]) - 1
            end    = int(fields[4])
            strand = fields[6]

            attrs: dict[str, str] = {}
            for token in fields[8].split(";"):
                token = token.strip()
                if not token:
                    continue
                kv = token.split(" ", 1)
                if len(kv) == 2:
                    attrs[kv[0]] = kv[1].strip().strip('"')

            if protein_coding_only and attrs.get("gene_type") != "protein_coding":
                continue

            gene_id   = attrs.get("gene_id", ".")
            gene_name = attrs.get("gene_name", gene_id)

            # Strip Ensembl version suffix
            gene_id = gene_id.split(".")[0]

            tss = start if strand == "+" else end - 1
            records.append((chrom, tss, tss + 1, gene_id, gene_name, strand))

    if not records:
        raise ValueError(
            f"No TSS records extracted from {gtf_path}. "
            "Check --no-protein-coding-filter and file format."
        )

    tss_df = (
        pd.DataFrame(records, columns=["chrom", "start", "end",
                                       "gene_id", "gene_name", "strand"])
        .drop_duplicates(subset=["gene_id"])
        .sort_values(["chrom", "start"])
        .reset_index(drop=True)
    )

    log.info(
        "Extracted %d%s gene TSSs",
        len(tss_df),
        " protein-coding" if protein_coding_only else "",
    )
    return pybedtools.BedTool.from_dataframe(tss_df)


def run_distance_linkage(
    peaks_bed: str,
    gtf_path: str,
    chrom_sizes: str,
    out_tsv: str,
    *,
    window_kb: int = 500,
    compartment: str = "all",
    protein_coding_only: bool = True,
) -> pd.DataFrame:
    """Link each ATAC peak to all gene TSSs within +/- window_kb.

    Steps
    -----
    1. Build 1-bp TSS BedTool from GTF.
    2. Slop each TSS symmetrically by window_bp using chrom_sizes.
    3. Intersect slopped-TSS windows with peak BED (wa=True, wb=True).
    4. Recover TSS from slop_start + window_bp (symmetric slop).
    5. distance_bp = |peak_midpoint - tss_pos|.
    """
    window_bp = window_kb * 1_000
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)

    peaks = pybedtools.BedTool(peaks_bed)
    n_peaks = peaks.count()
    log.info("[%s] %d peaks loaded from %s", compartment, n_peaks, peaks_bed)

    tss_bt = extract_tss_bed(gtf_path, protein_coding_only=protein_coding_only)

    log.info(
        "[%s] Slopping TSSs by +/-%d bp, then intersecting with peaks ...",
        compartment, window_bp,
    )
    tss_slopped = tss_bt.slop(b=window_bp, g=chrom_sizes)
    hits = tss_slopped.intersect(peaks, wa=True, wb=True)

    records: list[dict] = []
    for h in hits:
        slop_start = int(h[1])
        strand     = h[5]
        gene_id    = h[3]
        gene_name  = h[4]

        tss_pos = slop_start + window_bp

        p_chrom = h[6]
        p_start = int(h[7])
        p_end   = int(h[8])
        peak_id = f"{p_chrom}:{p_start}-{p_end}"

        peak_mid    = (p_start + p_end) // 2
        distance_bp = abs(peak_mid - tss_pos)

        records.append({
            "peak_id":     peak_id,
            "gene_id":     gene_id,
            "gene_name":   gene_name,
            "chrom":       p_chrom,
            "peak_start":  p_start,
            "peak_end":    p_end,
            "tss_pos":     tss_pos,
            "strand":      strand,
            "distance_bp": distance_bp,
            "compartment": compartment,
            "method":      "method1_distance",
            "score":       1.0,
        })

    df = (
        pd.DataFrame(records)
        .drop_duplicates(subset=["peak_id", "gene_id"])
        .sort_values(["chrom", "peak_start", "gene_id"])
        .reset_index(drop=True)
    )

    df.to_csv(out_tsv, sep="\t", index=False)
    log.info("[%s] %d peak-gene links written to %s", compartment, len(df), out_tsv)
    _print_summary(df, compartment=compartment)
    return df


def _print_summary(df: pd.DataFrame, *, compartment: str = "all") -> None:
    """Log descriptive statistics about a linkage result."""
    if df.empty:
        log.warning("[%s] Empty linkage result - nothing to summarise.", compartment)
        return

    n_links      = len(df)
    n_peaks      = df["peak_id"].nunique()
    n_genes      = df["gene_id"].nunique()
    avg_per_peak = n_links / n_peaks if n_peaks else float("nan")
    avg_per_gene = n_links / n_genes if n_genes else float("nan")

    dist = df["distance_bp"].values
    pcts = np.percentile(dist, [0, 10, 25, 50, 75, 90, 100])

    log.info(
        "[%s] -- Summary --\n"
        "  Total links   : %d\n"
        "  Unique peaks  : %d  (avg %.1f genes / peak)\n"
        "  Unique genes  : %d  (avg %.1f peaks / gene)\n"
        "  Distance dist : min=%d  p10=%d  p25=%d  p50=%d"
        "  p75=%d  p90=%d  max=%d bp",
        compartment,
        n_links,
        n_peaks, avg_per_peak,
        n_genes, avg_per_gene,
        int(pcts[0]), int(pcts[1]), int(pcts[2]), int(pcts[3]),
        int(pcts[4]), int(pcts[5]), int(pcts[6]),
    )


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])

    p.add_argument("--peaks", metavar="BED",
                   help="Tab-delimited ATAC peak BED (required without --all-compartments).")
    p.add_argument("--compartment", default="all", metavar="LABEL",
                   help='Compartment label for output column (default: "all").')
    p.add_argument("--out", metavar="TSV",
                   help="Output TSV path (required without --all-compartments).")

    p.add_argument("--gtf", required=True, metavar="GTF",
                   help="GENCODE annotation GTF (plain or .gz).")
    p.add_argument("--chrom-sizes", required=True, metavar="SIZES",
                   dest="chrom_sizes",
                   help="UCSC chrom sizes file (hg38.chrom.sizes).")

    p.add_argument("--window-kb", type=int, default=500, metavar="KB",
                   dest="window_kb",
                   help="Half-window in kilobases (default: 500).")
    p.add_argument("--no-protein-coding-filter", action="store_true",
                   dest="all_gene_types",
                   help="Include all gene biotypes (default: protein-coding only).")

    p.add_argument("--all-compartments", action="store_true", dest="all_compartments",
                   help=("Run epithelial, immune, and stromal compartments using "
                         "default BED paths. Requires --out-dir."))
    p.add_argument("--out-dir", metavar="DIR", dest="out_dir",
                   help="Output directory when using --all-compartments.")

    return p


def main(argv: Optional[list[str]] = None) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )

    args = _build_parser().parse_args(argv)
    protein_coding = not args.all_gene_types

    if args.all_compartments:
        if not args.out_dir:
            sys.exit("--all-compartments requires --out-dir")
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        for comp, bed_path in DEFAULT_PEAKS.items():
            if not Path(bed_path).exists():
                log.warning("BED not found, skipping %s: %s", comp, bed_path)
                continue
            run_distance_linkage(
                peaks_bed=bed_path,
                gtf_path=args.gtf,
                chrom_sizes=args.chrom_sizes,
                out_tsv=str(out_dir / f"{comp}.tsv"),
                window_kb=args.window_kb,
                compartment=comp,
                protein_coding_only=protein_coding,
            )

    else:
        if not args.peaks:
            sys.exit("--peaks is required unless --all-compartments is set")
        if not args.out:
            sys.exit("--out is required unless --all-compartments / --out-dir is set")

        run_distance_linkage(
            peaks_bed=args.peaks,
            gtf_path=args.gtf,
            chrom_sizes=args.chrom_sizes,
            out_tsv=args.out,
            window_kb=args.window_kb,
            compartment=args.compartment,
            protein_coding_only=protein_coding,
        )


if __name__ == "__main__":
    main()
