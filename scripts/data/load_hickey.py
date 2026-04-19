"""
Load Hickey 2023 gut data from Dryad-deposited files into AnnData objects.

Supports three modes: ``rna`` merges per-sample MTX trios for the RNA cohort
(B001/B004/B005); ``rna-multiome`` does the same for the Multiome cohort's
RNA half (B006/B008-B012); ``atac`` loads either the pre-computed peak
matrices (default) or builds per-sample ATAC via SnapATAC2 from fragments.

Usage
-----
    python scripts/data/load_hickey.py --mode rna
    python scripts/data/load_hickey.py --mode rna-multiome
    python scripts/data/load_hickey.py --mode atac --source peaks
    python scripts/data/load_hickey.py --mode atac --source fragments --samples B001-A-301
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io

log = logging.getLogger(__name__)

ROOT      = Path(__file__).resolve().parents[2]
RAW       = ROOT / "data/raw/hickey2023"
PROCESSED = ROOT / "data/processed/hickey2023"

# Colon samples from RNA Dryad deposit (Multiome=No, donors B001/B004/B005)
RNA_COLON_SAMPLES = [
    "B001-A-301",   # Sigmoid
    "B001-A-401",   # Transverse
    "B004-A-204",   # Ascending
    "B005-A-201",   # Ascending
]

# Colon samples from ATAC Dryad deposit (Multiome=Yes, donors B006/B008-B012).
# Pre-computed peak matrices aggregate these; prefer --source peaks.
ATAC_MULTIOME_COLON_SAMPLES = [
    "B006-A-001", "B006-A-002", "B006-A-101", "B006-A-201",
    "B008-A-001", "B008-A-002", "B008-A-101", "B008-A-201",
    "B009-A-001", "B009-A-101",
    "B010-A-001", "B010-A-002", "B010-A-101", "B010-A-201",
    "B011-A-001", "B011-A-002", "B011-A-101", "B011-A-201",
    "B012-A-001", "B012-A-002", "B012-A-101", "B012-A-201",
]

# Multiome cohort RNA samples (41 IDs matching atac_merged.h5ad cell-set;
# colon + small intestine). Downloaded from the RNA Dryad deposit.
MULTIOME_RNA_SAMPLES = [
    "B006-A-001", "B006-A-002", "B006-A-101", "B006-A-201",
    "B006-A-201-R2",
    "B006-A-301", "B006-A-402", "B006-A-501",
    "B008-A-001", "B008-A-002", "B008-A-101", "B008-A-201",
    "B008-A-301", "B008-A-401", "B008-A-402", "B008-A-501",
    "B009-A-101", "B009-A-405",
    "B010-A-001", "B010-A-002", "B010-A-101", "B010-A-201",
    "B010-A-301", "B010-A-401", "B010-A-405", "B010-A-501",
    "B011-A-001", "B011-A-002", "B011-A-101", "B011-A-201",
    "B011-A-301", "B011-A-401", "B011-A-405", "B011-A-501",
    "B012-A-001", "B012-A-002", "B012-A-101", "B012-A-201",
    "B012-A-301", "B012-A-401", "B012-A-405", "B012-A-501",
]


def _load_sample_meta() -> pd.DataFrame:
    """Load sample_location_metadata.csv indexed by SampleNameRNA."""
    meta_f = RAW / "metadata" / "sample_location_metadata.csv"
    if not meta_f.exists():
        log.warning(f"sample_location_metadata.csv not found at {meta_f}")
        return pd.DataFrame()
    meta = pd.read_csv(meta_f)
    meta.columns = [c.lower() for c in meta.columns]
    idx_col = next((c for c in meta.columns if "samplenamera" in c or c == "samplenameonly"),
                   meta.columns[0])
    meta = meta.set_index(idx_col)
    return meta


def load_mtx_sample(sample_id: str, rna_dir: Path,
                    valid_barcodes: set[str] | None = None) -> ad.AnnData:
    """Load one 10x MTX trio for a single sample into AnnData.

    Dryad deposits the unfiltered (raw) barcode matrix. Passing
    ``valid_barcodes`` (the annotated cell obs_names from the published TSVs)
    filters to real cells before loading the sparse matrix.
    """
    sample_dir = rna_dir / sample_id

    def _find(stem: str) -> Path:
        prefixed = sample_dir / f"{sample_id}_{stem}"
        plain    = sample_dir / stem
        if prefixed.exists():
            return prefixed
        if plain.exists():
            return plain
        raise FileNotFoundError(
            f"Missing {stem} for {sample_id} "
            f"(tried {prefixed.name} and {plain.name})\n"
            f"Run: python scripts/data/download_hickey_dryad.py rna --samples {sample_id}"
        )

    barcodes_f = _find("barcodes.tsv.gz")
    features_f = _find("features.tsv.gz")
    matrix_f   = _find("matrix.mtx.gz")

    barcodes = pd.read_csv(barcodes_f, header=None, names=["barcode"])

    # features.tsv layouts:
    #   RNA-only Cell Ranger (3 cols): gene_id, gene_name, feature_type
    #   10x Multiome (6 cols): + chrom, start, end, and mixes Gene Expression + Peaks rows
    # We keep only Gene Expression here; ATAC peaks come via the peak_matrix path.
    features = pd.read_csv(features_f, header=None, sep="\t")
    if features.shape[1] >= 3:
        base = ["gene_id", "gene_name", "feature_type"]
        extras = [f"extra{i}" for i in range(features.shape[1] - 3)]
        features.columns = base + extras
    else:
        raise RuntimeError(
            f"Unexpected features.tsv layout for {sample_id}: "
            f"{features.shape[1]} columns (expected >=3)"
        )

    if "feature_type" in features.columns:
        gene_mask = features["feature_type"].astype(str).str.startswith("Gene Expression")
        n_total = len(features)
        n_genes = int(gene_mask.sum())
        if n_genes != n_total:
            log.info(f"  {sample_id}: dropping {n_total - n_genes:,} non-Gene-Expression "
                     f"features; keeping {n_genes:,} RNA genes")
    else:
        gene_mask = pd.Series(True, index=features.index)

    # obs_names match published TSVs: "{sample_id}_{raw_barcode}"
    obs_names = np.array([f"{sample_id}_{bc}" for bc in barcodes["barcode"].values])

    if valid_barcodes is not None:
        keep_mask = np.array([n in valid_barcodes for n in obs_names])
        n_total   = len(obs_names)
        keep_idx  = np.where(keep_mask)[0]
        log.info(f"  {sample_id}: filtering {n_total:,} raw barcodes -> "
                 f"{len(keep_idx):,} annotated cells")
    else:
        keep_idx = None
        log.warning(f"  {sample_id}: no valid_barcodes supplied - loading ALL "
                    f"{len(obs_names):,} raw barcodes")

    mat = scipy.io.mmread(matrix_f).T.tocsr()  # genes x cells -> cells x genes

    if not gene_mask.all():
        mat = mat[:, gene_mask.values]
        features = features.loc[gene_mask].reset_index(drop=True)

    if keep_idx is not None:
        mat       = mat[keep_idx, :]
        obs_names = obs_names[keep_idx]

    # Use Ensembl gene_id as var index (stable across samples; avoids symbol
    # collision explosions on outer concat). gene_name kept as a var column.
    var_df = pd.DataFrame(
        {"gene_id":   features["gene_id"].values,
         "gene_name": features["gene_name"].values},
        index=features["gene_id"].values,
    )

    adata = ad.AnnData(
        X=mat,
        obs=pd.DataFrame(index=obs_names),
        var=var_df,
    )
    adata.var_names_make_unique()
    adata.obs["sample_id"] = sample_id
    adata.obs["donor"]     = sample_id.split("-")[0]
    log.info(f"  {sample_id}: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def transfer_rna_cell_types(adata: ad.AnnData, cell_types_dir: Path) -> ad.AnnData:
    """Transfer published RNA cell-type labels from the Hickey Dryad TSVs."""
    ct_files = {
        "colon_epithelial": "Colon_Epithelial_UMAP_CellType.tsv",
        "immune":           "Immune_UMAP_CellType.tsv",
        "stromal":          "Stromal_UMAP_CellType.tsv",
    }

    frames = []
    for compartment, fname in ct_files.items():
        fpath = cell_types_dir / fname
        if not fpath.exists():
            log.warning(f"Cell-type file not found: {fpath} - skipping {compartment}")
            continue
        df = pd.read_csv(fpath, sep="\t")
        df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

        bc_col = next((c for c in df.columns if c in {"cell_id", "cell"}), None)
        ct_col = next((c for c in df.columns
                       if "cell_type" in c or c == "celltype"), None)
        if bc_col is None or ct_col is None:
            log.warning(f"Cannot identify columns in {fname}: {list(df.columns)}")
            continue

        df = df[[bc_col, ct_col]].rename(columns={bc_col: "obs_name",
                                                    ct_col: "cell_type"})
        df["compartment"] = compartment
        frames.append(df)

    if not frames:
        log.warning("No cell-type TSVs found - cell_type set to 'Unknown'")
        adata.obs["cell_type"]   = "Unknown"
        adata.obs["compartment"] = "Unknown"
        return adata

    ct_df = pd.concat(frames, ignore_index=True).drop_duplicates("obs_name")
    lookup = ct_df.set_index("obs_name")

    adata.obs["cell_type"]   = lookup["cell_type"].reindex(adata.obs_names).fillna("Unknown").values
    adata.obs["compartment"] = lookup["compartment"].reindex(adata.obs_names).fillna("Unknown").values

    n = (adata.obs["cell_type"] != "Unknown").sum()
    log.info(f"  RNA cell-type transfer: {n}/{adata.n_obs} annotated ({100*n/adata.n_obs:.1f}%)")
    return adata


def run_rna(samples: list[str], output_name: str = "rna_merged.h5ad") -> None:
    """Load per-sample MTX files, merge, annotate with cell types, save h5ad."""
    rna_dir = RAW / "rna"
    ct_dir  = RAW / "cell_types"
    out_dir = PROCESSED
    out_dir.mkdir(parents=True, exist_ok=True)

    available = [s for s in samples
                 if (rna_dir / s / "matrix.mtx.gz").exists()
                 or (rna_dir / s / f"{s}_matrix.mtx.gz").exists()]
    missing   = sorted(set(samples) - set(available))
    if missing:
        log.warning(f"Samples not yet downloaded: {missing}")
        log.warning(f"Download with: python scripts/data/download_hickey_dryad.py rna "
                    f"--samples {' '.join(missing)}")
    if not available:
        raise RuntimeError("No RNA sample files found. Download first (see above).")

    log.info("Pre-loading valid barcodes from cell-type TSVs...")
    valid_barcodes: set[str] = set()
    for fname in ["Colon_Epithelial_UMAP_CellType.tsv",
                  "Immune_UMAP_CellType.tsv",
                  "Stromal_UMAP_CellType.tsv"]:
        fpath = ct_dir / fname
        if fpath.exists():
            df = pd.read_csv(fpath, sep="\t", usecols=[0])
            valid_barcodes.update(df.iloc[:, 0].astype(str))
    log.info(f"  {len(valid_barcodes):,} unique annotated barcodes across all compartments")

    log.info(f"Loading {len(available)}/{len(samples)} RNA samples: {available}")
    adatas = [load_mtx_sample(s, rna_dir, valid_barcodes=valid_barcodes)
              for s in available]

    merged = ad.concat(adatas, join="outer", index_unique=None)
    merged.var_names_make_unique()
    log.info(f"Merged RNA: {merged.n_obs} cells x {merged.n_vars} genes")

    merged = transfer_rna_cell_types(merged, ct_dir)

    meta = _load_sample_meta()
    if not meta.empty and "location" in meta.columns:
        loc_map = meta["location"][~meta.index.duplicated(keep="first")]
        merged.obs["region"] = (merged.obs["sample_id"]
                                .map(loc_map)
                                .fillna("unknown"))
        log.info("Region annotation attached from sample_location_metadata.csv")

    out_path = out_dir / output_name
    merged.write_h5ad(out_path)
    log.info(f"Saved: {out_path}")

    print(f"\n{'='*60}")
    print(f"RNA merged: {merged.n_obs} cells x {merged.n_vars} genes")
    print(f"Samples : {merged.obs['sample_id'].value_counts().to_dict()}")
    print(f"\nCell types:")
    print(merged.obs["cell_type"].value_counts().head(15).to_string())
    print(f"\n-> {out_path}")


def load_peak_matrix(compartment: str, atac_dir: Path) -> ad.AnnData:
    """Load a pre-computed cells x peaks matrix from the Dryad scATAC dataset."""
    mtx_f   = atac_dir / f"{compartment}_peak_matrix.mtx"
    cells_f = atac_dir / f"{compartment}_peak_matrix_cells.tsv"
    peaks_f = atac_dir / f"{compartment}_peaks.bed"

    for f in [mtx_f, cells_f, peaks_f]:
        if not f.exists():
            raise FileNotFoundError(
                f"Missing: {f}\n"
                f"Download with: python scripts/data/download_hickey_dryad.py peaks"
            )

    cells = pd.read_csv(cells_f, header=None, names=["barcode"])
    peaks = pd.read_csv(peaks_f, sep=r"\s+", header=None,
                        names=["chrom", "start", "end"])
    peak_ids = (peaks["chrom"] + ":" +
                peaks["start"].astype(str) + "-" +
                peaks["end"].astype(str))

    # MTX is peaks x cells; transpose to cells x peaks
    mat = scipy.io.mmread(mtx_f).T.tocsr()

    adata = ad.AnnData(
        X=mat,
        obs=pd.DataFrame({"barcode": cells["barcode"].values},
                         index=cells["barcode"].values),
        var=pd.DataFrame({"chrom": peaks["chrom"].values,
                          "start": peaks["start"].values,
                          "end":   peaks["end"].values},
                         index=peak_ids),
    )
    adata.obs["sample_id"] = [bc.split("#")[0] if "#" in bc else "unknown"
                              for bc in adata.obs_names]
    adata.obs["donor"]     = [sid.split("-")[0] for sid in adata.obs["sample_id"]]
    adata.obs["compartment"] = compartment
    log.info(f"  {compartment}: {adata.n_obs} cells x {adata.n_vars} peaks")
    return adata


def transfer_atac_cell_types(adata: ad.AnnData, cell_types_dir: Path,
                             compartment: str) -> ad.AnnData:
    """Transfer ATAC cell-type labels from the Dryad scATAC TSVs."""
    fname_map = {
        "colon_epithelial": "scATAC_multiome_cell_types_epithelial_colon.tsv",
        "immune":           "scATAC_multiome_cell_types_immune.tsv",
        "stromal":          "scATAC_multiome_cell_types_stromal.tsv",
    }
    fname = fname_map.get(compartment)
    if fname is None:
        log.warning(f"No cell-type file mapped for compartment: {compartment}")
        adata.obs["cell_type"] = "Unknown"
        return adata

    fpath = cell_types_dir / fname
    if not fpath.exists():
        log.warning(f"Cell-type file not found: {fpath}")
        adata.obs["cell_type"] = "Unknown"
        return adata

    df = pd.read_csv(fpath, sep="\t")
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    bc_col = next((c for c in df.columns if c in {"cell", "cell_id", "barcode"}), None)
    ct_col = next((c for c in df.columns if "cell_type" in c or c == "celltype"), None)

    if bc_col is None or ct_col is None:
        log.warning(f"Unexpected columns in {fname}: {list(df.columns)}")
        adata.obs["cell_type"] = "Unknown"
        return adata

    lookup = df.set_index(bc_col)[ct_col]
    adata.obs["cell_type"] = lookup.reindex(adata.obs_names).fillna("Unknown").values

    n = (adata.obs["cell_type"] != "Unknown").sum()
    log.info(f"  ATAC cell-type transfer ({compartment}): {n}/{adata.n_obs} annotated "
             f"({100*n/adata.n_obs:.1f}%)")
    return adata


def run_atac_peaks(compartments: list[str]) -> None:
    """Load pre-computed peak matrices by compartment, merge, annotate, save."""
    atac_dir = RAW / "atac"
    ct_dir   = RAW / "cell_types"
    out_dir  = PROCESSED
    out_dir.mkdir(parents=True, exist_ok=True)

    adatas = []
    for comp in compartments:
        mtx_f = atac_dir / f"{comp}_peak_matrix.mtx"
        if not mtx_f.exists():
            log.warning(f"Peak matrix not found for '{comp}': {mtx_f} - skipping")
            continue
        a = load_peak_matrix(comp, atac_dir)
        a = transfer_atac_cell_types(a, ct_dir, comp)
        adatas.append(a)

    if not adatas:
        raise RuntimeError(
            "No peak matrix files found.\n"
            "Download: python scripts/data/download_hickey_dryad.py peaks"
        )

    # Compartments have different peak sets; outer join fills 0
    merged = ad.concat(adatas, join="outer", index_unique="_")
    merged.X = merged.X.tocsr()
    merged.var_names_make_unique()
    log.info(f"Merged ATAC: {merged.n_obs} cells x {merged.n_vars} peaks")

    meta = _load_sample_meta()
    if not meta.empty and "location" in meta.columns:
        loc_map = meta["location"][~meta.index.duplicated(keep="first")]
        merged.obs["region"] = (merged.obs["sample_id"]
                                .map(loc_map)
                                .fillna("unknown"))

    out_path = out_dir / "atac_merged.h5ad"
    merged.write_h5ad(out_path)
    log.info(f"Saved: {out_path}")

    print(f"\n{'='*60}")
    print(f"ATAC merged: {merged.n_obs} cells x {merged.n_vars} peaks")
    print(f"Compartments loaded: {[a.obs['compartment'].iloc[0] for a in adatas]}")
    print(f"Donors: {sorted(merged.obs['donor'].unique())}")
    print(f"\nCell types:")
    print(merged.obs["cell_type"].value_counts().head(15).to_string())
    print(f"\n-> {out_path}")


def run_atac_fragments(samples: list[str]) -> None:
    """Build per-sample ATAC AnnData from fragment files using SnapATAC2."""
    try:
        import snapatac2 as snap
    except ImportError:
        raise ImportError("snapatac2 not installed. Run: pip install snapatac2")

    atac_dir = RAW / "atac"
    out_dir  = PROCESSED
    out_dir.mkdir(parents=True, exist_ok=True)

    available = [s for s in samples
                 if (atac_dir / s / "atac_fragments.tsv.gz").exists()]
    missing   = sorted(set(samples) - set(available))
    if missing:
        log.warning(f"Fragment files not found for: {missing}")
    if not available:
        raise RuntimeError(
            "No fragment files found. Download first, or use --source peaks instead."
        )

    log.info(f"Building ATAC from fragments for {len(available)} samples")
    for sample in available:
        frag = str(atac_dir / sample / "atac_fragments.tsv.gz")
        tmp  = str(out_dir / f"{sample}_snap.h5ad")
        data = snap.pp.import_data(
            fragment_file=frag,
            genome=snap.genome.hg38,
            file=tmp,
            min_num_fragments=100,
        )
        data.obs["sample_id"] = sample
        data.obs["donor"]     = sample.split("-")[0]
        log.info(f"  {sample}: {data.n_obs} cells imported -> {tmp}")

    log.info("Fragment import done.")


def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  %(levelname)-8s  %(message)s",
                        datefmt="%H:%M:%S")

    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--mode",
                        choices=["rna", "rna-multiome", "atac"],
                        required=True,
                        help="rna = B001/B004/B005 RNA-only cohort; "
                             "rna-multiome = B006/B008-B012 Multiome cohort's RNA half; "
                             "atac = Multiome ATAC (peaks or fragments)")
    parser.add_argument("--source", choices=["peaks", "fragments"], default="peaks",
                        help="ATAC source: pre-computed peak matrices (default) or raw fragments")
    parser.add_argument("--samples", nargs="+", default=None,
                        help="Sample IDs (RNA/fragments mode).")
    parser.add_argument("--compartments", nargs="+",
                        default=["colon_epithelial", "immune", "stromal"],
                        help="Compartments for peak-matrix ATAC mode")
    args = parser.parse_args()

    if args.mode == "rna":
        samples = args.samples or RNA_COLON_SAMPLES
        run_rna(samples, output_name="rna_merged.h5ad")
    elif args.mode == "rna-multiome":
        samples = args.samples or MULTIOME_RNA_SAMPLES
        run_rna(samples, output_name="multiome_rna_merged.h5ad")
    elif args.mode == "atac" and args.source == "peaks":
        run_atac_peaks(args.compartments)
    elif args.mode == "atac" and args.source == "fragments":
        samples = args.samples or ATAC_MULTIOME_COLON_SAMPLES
        run_atac_fragments(samples)


if __name__ == "__main__":
    main()
