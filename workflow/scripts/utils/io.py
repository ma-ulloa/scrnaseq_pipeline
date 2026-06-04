"""
utils/io.py
────────────────────────────────────────────────────────────────
Functions for the scRNA-seq pipeline.

Functions:
    load_data(input_path)        — auto-detects format and loads AnnData
    attach_metadata(adata, ...)  — attaches all metadata columns to adata.obs
    get_sample_ids(metadata)     — returns list of sample IDs (for Snakemake)
────────────────────────────────────────────────────────────────
"""

import os
import sys
import pandas as pd
import scanpy as sc


# Load data

def load_data(input_path: str) -> sc.AnnData:
    """
    Auto-detect file format and load into AnnData.
    # Note Ale: I think it makes the pipeline more flexible.

    Supported formats:
        .h5ad        -> AnnData (already processed)
        .h5          ->  tries 10x Genomics first, falls back to generic HDF5 (e.g. PARSE)
        .mtx / dir   ->  10x MEX format (matrix.mtx + barcodes + features)
        .csv         ->  dense count matrix
        .txt / .tsv  ->  tab-separated count matrix
    """
    if not os.path.exists(input_path):
        sys.exit(f"[ERROR] Input not found: {input_path}")

    # Directory <> assume 10x MEX format
    if os.path.isdir(input_path):
        print(f"[INFO] Detected directory — loading as 10x MEX: {input_path}")
        adata = sc.read_10x_mtx(input_path, var_names="gene_symbols")
        adata.var_names_make_unique()
        return adata

    ext = os.path.splitext(input_path)[-1].lower()

    if ext == ".h5ad":
        print(f"[INFO] Detected .h5ad — loading AnnData: {input_path}")
        adata = sc.read_h5ad(input_path)

    elif ext == ".h5":
        adata = _load_h5(input_path)

    elif ext == ".mtx":
        parent_dir = os.path.dirname(input_path)
        print(f"[INFO] Detected .mtx — loading as 10x MEX from: {parent_dir}")
        adata = sc.read_10x_mtx(parent_dir, var_names="gene_symbols")

    elif ext == ".csv":
        print(f"[INFO] Detected .csv — loading dense matrix: {input_path}")
        adata = sc.read_csv(input_path).T   # transpose to cells x genes

    elif ext in (".txt", ".tsv"):
        print(f"[INFO] Detected {ext} — loading tab-separated matrix: {input_path}")
        adata = sc.read_text(input_path).T

    else:
        sys.exit(
            f"[ERROR] Unrecognised file format: '{ext}'.\n"
            f"        Supported: .h5ad, .h5, .mtx, .csv, .txt, .tsv, or a 10x MEX directory."
        )

    adata.var_names_make_unique()
    return adata


def _load_h5(path: str) -> sc.AnnData:
    """
    Try loading a .h5 as 10x Genomics; fall back to generic HDF5 (PARSE/other).
    """
    try:
        print(f"[INFO] Detected .h5 — trying 10x Genomics format: {path}")
        adata = sc.read_10x_h5(path)
        print("[INFO] Loaded successfully as 10x .h5")
        return adata
    except Exception as e_10x:
        print(f"[WARN] 10x loading failed ({e_10x}) — trying generic HDF5 (PARSE/other)")
        try:
            adata = sc.read_h5ad(path)
            print("[INFO] Loaded successfully as generic HDF5")
            return adata
        except Exception as e_generic:
            sys.exit(
                f"[ERROR] Could not load .h5 file.\n"
                f"        10x error:     {e_10x}\n"
                f"        Generic error: {e_generic}"
            )


# Add Metadata 

def attach_metadata(adata: sc.AnnData, metadata_path: str, sample_id: str) -> sc.AnnData:
    """
    Attach all metadata columns to adata.obs dynamically.

    - First column is the sample ID index (Note Ale: idk if this is the best way tho)
    - All remaining columns are attached automatically
    - A consistent 'sample_id' key is always written to adata.obs
    """
    if not os.path.exists(metadata_path):
        sys.exit(f"[ERROR] Metadata file not found: {metadata_path}")

    meta   = pd.read_csv(metadata_path)
    id_col = meta.columns[0]
    meta[id_col] = meta[id_col].astype(str)

    row = meta[meta[id_col] == str(sample_id)]
    if row.empty:
        sys.exit(
            f"[ERROR] sample_id '{sample_id}' not found in metadata column '{id_col}'.\n"
            f"        Available IDs: {meta[id_col].tolist()}"
        )

    extra_cols = [c for c in meta.columns if c != id_col]
    for col in extra_cols:
        adata.obs[col] = str(row[col].values[0])

    adata.obs["sample_id"] = sample_id
    print(f"[INFO] Metadata attached for '{sample_id}': {extra_cols}")
    return adata


