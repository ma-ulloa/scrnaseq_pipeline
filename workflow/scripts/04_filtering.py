#!/usr/bin/env python3
"""
03_filtering.py
────────────────────────────────────────────────────────────────
Filter a single sample using the outlier flags written by
01_qc_report.py.

The outlier columns (outlier_1mad … outlier_5mad) are already in
the per-cell CSV produced by 01_qc_report.py.  This script reads
the chosen n_mads from config.yaml, selects the matching column,
and filters the AnnData — no MAD maths here.

Usage:
    python 03_filtering.py \
        --input     data/raw/sample.h5 \
        --metrics   results/qc/pre_KS2204T1_2_cell_qc_metrics.csv \
        --metadata  config/metadata.csv \
        --cfg       config/config.yaml \
        --sample_id KS2204T1_2 \
        --output    results/filtered/KS2204T1_2_filtered.h5ad

Filtering logic:
    Keeps cells where outlier_{n_mads}mad == False.
    If that column is absent (e.g. n_mads not in 1-5) the script
    falls back to recomputing thresholds and exits with a warning.
────────────────────────────────────────────────────────────────
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

sys.path.insert(0, os.path.dirname(__file__))
from utils.io       import load_data, attach_metadata
from utils.qc_metrics import compute_qc_metrics, thresholds_from_mad, add_outlier_flags

sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser(description="Filter single sample by MAD-derived thresholds")
    p.add_argument("--input",     required=True, help="Count matrix (any supported format)")
    p.add_argument("--metrics",   required=True, help="Per-cell QC metrics CSV from 01_qc_report.py")
    p.add_argument("--samples",  required=True, help="Metadata CSV")
    p.add_argument("--cfg",       required=True, help="config.yaml")
    p.add_argument("--sample_id", required=True, help="Sample ID")
    p.add_argument("--output",    required=True, help="Output .h5ad path")
    return p.parse_args()


def main():
    args = parse_args()

    with open(args.cfg) as f:
        config = yaml.safe_load(f)

    n_mads = config.get("qc_thresholds", {}).get("n_mads", 3)
    flag_col = f"outlier_{int(n_mads)}mad"

    # ── Load AnnData ──────────────────────────────────────────────────────────
    adata = load_data(args.input)
    adata = attach_metadata(adata, args.samples, args.sample_id)

    # ── Attach pre-computed outlier flags from the QC metrics CSV ─────────────
    metrics_df = pd.read_csv(args.metrics, index_col="cell_id")

    outlier_cols = [c for c in metrics_df.columns if c.startswith("outlier_")]
    if not outlier_cols:
        raise ValueError(
            f"No outlier_* columns found in {args.metrics}. "
            "Re-run 01_qc_report.py to regenerate the metrics CSV."
        )

    if flag_col not in metrics_df.columns:
        print(
            f"[WARN] Column '{flag_col}' not found in metrics CSV "
            f"(available: {outlier_cols}). "
            f"Recomputing outlier flags for n_mads={n_mads}."
        )
        adata = compute_qc_metrics(adata)
        adata = add_outlier_flags(adata, mad_values=[n_mads])
    else:
        # Align metrics to adata cell order and attach
        shared_cells = adata.obs_names.intersection(metrics_df.index)
        if len(shared_cells) < adata.n_obs:
            print(
                f"[WARN] {adata.n_obs - len(shared_cells)} cells in AnnData "
                "not found in metrics CSV — they will be treated as outliers."
            )
        for col in outlier_cols:
            adata.obs[col] = metrics_df.reindex(adata.obs_names)[col].fillna(True).values

    # ── Filter ────────────────────────────────────────────────────────────────
    n_before = adata.n_obs
    keep     = ~adata.obs[flag_col].values
    adata    = adata[keep].copy()
    n_after  = adata.n_obs
    n_removed = n_before - n_after

    print(
        f"[INFO] Filtering with {flag_col}: "
        f"{n_before:,} → {n_after:,} cells "
        f"({n_removed:,} removed, {n_after/n_before*100:.1f}% retained)"
    )

    # ── Write output ──────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"[INFO] Filtered AnnData written → {args.output}")


if __name__ == "__main__":
    main()