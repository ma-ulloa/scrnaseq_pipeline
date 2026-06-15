#!/usr/bin/env python3
"""
04_filtering.py
────────────────────────────────────────────────────────────────
Filter a single sample using per-sample QC thresholds from samples.csv.

Thresholds are read from the samples CSV columns:
    mad_mt, mad_counts_lower, mad_counts_upper,
    mad_genes_lower, mad_genes_upper,
    threshold_mt_upper, threshold_counts_upper, threshold_counts_lower,
    threshold_genes_lower, threshold_genes_upper, min_cells_per_gene

An empty cell in any column means that filter is not applied.
MAD-based thresholds are computed on log1p-transformed values for
genes and counts (back-transformed via expm1 for the final cutoff).
────────────────────────────────────────────────────────────────
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(__file__))
from utils.io       import load_data, attach_metadata
from utils.qc_metrics import (
    get_sample_thresholds,
    compute_outlier_flags,
)

sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser(description="Filter single sample by per-sample QC thresholds")
    p.add_argument("--input",     required=True, help="Count matrix (any supported format)")
    p.add_argument("--metrics",   required=True, help="Per-cell QC metrics CSV from 02_qc_stats_report.py")
    p.add_argument("--samples",   required=True, help="Samples CSV containing per-sample threshold columns")
    p.add_argument("--sample_id",       required=True, help="Sample ID")
    p.add_argument("--remove_doublets", type=lambda x: str(x).lower() == "true",
                   default=False, help="Remove predicted doublets (True/False)")
    p.add_argument("--output",          required=True, help="Output .h5ad path")
    return p.parse_args()


def main():
    args = parse_args()

    remove_doublets = args.remove_doublets

    # ── Load AnnData and metadata ─────────────────────────────────────────────
    adata = load_data(args.input)
    adata = attach_metadata(adata, args.samples, args.sample_id)

    # ── Attach QC metric values from pre-computed CSV ─────────────────────────
    # (avoids re-running compute_qc_metrics; compute_outlier_flags needs these cols)
    metrics_df = pd.read_csv(args.metrics, index_col="cell_id")

    metric_cols = [
        "n_genes_by_counts", "log1p_n_genes_by_counts",
        "total_counts",      "log1p_total_counts",
        "pct_counts_mt",     "pct_counts_ribo",
    ]
    for col in metric_cols:
        if col in metrics_df.columns:
            adata.obs[col] = metrics_df.reindex(adata.obs_names)[col].values

    missing = [c for c in metric_cols if c not in adata.obs.columns]
    if missing:
        sys.exit(
            f"[ERROR] Required QC metric columns missing in {args.metrics}: {missing}\n"
            "Re-run 02_qc_stats_report.py to regenerate the metrics CSV."
        )

    # ── Read per-sample thresholds from samples.csv ───────────────────────────
    thresholds = get_sample_thresholds(args.samples, args.sample_id)
    active = {k: v for k, v in thresholds.items() if v is not None}
    print(f"[INFO] Active thresholds for {args.sample_id}: {active}")

    # ── Compute outlier flags ─────────────────────────────────────────────────
    adata = compute_outlier_flags(adata, thresholds)

    # ── Filter cells ──────────────────────────────────────────────────────────
    n_before = adata.n_obs
    adata    = adata[~adata.obs["outlier"]].copy()
    n_after  = adata.n_obs
    print(
        f"[INFO] Cell filter: {n_before:,} → {n_after:,} "
        f"({n_before - n_after:,} removed, {n_after / n_before * 100:.1f}% retained)"
    )

    # ── Doublet removal ───────────────────────────────────────────────────────
    if remove_doublets:
        if "predicted_doublet" in metrics_df.columns:
            adata.obs["predicted_doublet"] = (
                metrics_df.reindex(adata.obs_names)["predicted_doublet"]
                .fillna(False)
                .astype(bool)
                .values
            )
            n_before = adata.n_obs
            adata    = adata[~adata.obs["predicted_doublet"]].copy()
            print(
                f"[INFO] Doublet removal: {n_before:,} → {adata.n_obs:,} "
                f"({n_before - adata.n_obs:,} doublets removed)"
            )
        else:
            print(
                "[WARN] remove_doublets=true but 'predicted_doublet' column is absent "
                "from the metrics CSV — skipping doublet removal. "
                "Re-run 02_qc_stats_report.py with run_scrublet: true."
            )

    # ── Filter genes ──────────────────────────────────────────────────────────
    min_cells = thresholds.get("min_cells_per_gene")
    if min_cells is not None:
        n_genes_before = adata.n_vars
        sc.pp.filter_genes(adata, min_cells=int(min_cells))
        print(
            f"[INFO] Gene filter (min_cells={int(min_cells)}): "
            f"{n_genes_before:,} → {adata.n_vars:,} genes"
        )

    # ── Write output ──────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"[INFO] Filtered AnnData written → {args.output}")


if __name__ == "__main__":
    main()
