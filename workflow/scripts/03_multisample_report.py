#!/usr/bin/env python3
"""
02_multisample_report.py
────────────────────────────────────────────────────────────────
Joins per-sample QC metrics CSVs into one cohort-level report.

Loads the per-cell CSV files produced by 01_qc_stats_report.py,
optionally joins metadata columns (patient, timepoint, cells,
response), builds cohort-level figures, and writes a combined
HTML + PDF report plus a master CSV.

Usage (called by Snakemake rule multiqc_summary):
    python 02_multisample_report.py \
        --input_files results/qc/pre_KS2103T1_1_cell_qc_metrics.csv \
                      results/qc/pre_KS2103T1_2_cell_qc_metrics.csv ... \
        --metadata    config/metadata.csv \
        --output_csv  results/qc/pre_cohort_qc_metrics.csv \
        --plot_dir    results/qc/cohort_plots_pre \
        --output_html results/qc/pre_cohort_qc_report.html \
        --stage       pre
────────────────────────────────────────────────────────────────
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import pandas as pd

# ── make utils importable regardless of working directory ──────────────────
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from utils.report import build_html_report, export_pdf
from utils.qc_plots import (
    fig_cohort_bar_plots,
    fig_cohort_scatters,
    fig_cohort_violins,
    fig_metadata_breakdown,
)


# ── CLI ────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Merge per-cell QC metrics and generate cohort-level report"
    )
    p.add_argument(
        "--input_files", nargs="+", required=True,
        help="Space-separated list of per-sample cell_qc_metrics CSV files"
    )
    p.add_argument(
        "--samples", required=False, default=None,
        help="Optional: metadata CSV (first column = sample_id) to join onto merged table"
    )
    p.add_argument(
        "--output_csv", required=True,
        help="Path for the master merged CSV (one row per cell)"
    )
    p.add_argument(
        "--plot_dir", required=True,
        help="Directory where individual plot PNGs are saved"
    )
    p.add_argument(
        "--output_html", required=True,
        help="Path for the combined HTML report"
    )
    p.add_argument(
        "--stage", required=True, choices=["pre", "post"],
        help="Pipeline stage — controls report header colour"
    )
    return p.parse_args()


# ── data loading ────────────────────────────────────────────────────────────

def merge_qc_metrics(file_list: list, stage: str) -> pd.DataFrame:
    """
    Concatenate per-sample cell QC CSVs and tag with the current stage.

    Parameters
    ----------
    file_list : list of paths to per-sample *_cell_qc_metrics.csv files
    stage     : 'pre' or 'post'

    Returns
    -------
    pd.DataFrame with one row per cell, columns:
        sample_id, n_genes_by_counts, log1p_n_genes_by_counts,
        total_counts, log1p_total_counts, pct_counts_mt,
        pct_counts_ribo, pass_mad_filters, stage
    """
    dfs = []
    for path in file_list:
        if not os.path.exists(path):
            print(f"[WARN] File not found, skipping: {path}")
            continue
        df = pd.read_csv(path)
        df["stage"] = stage
        dfs.append(df)

    if not dfs:
        sys.exit("[ERROR] No valid input CSV files found.")

    merged = pd.concat(dfs, ignore_index=True)
    print(f"[INFO] Merged {len(dfs)} files — {len(merged):,} total cells.")
    return merged


def attach_metadata_to_merged(merged_df: pd.DataFrame,
                               metadata_path: str) -> pd.DataFrame:
    """
    Left-join metadata columns (patient, timepoint, cells, response, …)
    onto the merged per-cell DataFrame using sample_id as the key.

    Columns already present in merged_df are not overwritten.
    """
    meta   = pd.read_csv(metadata_path)
    id_col = meta.columns[0]                 # first column is always sample_id
    meta   = meta.rename(columns={id_col: "sample_id"})
    meta["sample_id"] = meta["sample_id"].astype(str)

    extra_cols = [c for c in meta.columns if c != "sample_id"
                  and c not in merged_df.columns]
    if not extra_cols:
        print("[INFO] No new metadata columns to join.")
        return merged_df

    merged_df = merged_df.merge(
        meta[["sample_id"] + extra_cols],
        on="sample_id", how="left"
    )
    print(f"[INFO] Joined metadata columns: {extra_cols}")
    return merged_df


# ── figure naming ───────────────────────────────────────────────────────────

def _save_figures(figures: list, names: list, plot_dir: str):
    """Save each figure to plot_dir/<name>.png. Names and figures are zipped."""
    os.makedirs(plot_dir, exist_ok=True)
    saved = []
    for fig, name in zip(figures, names):
        path = os.path.join(plot_dir, f"{name}.png")
        fig.savefig(path, dpi=200, bbox_inches="tight")
        saved.append(path)
        print(f"[INFO] Saved plot: {path}")
    return saved


# ── main ────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    plt.rcParams["font.family"]     = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

    # 1. Merge per-sample CSVs
    master_df = merge_qc_metrics(args.input_files, args.stage)

    # 2. Optionally join metadata
    if args.samples and os.path.exists(args.samples):
        master_df = attach_metadata_to_merged(master_df, args.samples)
    else:
        print("[INFO] No metadata file provided — skipping metadata join.")

    # 3. Save master CSV
    os.makedirs(os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True)
    master_df.to_csv(args.output_csv, index=False)
    print(f"[INFO] Master CSV saved: {args.output_csv}")

    # 4. Generate figures
    bar_figs      = fig_cohort_bar_plots(master_df)
    scatter_figs  = fig_cohort_scatters(master_df)
    violin_figs   = fig_cohort_violins(master_df)
    meta_fig      = fig_metadata_breakdown(master_df)

    # Build named list — violins and bar plots are per-metric
    bar_names     = ["bar_total_counts", "bar_n_genes", "bar_pct_mt", "bar_pct_ribo"]
    scatter_names = ["scatter_counts_vs_genes", "scatter_counts_vs_mt"]
    violin_names  = ["violin_n_genes", "violin_total_counts", "violin_pct_mt", "violin_pct_ribo"]

    all_figures = bar_figs + scatter_figs + violin_figs
    all_names   = (bar_names[:len(bar_figs)]
                   + scatter_names[:len(scatter_figs)]
                   + violin_names[:len(violin_figs)])

    if meta_fig is not None:
        all_figures.append(meta_fig)
        all_names.append("metadata_breakdown")

    # 5. Save individual PNGs
    _save_figures(all_figures, all_names, args.plot_dir)

    # 6. Build HTML report
    n_samples = master_df["sample_id"].nunique()
    html = build_html_report(
        figures=all_figures,
        sample_id=f"Cohort ({n_samples} samples)",
        stage=args.stage,
        thresholds={
            "Total samples evaluated": n_samples,
            "Total cells evaluated":   f"{len(master_df):,}",
            "Stage":                   args.stage.capitalize(),
        },
        n_cells=len(master_df),
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.output_html)), exist_ok=True)
    with open(args.output_html, "w") as f:
        f.write(html)
    print(f"[INFO] HTML report saved: {args.output_html}")

    # 7. Export PDF
    pdf_path = os.path.splitext(args.output_html)[0] + ".pdf"
    export_pdf(all_figures, pdf_path,
               sample_id=f"Cohort_{args.stage.upper()}", stage=args.stage)


if __name__ == "__main__":
    main()
