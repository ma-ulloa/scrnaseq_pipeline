#!/usr/bin/env python3
"""
01_qc_report.py
────────────────────────────────────────────────────────────────
QC report orchestrator for the scRNA-seq pipeline.

Loads data, computes QC metrics, builds figures, and writes
an HTML + PDF report. Designed to run pre- and post-filtering.

Usage:
    python 01_qc_report.py \
        --input     data/raw/sample.h5 \
        --metadata  config/metadata.csv \
        --config    config/config.yaml \
        --sample_id sample \
        --stage     pre \
        --output    results/qc/pre_sample_qc.html

Input format is auto-detected from the file extension:
    .h5ad → AnnData  |  .h5 → 10x or generic HDF5
    .mtx / dir → 10x MEX  |  .csv / .txt / .tsv → dense matrix
────────────────────────────────────────────────────────────────
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import scanpy as sc
import yaml
import numpy as np

# Pipeline modules
sys.path.insert(0, os.path.dirname(__file__))
from utils.io    import load_data, attach_metadata
from utils.report import build_html_report, export_pdf
from utils.qc_plots import (
    fig_summary_table,
    fig_violin,
    fig_scatter,
    fig_histogram_genes,
    fig_histogram_mt,
    fig_multi_mad_histograms,
    fig_mt_decay_curve,
)

sc.settings.verbosity = 1


# ── Mitochondrial / Ribosomal prefixes ───────────────────────
# Human: MT-  |  Mouse: mt-  |  other conventions: MT_
MT_PREFIXES   = ("MT-", "mt-", "MT_")
RIBO_PREFIXES = ("RPS", "RPL", "Rps", "Rpl")


# Arguments to parse file

def parse_args():
    p = argparse.ArgumentParser(description="Generate scRNA-seq QC report")
    p.add_argument("--input",     required=True, help="Count matrix (any supported format)")
    p.add_argument("--metadata",  required=True, help="Metadata CSV (first column = sample ID)")
    p.add_argument("--cfg",       required=True, help="config.yaml")
    p.add_argument("--sample_id", required=True, help="Sample ID to look up in metadata")
    p.add_argument("--stage",     required=True, choices=["pre", "post"],
                   help="'pre' = before filtering, 'post' = after filtering")
    p.add_argument("--output",    required=True, help="Output HTML path (PDF saved alongside)")
    return p.parse_args()


# QC Metrics 

def compute_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    """Flag MT / ribo genes and compute per-cell QC metrics."""
    adata.var["mt"]   = adata.var_names.str.startswith(MT_PREFIXES)
    adata.var["ribo"] = adata.var_names.str.startswith(RIBO_PREFIXES)

    n_mt   = adata.var["mt"].sum()
    n_ribo = adata.var["ribo"].sum()
    print(f"[INFO] MT genes: {n_mt} | Ribo genes: {n_ribo}")
    if n_mt == 0:
        print("[WARN] No mitochondrial genes detected — check gene naming convention.")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"],
        percent_top=None, log1p=True, inplace=True,
    )
    print(f"[INFO] QC metrics computed. Cells: {adata.n_obs} | Genes: {adata.n_vars}")
    return adata


def calculate_mad_thresholds(series: pd.Series, n_mads: float = 3.0) -> tuple:
    """
    Thresholds usingf the Median Absolute Deviation(a.k.a MAD)
    """
    median = series.median()
    mad = np.median(np.abs(series - median))
    
    scaled_mad = mad * 1.4826
    
    lower_limit = median - (n_mads * scaled_mad)
    upper_limit = median + (n_mads * scaled_mad)
    
    return max(0, lower_limit), upper_limit

# QC summary tables
def save_qc_metrics(adata: sc.AnnData, sample_id: str) -> pd.DataFrame:
    """
    Extracts the computed QC metrics from adata.obs into a standalone pandas DataFrame.
    """
    qc_columns = [
        "n_genes_by_counts",
        "log1p_n_genes_by_counts",
        "total_counts",
        "log1p_total_counts",
        "pct_counts_mt",
        "pct_counts_ribo"
    ]
    
    available_columns = [col for col in qc_columns if col in adata.obs.columns]
    metrics_df = adata.obs[available_columns].copy()
    metrics_df.insert(0, "sample_id", sample_id)

    return metrics_df

#  Main 

def main():
    args = parse_args()

    with open(args.cfg) as f:
        config = yaml.safe_load(f)
    
    # Let's see if n_mads is configured, default to 3 if not
    n_mads = config.get("qc_thresholds", {}).get("n_mads", 3)

    # Load + annotate
    adata = load_data(args.input)
    adata = attach_metadata(adata, args.metadata, args.sample_id)
    adata = compute_qc_metrics(adata)

    log_genes = adata.obs["log1p_n_genes_by_counts"]
    log_counts = adata.obs["log1p_total_counts"]
    log_gene_low, log_gene_high = calculate_mad_thresholds(log_genes, n_mads=n_mads)
    log_count_low, _ = calculate_mad_thresholds(log_counts, n_mads=n_mads)
    
    # Convert log-space thresholds back to natural linear scale space f
    mad_thresholds = {
        "min_genes": int(np.expm1(log_gene_low)),
        "max_genes": int(np.expm1(log_gene_high)),
        "min_counts": int(np.expm1(log_count_low)),
        "max_pct_mt": float(calculate_mad_thresholds(adata.obs["pct_counts_mt"], n_mads=n_mads)[1])
    }
    
    # Figures
    figures = [
        fig_summary_table(adata, args.sample_id, mad_thresholds),
        fig_violin(adata, args.sample_id),
        fig_scatter(adata, args.sample_id),
        fig_histogram_genes(adata, args.sample_id,
                            mad_thresholds["min_genes"], mad_thresholds["max_genes"]),
        fig_histogram_mt(adata, args.sample_id, mad_thresholds["max_pct_mt"]),
        fig_multi_mad_histograms(adata, args.sample_id),
        fig_mt_decay_curve(adata, args.sample_id),
    ]
    
    # Save qc metrics summary table
    qc_df = save_qc_metrics(adata, args.sample_id)
    
    # Label cell-specific pass/fail flags directly inside your output CSV for downstream tracking
    qc_df["pass_mad_filters"] = (
        (adata.obs["n_genes_by_counts"] >= mad_thresholds["min_genes"]) &
        (adata.obs["n_genes_by_counts"] <= mad_thresholds["max_genes"]) &
        (adata.obs["total_counts"] >= mad_thresholds["min_counts"]) &
        (adata.obs["pct_counts_mt"] <= mad_thresholds["max_pct_mt"])
    )

    # HTML
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    html = build_html_report(
        figures=figures,
        sample_id=args.sample_id,
        stage=args.stage,
        thresholds=mad_thresholds, 
        n_cells=adata.n_obs,
    )
    with open(args.output, "w") as f:
        f.write(html)

    # PDF
    pdf_path = os.path.splitext(args.output)[0] + ".pdf"
    export_pdf(figures, pdf_path, args.sample_id, args.stage)

    # CSV
    csv_path = args.output.replace("_qc.html", "_cell_qc_metrics.csv")
    qc_df.to_csv(csv_path, index=False)

if __name__ == "__main__":
    main()