#!/usr/bin/env python3
"""
01_qc_report.py
────────────────────────────────────────────────────────────────
QC report orchestrator for the scRNA-seq pipeline.

Loads data, computes QC metrics, builds figures, and writes
an HTML + PDF report. Designed to run pre- and post-filtering.

Usage:
    python 01_qc_report.py \
        --input     data/raw/KS2103T1_1.h5 \
        --metadata  config/metadata.csv \
        --config    config/config.yaml \
        --sample_id KS2103T1_1 \
        --stage     pre \
        --output    results/qc/pre_KS2103T1_1_qc.html

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


# QC summary tables
def save_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    """
    creates csv dataframe with general metrics 
    """
    data = 
    metrics_df = pd.DataFrame(data)


#  Main 

def main():
    args = parse_args()

    # Config
    with open(args.cfg) as f:
        config = yaml.safe_load(f)
    thresholds = config["qc_thresholds"]

    # Load + annotate
    adata = load_data(args.input)
    adata = attach_metadata(adata, args.metadata, args.sample_id)
    adata = compute_qc_metrics(adata)

    # Build figures
    figures = [
        fig_summary_table(adata, args.sample_id, thresholds),
        fig_violin(adata, args.sample_id),
        fig_scatter(adata, args.sample_id),
        fig_histogram_genes(adata, args.sample_id,
                            thresholds["min_genes"], thresholds["max_genes"]),
        fig_histogram_mt(adata, args.sample_id, thresholds["max_pct_mt"]),
    ]
    
    # Build qc metrics summary table



    # HTML
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    html = build_html_report(
        figures=figures,
        sample_id=args.sample_id,
        stage=args.stage,
        thresholds=thresholds,
        n_cells=adata.n_obs,
    )
    with open(args.output, "w") as f:
        f.write(html)
    print(f"[INFO] HTML report saved to: {args.output}")

    # PDF (same base path, .pdf extension)
    pdf_path = os.path.splitext(args.output)[0] + ".pdf"
    export_pdf(figures, pdf_path, args.sample_id, args.stage)


if __name__ == "__main__":
    main()