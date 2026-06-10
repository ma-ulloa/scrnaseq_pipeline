#!/usr/bin/env python3
"""
04_qc_post_report.py
────────────────────────────────────────────────────────────────
POST-FILTER QC report for a single sample.

Usage:
    python 04_qc_post_report.py \
        --input     results/filtered/KS2204T1_2_filtered.h5ad \
        --metadata  config/metadata.csv \
        --cfg       config/config.yaml \
        --sample_id KS2204T1_2 \
        --output    results/qc/post_KS2204T1_2_qc.html

Output files (derived from --output path):
    post_<sid>_qc.html
    post_<sid>_qc.pdf
    post_<sid>_cell_qc_metrics.csv
────────────────────────────────────────────────────────────────
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import yaml
import scanpy as sc

sys.path.insert(0, os.path.dirname(__file__))
from utils.io       import load_data, attach_metadata
from utils.report   import build_html_report, export_pdf
from utils.qc_metrics import (
    compute_qc_metrics,
    thresholds_from_mad,
    save_qc_metrics,
)
from utils.qc_plots import (
    fig_summary_table,
    fig_violin,
    fig_scatter,
    fig_histogram_genes,
    fig_histogram_mt,
    fig_mt_decay_curve,
)

sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser(description="Post-filter QC report — single sample")
    p.add_argument("--input",     required=True, help="Filtered .h5ad from 03_filtering.py")
    p.add_argument("--metadata",  required=True, help="Metadata CSV (first column = sample ID)")
    p.add_argument("--cfg",       required=True, help="config.yaml")
    p.add_argument("--sample_id", required=True, help="Sample ID to look up in metadata")
    p.add_argument("--output",    required=True, help="Output HTML path")
    return p.parse_args()


def main():
    args = parse_args()

    with open(args.cfg) as f:
        config = yaml.safe_load(f)

    n_mads = config.get("qc_thresholds", {}).get("n_mads", 3)

    # ── Load filtered data ────────────────────────────────────────────────────
    adata = load_data(args.input)
    adata = attach_metadata(adata, args.metadata, args.sample_id)

    # Recompute QC metrics on the filtered cell population.
    # Gene/count distributions shift post-filter, so fresh metrics are needed.
    adata = compute_qc_metrics(adata)

    # Thresholds are shown as reference lines on the histogram plots.
    # Post-filter they describe what was applied, not what will be applied.
    thresholds = thresholds_from_mad(adata, n_mads)

    # ── Figures (standard plots only — no MAD sensitivity) ────────────────────
    figures = [
        fig_summary_table(adata, args.sample_id, thresholds),
        fig_violin(adata, args.sample_id),
        fig_scatter(adata, args.sample_id),
        fig_histogram_genes(adata, args.sample_id,
                            thresholds["min_genes"], thresholds["max_genes"]),
        fig_histogram_mt(adata, args.sample_id, thresholds["max_pct_mt"]),
    ]

    # ── CSV ───────────────────────────────────────────────────────────────────
    # No outlier flags added here — the post-filter CSV is a clean record
    # of surviving cells and their QC values.
    qc_df = save_qc_metrics(adata, args.sample_id)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    csv_path = args.output.replace("_qc.html", "_cell_qc_metrics.csv")
    qc_df.to_csv(csv_path, index=False)
    print(f"[INFO] Post-filter cell QC metrics written → {csv_path}")

    # ── HTML report ───────────────────────────────────────────────────────────
    html = build_html_report(
        figures=figures,
        sample_id=args.sample_id,
        stage="post",
        thresholds=thresholds,
        n_cells=adata.n_obs,
    )
    with open(args.output, "w") as f:
        f.write(html)
    print(f"[INFO] HTML report written → {args.output}")

    # ── PDF report ────────────────────────────────────────────────────────────
    pdf_path = args.output.replace(".html", ".pdf")
    export_pdf(figures, pdf_path, args.sample_id, stage="post")
    print(f"[INFO] PDF report written  → {pdf_path}")


if __name__ == "__main__":
    main()