#!/usr/bin/env python3
"""
02_qc_report.py
────────────────────────────────────────────────────────────────
PRE-FILTER QC report for a single sample.

Loads raw data, computes QC metrics, writes per-MAD outlier flags
to both the AnnData and the output CSV, and produces an HTML + PDF
report that includes a MAD sensitivity section so we can
decide which MAD multiplier to use before running 03_filtering.py.

Output files (derived from --output path):
    pre_<sid>_qc.html
    pre_<sid>_qc.pdf
    pre_<sid>_cell_qc_metrics.csv   ← one row per cell, includes outlier_* columns
────────────────────────────────────────────────────────────────
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import yaml
import matplotlib.pyplot as plt
import scanpy as sc

sys.path.insert(0, os.path.dirname(__file__))
from utils.io       import load_data, attach_metadata
from utils.report   import build_html_report, export_pdf
from utils.qc_metrics import (
    compute_qc_metrics,
    thresholds_from_mad,
    add_outlier_flags,
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
from utils.mad_plots import (
    fig_mad_retention_curve,
    fig_mad_sensitivity_table,
)

sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser(description="Pre-filter QC report — single sample")
    p.add_argument("--input",     required=True, help="Count matrix (any supported format)")
    p.add_argument("--samples",  required=True, help="samples CSV (first column = sample ID)")
    p.add_argument("--cfg",       required=True, help="config.yaml")
    p.add_argument("--sample_id", required=True, help="Sample ID to look up in metadata")
    p.add_argument("--output",    required=True, help="Output HTML path")
    return p.parse_args()


def main():
    args = parse_args()

    with open(args.cfg) as f:
        config = yaml.safe_load(f)

    n_mads = config.get("qc_thresholds", {}).get("n_mads", 3)

    plt.rcParams["font.family"]     = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

    # ── Load and annotate ─────────────────────────────────────────────────────
    adata = load_data(args.input)
    adata = attach_metadata(adata, args.samples, args.sample_id)
    adata = compute_qc_metrics(adata, species=config.get("species"))

    # Write outlier flags for every MAD 1-5 into adata.obs (and later the CSV)
    adata = add_outlier_flags(adata, mad_values=[1, 2, 3, 4, 5])

    # Thresholds at the configured MAD — used for histogram annotation lines
    thresholds = thresholds_from_mad(adata, n_mads)

    # ── Figures ───────────────────────────────────────────────────────────────
    figures = [
        # Standard QC plots
        fig_summary_table(adata, args.sample_id, thresholds),
        fig_violin(adata, args.sample_id),
        fig_scatter(adata, args.sample_id),
        fig_histogram_genes(adata, args.sample_id,
                            thresholds["min_genes"], thresholds["max_genes"]),
        fig_histogram_mt(adata, args.sample_id, thresholds["max_pct_mt"]),
        fig_mt_decay_curve(adata, args.sample_id),
        # MAD 
        fig_mad_retention_curve(adata, args.sample_id, highlight_mad=n_mads),
        fig_mad_sensitivity_table(adata, args.sample_id, highlight_mad=n_mads),
    ]

    # ── CSV ───────────────────────────────────────────────────────────────────
    # save_qc_metrics auto-includes all outlier_* columns from adata.obs
    qc_df = save_qc_metrics(adata, args.sample_id)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    csv_path = args.output.replace("_qc.html", "_cell_qc_metrics.csv")
    qc_df.to_csv(csv_path, index=False)
    print(f"[INFO] Cell QC metrics written → {csv_path}")

    # ── HTML report ───────────────────────────────────────────────────────────
    html = build_html_report(
        figures=figures,
        sample_id=args.sample_id,
        stage="pre",
        thresholds=thresholds,
        n_cells=adata.n_obs,
    )
    with open(args.output, "w") as f:
        f.write(html)
    print(f"[INFO] HTML report written → {args.output}")

    # ── PDF report ────────────────────────────────────────────────────────────
    pdf_path = args.output.replace(".html", ".pdf")
    export_pdf(figures, pdf_path, args.sample_id, stage="pre")
    print(f"[INFO] PDF report written  → {pdf_path}")


if __name__ == "__main__":
    main()