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
    get_sample_thresholds,
    compute_outlier_flags,
    resolve_display_thresholds,
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

    plt.rcParams["font.family"]     = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

    # ── Load and annotate ─────────────────────────────────────────────────────
    adata = load_data(args.input)
    adata = attach_metadata(adata, args.samples, args.sample_id)
    adata = compute_qc_metrics(adata, species=config.get("species"))

    # Per-sample thresholds from samples.csv
    sample_thr = get_sample_thresholds(args.samples, args.sample_id)

    # Compute the actual outlier flag used for filtering
    adata = compute_outlier_flags(adata, sample_thr)

    # Also write per-MAD flags for the sensitivity diagnostic plots
    adata = add_outlier_flags(adata, mad_values=[1, 2, 3, 4, 5])

    # Effective threshold values for histogram annotation lines
    display_thr = resolve_display_thresholds(adata, sample_thr)

    # Pick a highlight MAD for the sensitivity plot (first MAD setting found, or 3)
    highlight_mad = next(
        (sample_thr[k] for k in ("mad_genes_lower", "mad_counts_lower", "mad_mt")
         if sample_thr.get(k) is not None),
        3,
    )

    # ── Figures ───────────────────────────────────────────────────────────────
    figures = [
        # Standard QC plots
        fig_summary_table(adata, args.sample_id, display_thr),
        fig_violin(adata, args.sample_id),
        fig_scatter(adata, args.sample_id),
        fig_histogram_genes(adata, args.sample_id,
                            display_thr.get("min_genes"),
                            display_thr.get("max_genes")),
        fig_histogram_mt(adata, args.sample_id, display_thr.get("max_pct_mt")),
        fig_mt_decay_curve(adata, args.sample_id),
        # MAD sensitivity (diagnostic only — actual filter uses per-sample thresholds)
        fig_mad_retention_curve(adata, args.sample_id, highlight_mad=highlight_mad),
        fig_mad_sensitivity_table(adata, args.sample_id, highlight_mad=highlight_mad),
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
        thresholds=display_thr,
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