#!/usr/bin/env python3
"""
02_filtering.py
────────────────────────────────────────────────────────────────
Applies QC thresholds from config.yaml to a single sample,
optionally runs Scrublet for doublet detection,
and saves the filtered AnnData as .h5ad.

Usage:
    python qc_filtering.py \
        --input     data/raw/matrix.h5 \
        --metadata  config/metadata.csv \
        --config    config/config.yaml \
        --sample_id sample \
        --output    results/filtered/sample_filtered.h5ad
────────────────────────────────────────────────────────────────
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import yaml

# ── Imports from utils/ ───────────────────────────────────────
# sys.path is set to workflow/scripts/ so utils/ is always found,
# regardless of which directory you launch snakemake from.
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from utils.io import load_data, attach_metadata

sc.settings.verbosity = 0

MT_PREFIXES   = ("MT-", "mt-", "MT_")
RIBO_PREFIXES = ("RPS", "RPL", "Rps", "Rpl")


# ── CLI ──────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Filter cells by QC thresholds")
    p.add_argument("--input",     required=True,  help="Count matrix (any supported format)")
    p.add_argument("--metadata",  required=True,  help="Metadata CSV")
    p.add_argument("--config",    required=True,  help="config.yaml")
    p.add_argument("--sample_id", required=True,  help="Sample ID")
    p.add_argument("--output",    required=True,  help="Output .h5ad path")
    return p.parse_args()


#  QC metrics 

def compute_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    adata.var["mt"]   = adata.var_names.str.startswith(MT_PREFIXES)
    adata.var["ribo"] = adata.var_names.str.startswith(RIBO_PREFIXES)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"],
        percent_top=None, log1p=False, inplace=True,
    )
    print(f"[INFO] Before filtering | Cells: {adata.n_obs:,} | Genes: {adata.n_vars:,}")
    return adata


# Filtering

def filter_cells(adata: sc.AnnData, thresholds: dict) -> sc.AnnData:
    """
    Apply QC thresholds from config.yaml.

    Order matters:
      1. Remove cells with too few genes     → empty droplets
      2. Remove cells with too many genes    → likely doublets
      3. Remove cells with too few counts    → low quality
      4. Remove cells with too high % MT     → dying cells
      5. Remove genes in too few cells       → noise reduction
    """
    n_start = adata.n_obs

    # Too few genes
    mask    = adata.obs["n_genes_by_counts"] >= thresholds["min_genes"]
    adata   = adata[mask].copy()
    print(f"[INFO] Removed {(~mask).sum():,} cells with < {thresholds['min_genes']} genes")

    # Too many genes
    mask    = adata.obs["n_genes_by_counts"] <= thresholds["max_genes"]
    adata   = adata[mask].copy()
    print(f"[INFO] Removed {(~mask).sum():,} cells with > {thresholds['max_genes']} genes")

    # Too few counts
    mask    = adata.obs["total_counts"] >= thresholds["min_counts"]
    adata   = adata[mask].copy()
    print(f"[INFO] Removed {(~mask).sum():,} cells with < {thresholds['min_counts']} counts")

    # High % MT
    mask    = adata.obs["pct_counts_mt"] <= thresholds["max_pct_mt"]
    adata   = adata[mask].copy()
    print(f"[INFO] Removed {(~mask).sum():,} cells with > {thresholds['max_pct_mt']}% MT")

    # Gene filter
    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=thresholds["min_cells_per_gene"])
    print(f"[INFO] Removed {n_genes_before - adata.n_vars:,} genes in < {thresholds['min_cells_per_gene']} cells")

    # Summary
    pct_kept = (adata.n_obs / n_start) * 100
    print(f"[INFO] Kept {adata.n_obs:,} / {n_start:,} cells ({pct_kept:.1f}%)")

    return adata


#  Doublet detection 

def run_scrublet(adata: sc.AnnData, expected_rate: float) -> sc.AnnData:
    """
    Flag probable doublets using Scrublet.
    """

    print(f"[INFO] Running Scrublet")
    # 'doublet_score' and 'predicted_doublet' to adata.obs
    sc.pp.scrublet(adata, expected_doublet_rate=expected_rate)

    n   = int(adata.obs["predicted_doublet"].sum())
    pct = (n / adata.n_obs) * 100
    print(f"[INFO] Predicted doublets: {n:,} ({pct:.1f}%) — flagged, not removed")
    return adata


#  Main 

def main():
    args = parse_args()

    print(f"[INFO] {'='*55}")
    print(f"[INFO] Filtering | sample: {args.sample_id}")
    print(f"[INFO] {'='*55}")

    with open(args.config) as f:
        config = yaml.safe_load(f)
    thresholds     = config["qc_thresholds"]
    doublet_config = config.get("doublets", {})
    print(f"[INFO] Thresholds: {thresholds}")

    print("[INFO] Loading data...")
    adata = load_data(args.input)

    print("[INFO] Attaching metadata...")
    adata = attach_metadata(adata, args.metadata, args.sample_id)

    print("[INFO] Computing QC metrics...")
    adata = compute_qc_metrics(adata)

    # Save raw counts before any modification — needed for DEG analysis later
    adata.layers["counts"] = adata.X.copy()
    print("[INFO] Raw counts saved to adata.layers['counts']")

    print("[INFO] Applying filters...")
    adata = filter_cells(adata, thresholds)

    if doublet_config.get("run_scrublet", True):
        adata = run_scrublet(adata)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"[INFO] Saved: {args.output}")
    print("[INFO] Done.")


if __name__ == "__main__":
    main()