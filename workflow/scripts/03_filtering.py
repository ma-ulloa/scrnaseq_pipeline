#!/usr/bin/env python3
import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import scanpy as sc
import yaml

sys.path.insert(0, os.path.dirname(__file__))
from utils.io import load_data, attach_metadata

def parse_args():
    p = argparse.ArgumentParser(description="Filter cells by pre-computed MAD thresholds & optional doublets")
    p.add_argument("--input",     required=True, help="Raw count matrix (.h5ad or .h5)")
    p.add_argument("--metadata",  required=True, help="Metadata CSV")
    p.add_argument("--qc_csv",    required=True, help="Pre-computed cell_qc_metrics.csv containing MAD flags")
    p.add_argument("--config",    required=True, help="config.yaml")
    p.add_argument("--sample_id", required=True, help="Sample ID")
    p.add_argument("--output",    required=True, help="Output filtered .h5ad path")
    return p.parse_args()

def main():
    args = parse_args()
    
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    min_cells_per_gene = config["qc_thresholds"].get("min_cells_per_gene", 3)
    
    # Define if scrublet needs to be run
    doublet_config = config.get("doublets", {})
    run_doublet_filter = doublet_config.get("run_scrublet", True)
    expected_doublet_rate = doublet_config.get("expected_rate", 0.08)

    print(f"[INFO] Processing Sample: {args.sample_id}")
    adata = load_data(args.input)
    adata = attach_metadata(adata, args.metadata, args.sample_id)
    
    # Save raw counts matrix to a layer for downstream differential expression
    adata.layers["counts"] = adata.X.copy()
    n_cells_start = adata.n_obs

    # Apply MAD filters
    qc_metrics = pd.read_csv(args.qc_csv)
    pass_mad = qc_metrics["pass_mad_filters"].values
    adata = adata[pass_mad].copy()
    print(f"[INFO] Removed {n_cells_start - adata.n_obs:,} cells failing adaptive MAD thresholds.")

    # Optional double detection
    if run_doublet_filter:
        print("[INFO] Running Scrublet...")
        sc.pp.scrublet(adata, expected_doublet_rate=expected_doublet_rate)
        
        is_singlet = adata.obs["predicted_doublet"] == False
        n_doublets = (~is_singlet).sum()
        adata = adata[is_singlet].copy()
        print(f"[INFO] Removed {n_doublets:,} predicted doublets.")
    else:
        print("[INFO] Skipping Scrublet doublet detection (disabled in config).")

    # 3. remove low representation genes
    n_genes_start = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    print(f"[INFO] Removed {n_genes_start - adata.n_vars:,} noise genes present in < {min_cells_per_gene} cells.")

    # Final pass stats summary
    pct_kept = (adata.n_obs / n_cells_start) * 100
    print(f"[INFO] Final Clean Dataset: {adata.n_obs:,} / {n_cells_start:,} cells kept ({pct_kept:.1f}%)")

    # Save output
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"[INFO] Saved clean dataset to: {args.output}\n")

if __name__ == "__main__":
    main()