#!/usr/bin/env python3
"""
03_merge_adatas.py
────────────────────────────────────────────────────────────────
Merges the per-sample AnnData objects produced by the QC/filtering
step into a single AnnData.

Usage:
    python 03_merge_adatas.py \
        --inputs results/filtered/sample1_filtered.h5ad results/filtered/sample2_filtered.h5ad \
        --join   union \
        --output results/integrated/merged.h5ad
────────────────────────────────────────────────────────────────
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import anndata as ad
import scanpy as sc

# ── Imports from utils/ ───────────────────────────────────────
# sys.path is set to workflow/scripts/ so utils/ is always found,
# regardless of which directory you launch snakemake from.
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from utils.io import load_data

sc.settings.verbosity = 0

# Map the user-facing flag to the join mode expected by anndata.concat
JOIN_MODES = {"union": "outer", "intersection": "inner"}


# ── CLI ──────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Merge per-sample AnnData objects into one")
    p.add_argument("--inputs", required=True, nargs="+",
                   help="Filtered .h5ad files to merge (one per sample)")
    p.add_argument("--join", choices=JOIN_MODES.keys(), default="union",
                   help="How to combine genes across samples: 'union' keeps every gene "
                        "seen in any sample (missing entries filled with zeros), "
                        "'intersection' keeps only genes shared by all samples. "
                        "Default: union")
    p.add_argument("--output", required=True, help="Output merged .h5ad path")
    return p.parse_args()


# ── Main ─────────────────────────────────────────────────────

def main():
    args = parse_args()
    join_mode = JOIN_MODES[args.join]

    print(f"[INFO] {'='*55}")
    print(f"[INFO] Merging {len(args.inputs)} samples | join: {args.join} ({join_mode})")
    print(f"[INFO] {'='*55}")

    adatas = []
    for path in args.inputs:
        print(f"[INFO] Loading {path}")
        adata = load_data(path)
        print(f"[INFO]   {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        adatas.append(adata)

    print(f"[INFO] Concatenating with join='{join_mode}'...")
    merged = ad.concat(adatas, join=join_mode, index_unique="-")
    print(f"[INFO] Merged: {merged.n_obs:,} cells x {merged.n_vars:,} genes")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    merged.write_h5ad(args.output)
    print(f"[INFO] Saved: {args.output}")
    print("[INFO] Done.")


if __name__ == "__main__":
    main()
