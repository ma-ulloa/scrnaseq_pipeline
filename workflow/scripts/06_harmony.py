#!/usr/bin/env python3
"""
06_harmony.py
────────────────────────────────────────────────────────────────
Pre-processes the merged dataset (normalization + optional scaling),
computes PCA, clusters and embeds it with UMAP, plots the result
coloured by sample and cluster, integrates samples with Harmony,
and repeats clustering/UMAP on the Harmony embedding for comparison.

Usage:
    python 05_harmony.py \
        --input            results/03_merge_adatas/files/sc.h5ad \
        --config           config/config.yaml \
        --integration_keys sample_id fraction \
        --output           results/05_harmony/sc.h5ad \
        --plots            results/05_harmony/sc_harmony.pdf
────────────────────────────────────────────────────────────────
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import harmonypy as hm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc
import yaml

sc.settings.verbosity = 0


# ── CLI ──────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Pre-process, cluster and integrate samples with Harmony")
    p.add_argument("--input",            required=True, help="Merged .h5ad file")
    p.add_argument("--config",           required=True, help="config.yaml")
    p.add_argument("--integration_keys", required=True, nargs="+",
                   help="adata.obs column(s) used as Harmony batch keys, e.g. sample_id fraction")
    p.add_argument("--n_pcs",  type=int, default=20, help="Number of leading PCs handed to Harmony (default: 20)")
    p.add_argument("--output",           required=True, help="Output .h5ad path")
    p.add_argument("--plots",            required=True, help="Output PDF path for PCA/UMAP plots")
    return p.parse_args()


# ── Pre-processing ───────────────────────────────────────────

def normalize(adata: sc.AnnData, method: str, factor: float) -> sc.AnnData:
    """Log-normalize counts using the method given in config['pp']['normalize']."""
    if method == "total_counts":
        sc.pp.normalize_total(adata, target_sum=factor)
    elif method == "median":
        sc.pp.normalize_total(adata)
    else:
        sys.exit(f"[ERROR] Unsupported normalization method: '{method}'")
    sc.pp.log1p(adata)
    print(f"[INFO] Normalized (method='{method}') and log1p-transformed")
    return adata


# ── Clustering / UMAP ────────────────────────────────────────

def cluster_and_embed(adata: sc.AnnData, use_rep: str, cluster_key: str, umap_key: str) -> sc.AnnData:
    """Compute neighbors, Leiden clusters and a UMAP from a given embedding."""
    sc.pp.neighbors(adata, use_rep=use_rep)
    sc.tl.leiden(adata, key_added=cluster_key)
    sc.tl.umap(adata)
    adata.obsm[umap_key] = adata.obsm["X_umap"].copy()
    print(f"[INFO] Clustering on '{use_rep}': {adata.obs[cluster_key].nunique()} clusters")
    return adata


# ── Plotting ─────────────────────────────────────────────────

def fig_embedding(adata: sc.AnnData, basis: str, color: str, title: str) -> plt.Figure:
    ax = sc.pl.embedding(adata, basis=basis, color=color, show=False, title=title)
    fig = ax.get_figure()
    plt.tight_layout()
    return fig


def fig_pca_variance_ratio(adata: sc.AnnData, n_pcs: int) -> plt.Figure:
    sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, log=True, show=False)
    fig = plt.gcf()
    plt.tight_layout()
    return fig


def save_plots(figures: list, pdf_path: str):
    os.makedirs(os.path.dirname(os.path.abspath(pdf_path)), exist_ok=True)
    with PdfPages(pdf_path) as pdf:
        for fig in figures:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"[INFO] Plots saved to: {pdf_path}")


# ── Main ─────────────────────────────────────────────────────

def main():
    args = parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    pp_config = config["pp"]

    print(f"[INFO] {'='*55}")
    print(f"[INFO] Harmony integration | keys: {args.integration_keys}")
    print(f"[INFO] {'='*55}")

    print(f"[INFO] Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[INFO] {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Pre-processing
    adata = normalize(adata, pp_config["normalize"], pp_config["factor"])
    if pp_config["scale"]:
        print("[INFO] Scaling...")
        sc.pp.scale(adata)
    else:
        print("[INFO] Scaling disabled (config['pp']['scale'] = false)")

    # PCA, clustering and UMAP — before integration
    print("[INFO] Computing PCA...")
    sc.tl.pca(adata)
    adata = cluster_and_embed(adata, use_rep="X_pca", cluster_key="leiden", umap_key="X_umap_pca")

    figures = [
        fig_pca_variance_ratio(adata, n_pcs=args.n_pcs),
        fig_embedding(adata, "pca",      "sample_id", "PCA — coloured by sample"),
        fig_embedding(adata, "pca",      "leiden",    "PCA — coloured by cluster"),
        fig_embedding(adata, "umap_pca", "sample_id", "UMAP (PCA) — coloured by sample"),
        fig_embedding(adata, "umap_pca", "leiden",    "UMAP (PCA) — coloured by cluster"),
    ]

    # Harmony integration — run harmonypy directly on the leading PCs,
    # mirroring the approach used in spatial/XeniumKKIK/05integration_instanseg
    thetas = config["harmony"].get("thetas")
    if thetas is not None and len(thetas) != len(args.integration_keys):
        sys.exit(
            f"[ERROR] harmony.thetas has {len(thetas)} value(s) but "
            f"{len(args.integration_keys)} integration key(s) were provided — must match."
        )

    print(f"[INFO] Running Harmony integration on keys: {args.integration_keys} "
          f"(n_pcs={args.n_pcs}, thetas={thetas})")
    pcs         = adata.obsm["X_pca"][:, :args.n_pcs]
    harmony_res = hm.run_harmony(pcs, adata.obs, args.integration_keys, theta=thetas)
    adata.obsm["X_pca_harmony"] = harmony_res.Z_corr.T

    # Clustering and UMAP — repeated on the Harmony embedding
    adata = cluster_and_embed(adata, use_rep="X_pca_harmony", cluster_key="leiden_harmony", umap_key="X_umap_harmony")

    figures += [
        fig_embedding(adata, "pca_harmony",  "sample_id",      "Harmony PCA — coloured by sample"),
        fig_embedding(adata, "pca_harmony",  "leiden_harmony", "Harmony PCA — coloured by cluster"),
        fig_embedding(adata, "umap_harmony", "sample_id",      "UMAP (Harmony) — coloured by sample"),
        fig_embedding(adata, "umap_harmony", "leiden_harmony", "UMAP (Harmony) — coloured by cluster"),
    ]

    save_plots(figures, args.plots)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"[INFO] Saved: {args.output}")
    print("[INFO] Done.")


if __name__ == "__main__":
    main()
