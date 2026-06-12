"""
Step 10 – CellTypist Annotation
==============================================
Automated cell type annotation using CellTypist with majority voting.

Inputs
------
    h5ad : clustered AnnData from step 09

Outputs
-------
    results/10_celltypist/files/sc.h5ad
    results/10_celltypist/figures/
"""
import argparse
import os
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── CLI ────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="CellTypist annotation with majority voting")
    p.add_argument("--input",      required=True,          help="Clustered .h5ad from step 09")
    p.add_argument("--output",     required=True,          help="Output .h5ad path")
    p.add_argument("--model",      required=True,          help="CellTypist model name (e.g. Cells_Intestinal_Tract)")
    p.add_argument("--leiden_key", default=None,           help="Leiden cluster key for crosstab (e.g. leiden_0.5)")
    p.add_argument("--basis",      default="umap_harmony", help="UMAP embedding basis for plots (default: umap_harmony)")
    return p.parse_args()


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    FILES   = os.path.dirname(args.output)
    FIGURES = os.path.join(os.path.dirname(FILES), "figures")
    os.makedirs(FILES,   exist_ok=True)
    os.makedirs(FIGURES, exist_ok=True)

    sc.settings.figdir    = FIGURES
    sc.settings.verbosity = 1
    plt.rcParams["font.family"]     = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

    # ── Load ──────────────────────────────────────────────────────────────────
    adata = sc.read_h5ad(args.input)
    print(f"[INFO] Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ── Normalize for CellTypist (10k counts + log1p) ─────────────────────────
    adata.layers["X_before_celltypist"] = adata.X.copy()
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ── Load model and annotate ───────────────────────────────────────────────
    model = models.Model.load(model=f"{args.model}.pkl")
    print(f"[INFO] Loaded model: {args.model} ({len(model.cell_types)} cell types)")

    predictions = celltypist.annotate(adata, model=model, majority_voting=True)
    adata = predictions.to_adata()

    # Restore X to what it was before normalization
    adata.X = adata.layers["X_before_celltypist"].copy()
    del adata.layers["X_before_celltypist"]

    # ── UMAP: celltypist labels ───────────────────────────────────────────────

    sc.pl.embedding(
        adata,
        basis=args.basis,
        color=[key for key in ["majority_voting", "predicted_labels"]],
        save=f"_predictions.pdf",
        legend_loc="right margin",
    )

    # ── Crosstab: leiden clusters vs majority_voting ──────────────────────────
    if args.leiden_key and args.leiden_key in adata.obs:
        ax = (
            pd.crosstab(
                adata.obs[args.leiden_key],
                adata.obs["majority_voting"],
                normalize="index",
            )
            .plot.bar(stacked=True, figsize=(12, 5))
        )
        ax.legend(bbox_to_anchor=(1.05, 1), frameon=False, fontsize=6)
        plt.gcf().tight_layout()
        plt.gcf().savefig(
            os.path.join(FIGURES, f"barplot_{args.leiden_key}_vs_celltypist.pdf"),
            dpi=500,
        )
        plt.close()
    elif args.leiden_key:
        print(f"[WARNING] leiden_key '{args.leiden_key}' not found in adata.obs — skipping crosstab")

    # ── Save ──────────────────────────────────────────────────────────────────
    adata.write_h5ad(args.output)
    print(f"[INFO] Saved: {args.output}")


if __name__ == "__main__":
    main()
