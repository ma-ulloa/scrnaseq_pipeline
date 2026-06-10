"""
Step 07 – Clustering
==============================================
Leiden clustering at multiple resolutions followed by marker gene scoring
to support manual cell type annotation.

Inputs
------
    h5ad  :  sc.h5ad from step 06 (Harmony-integrated)

Outputs
-------
    results/07_clustering/files/sc.h5ad
    results/07_clustering/figures/  (UMAPs, dotplots, marker plots)
"""
import argparse
import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


# ── CLI ────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Leiden clustering at multiple resolutions + marker gene scoring")
    p.add_argument("--input",       required=True,              help="Harmony-integrated .h5ad file")
    p.add_argument("--markers",     required=True,              help="CSV with columns: gene, cell_state")
    p.add_argument("--output",      required=True,              help="Output .h5ad path")
    p.add_argument("--resolutions", required=True, nargs="+", type=float,
                   help="Leiden resolutions to test, e.g. 0.3 0.5 0.8")
    p.add_argument("--seed",        type=int, default=42,          help="Random seed (default: 42)")
    p.add_argument("--basis",       default="umap_harmony",        help="Embedding basis for UMAP plots (default: umap_harmony)")
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

    # ── Load ───────────────────────────────────────────────────────────────
    adata = sc.read_h5ad(args.input)

    # ── Leiden clustering ──────────────────────────────────────────────────
    for res in args.resolutions:
        sc.tl.leiden(adata,
                     resolution=res,
                     key_added=f"leiden_{res}",
                     random_state=args.seed,
                     flavor="igraph",
                     n_iterations=-1)

    sc.pl.embedding(
        adata,
        basis=args.basis,
        color=[f"leiden_{res}" for res in args.resolutions],
        ncols=3,
        save="_clusters.pdf",
        legend_loc="on data",
    )

    # ── Marker genes per cluster ───────────────────────────────────────────
    for res in args.resolutions:
        sc.tl.rank_genes_groups(
            adata,
            groupby=f"leiden_{res}",
            method="wilcoxon",
            use_raw=False,
        )
        sc.pl.rank_genes_groups(adata, n_genes=20, save=f"_marker_genes_{res}.pdf")
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, save=f"_marker_dotplot_{res}.pdf")

    # ── Plot known markers in UMAP ─────────────────────────────────────────
    markers = pd.read_csv(args.markers)
    for cell_state, group in markers.groupby("cell_state"):
        present = [g for g in group["gene"] if g in adata.var_names]
        if not present:
            print(f"Warning: no markers found for {cell_state}, skipping")
            continue
        sc.pl.embedding(
            adata,
            basis=args.basis,
            color=present,
            ncols=3,
            save=f"_{cell_state}.pdf",
            cmap="Oranges",
        )

    # ── Save ───────────────────────────────────────────────────────────────
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()
