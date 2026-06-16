"""
Step 07 – FACS Contamination Removal
==============================================
Because samples are FACS-sorted, there is cross-fraction contamination:
  - Epcam+ (epithelial) cells present in the Immune fraction
  - Ptprc+ (CD45, immune) cells present in the EPI+DN fraction

Strategy:
  1. Concatenate QC-passed per-sample AnnData objects
  2. Remove cells expressing Epcam in the Immune fraction
  3. Remove cells expressing Ptprc in the EPI+DN fraction
  4. Normalize, PCA, UMAP, low-resolution Leiden clustering
  5. For each cluster, remove minority-fraction cells
     (cells from "Non-sorted" are always kept)

Inputs
------
    h5ads : per-sample QC-filtered AnnData files

Outputs
-------
    results/07_remove_facs_contamination/files/sc.h5ad
    results/07_remove_facs_contamination/figures/
"""
import argparse
import os
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


# ── CLI ────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Remove FACS cross-fraction contamination")
    p.add_argument("--input",     required=True, help="Merged .h5ad file")
    p.add_argument("--output",     required=True,            help="Output .h5ad path")
    p.add_argument("--seed",       type=int, default=42,     help="Random seed (default: 42)")
    p.add_argument("--resolution", type=float, default=1.5,  help="Leiden resolution for contamination clustering (default: 1.5)")
    p.add_argument("--species" , required=True,type=str,help="Species the data was obtained from"),
    p.add_argument("--non_sorted_flag", required=True, help="How cells coming from non-FACS sorted samples are specified in the sample file")
    p.add_argument("--immune_label",    required=True, help="Value of 'fraction' column for the immune (CD45+) sort")
    p.add_argument("--epidn_label",     required=True, help="Value of 'fraction' column for the epithelial/double-negative sort")
    p.add_argument("--remove_clusters", nargs="*", default=[], help="leiden_contam cluster IDs to drop entirely (e.g. 3 7)")
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

    # ── Concatenate samples ────────────────────────────────────────────────
    adata  = sc.read_h5ad(args.input)
    print(f"[INFO] Concatenated: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ── Step 1: Remove marker-positive cells from wrong fraction ───────────
    
    if args.species == "human":
        immune_marker = "PTPRC"
        epi_marker    = "EPCAM"
    elif args.species == "mouse":
        immune_marker = "Ptprc"
        epi_marker    = "Epcam"        
    
    epcam_in_immune = (
        (adata.obs["fraction"] == args.immune_label) &
        (adata[:, epi_marker].X.toarray().ravel() > 0)
    )
    ptprc_in_epidn = (
        (adata.obs["fraction"] == args.epidn_label) &
        (adata[:, immune_marker].X.toarray().ravel() > 0)
    )

    print(f"[INFO] Removing {epcam_in_immune.sum():,} {epi_marker}+ cells from {args.immune_label} fraction")
    print(f"[INFO] Removing {ptprc_in_epidn.sum():,} {immune_marker}+ cells from {args.epidn_label} fraction")

    adata = adata[~(epcam_in_immune | ptprc_in_epidn)].copy()

    # ── Normalize and embed ────────────────────────────────────────────────
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm_counts"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    sc.tl.pca(adata, n_comps=15, mask_var="highly_variable", random_state=args.seed)
    sc.pp.neighbors(adata, n_neighbors=40, n_pcs=15, random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)

    sc.pl.umap(
        adata,
        color=[epi_marker, immune_marker, "fraction", "sample_id"],
        ncols=2,
        save="_markers_and_fraction.pdf",
    )

    # ── Step 2: Cluster to catch remaining contamination ───────────────────
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        random_state=args.seed,
        key_added="leiden_contam",
        flavor="igraph",
        n_iterations=-1,
    )

    sc.pl.umap(
        adata,
        color=["leiden_contam", "fraction"],
        ncols=2,
        save="_clusters_fraction.pdf",
    )

    sc.tl.rank_genes_groups(adata, groupby="leiden_contam", use_raw=False)
    sc.pl.rank_genes_groups(adata, n_genes=20, save="_leiden_contam.pdf")

    # ── Step 3: Drop user-specified clusters entirely ──────────────────────
    if args.remove_clusters:
        to_drop = [str(c) for c in args.remove_clusters]
        mask = adata.obs["leiden_contam"].isin(to_drop)
        print(f"[INFO] Dropping clusters {to_drop}: {mask.sum():,} cells removed")
        adata = adata[~mask].copy()
        print(f"[INFO] Retained: {adata.n_obs:,} cells")

    # ── Step 4: Remove minority-fraction cells per cluster ─────────────────
    # Non-sorted cells are always kept regardless of cluster composition.
    sorted_mask = adata.obs["fraction"] != args.non_sorted_flag

    if not sorted_mask.any():
        print("[INFO] No FACS-sorted cells found — skipping minority-fraction removal")
    else:
        # Crosstab computed only on sorted cells to find the dominant fraction per cluster
        crosstab = pd.crosstab(
            adata.obs.loc[sorted_mask, "leiden_contam"],
            adata.obs.loc[sorted_mask, "fraction"],
        )
        crosstab.to_csv(os.path.join(FILES, "cluster_fraction_table.csv"))
        print(crosstab)

        dominant      = crosstab.idxmax(axis=1)                    # cluster → dominant fraction
        cell_dominant = adata.obs["leiden_contam"].map(dominant)   # NaN for clusters with no sorted cells

        # Minority = sorted cell whose fraction differs from its cluster's dominant
        is_minority = (
            sorted_mask &
            cell_dominant.notna() &
            (adata.obs["fraction"] != cell_dominant)
        )
        print(f"[INFO] Removing {is_minority.sum():,} minority-fraction cells")
        adata = adata[~is_minority].copy()
        print(f"[INFO] Retained: {adata.n_obs:,} cells")

    # ── Save ───────────────────────────────────────────────────────────────
    adata.write_h5ad(args.output)
    print(f"[INFO] Saved: {args.output}")


if __name__ == "__main__":
    main()
