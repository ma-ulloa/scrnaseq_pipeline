"""
qc_plots.py

All QC plot functions 
"""

import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import pandas as pd

# Set a clean, publication-ready style globally
sns.set_theme(style="whitegrid")
plt.rcParams['figure.titlesize'] = 14
plt.rcParams['axes.titlesize'] = 12

# Summary Table (Rendered as a clean visual table)
def fig_summary_table(adata: sc.AnnData, sample_id: str, thresholds: dict) -> plt.Figure:
    obs = adata.obs
    stats = {
        "Metric": [
            "Total cells", "Total genes", "Median genes / cell", 
            "Median counts / cell", "Median % MT", "Median % Ribo",
            "Cells > max_genes", "Cells < min_genes", "Cells > max_pct_mt"
        ],
        "Value": [
            f"{adata.n_obs:,}", f"{adata.n_vars:,}", 
            f"{obs['n_genes_by_counts'].median():.0f}", f"{obs['total_counts'].median():.0f}",
            f"{obs['pct_counts_mt'].median():.2f}%", f"{obs['pct_counts_ribo'].median():.2f}%",
            f"{(obs['n_genes_by_counts'] > thresholds['max_genes']).sum():,}",
            f"{(obs['n_genes_by_counts'] < thresholds['min_genes']).sum():,}",
            f"{(obs['pct_counts_mt'] > thresholds['max_pct_mt']).sum():,}"
        ]
    }
    df = pd.DataFrame(stats)

    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.axis('off')
    
    # Render table data
    table = ax.table(cellText=df.values, colLabels=df.columns, loc='center', cellLoc='left')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.4)
    
    # Format header styling
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight='bold', color='white')
            cell.set_facecolor('#4C72B0')
            
    fig.suptitle(f"Summary Statistics — {sample_id}", y=0.95)
    plt.tight_layout()
    return fig

# Violin Plots
def fig_violin(adata: sc.AnnData, sample_id: str) -> plt.Figure:
    metrics = [
        ("n_genes_by_counts", "N Genes"),
        ("total_counts",      "Total Counts (UMIs)"),
        ("pct_counts_mt",     "% MT Counts"),
        ("pct_counts_ribo",   "% Ribo Counts")
    ]
    
    fig, axes = plt.subplots(1, 4, figsize=(15, 5))
    colors = sns.color_palette("Set2", 4)
    
    for i, (col, label) in enumerate(metrics):
        sns.violinplot(y=adata.obs[col], ax=axes[i], color=colors[i], inner="box")
        axes[i].set_title(label)
        axes[i].set_ylabel("")
        
    fig.suptitle(f"QC Metrics — {sample_id}", y=1.02)
    plt.tight_layout()
    return fig

# Scatter Plot
def fig_scatter(adata: sc.AnnData, sample_id: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7, 5))
    
    scat = ax.scatter(
        x=adata.obs["total_counts"],
        y=adata.obs["n_genes_by_counts"],
        c=adata.obs["pct_counts_mt"],
        cmap="RdYlBu_r",
        alpha=0.6,
        edgecolors="none",
        s=10
    )
    
    ax.set_xlabel("Total Counts (UMIs)")
    ax.set_ylabel("N Genes")
    ax.set_title(f"Counts vs Genes (coloured by % MT) — {sample_id}")
    
    cbar = fig.colorbar(scat, ax=ax)
    cbar.set_label("% MT")
    
    plt.tight_layout()
    return fig

# Histogram: N Genes
def fig_histogram_genes(adata: sc.AnnData, sample_id: str, min_genes: int, max_genes: int) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(6, 4))
    
    sns.histplot(data=adata.obs, x="n_genes_by_counts", bins=50, ax=ax, color=sns.color_palette("Set2")[0], kde=False)
    ax.axvline(min_genes, color="red", linestyle="--", label=f"Min: {min_genes}")
    ax.axvline(max_genes, color="red", linestyle="--", label=f"Max: {max_genes}")
    
    ax.set_xlabel("N Genes")
    ax.set_ylabel("Cell Count")
    ax.set_title(f"N Genes Distribution — {sample_id}")
    ax.legend()
    
    plt.tight_layout()
    return fig

# Histogram: % MT
def fig_histogram_mt(adata: sc.AnnData, sample_id: str, max_pct_mt: float) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(6, 4))
    
    sns.histplot(data=adata.obs, x="pct_counts_mt", bins=50, ax=ax, color=sns.color_palette("Set2")[2], kde=False)
    ax.axvline(max_pct_mt, color="red", linestyle="--", label=f"Threshold: {max_pct_mt}%")
    
    ax.set_xlabel("% MT Counts")
    ax.set_ylabel("Cell Count")
    ax.set_title(f"% MT Counts Distribution — {sample_id}")
    ax.legend()
    
    plt.tight_layout()
    return fig