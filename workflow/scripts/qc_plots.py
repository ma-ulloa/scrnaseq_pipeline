"""
qc_plots.py
────────────────────────────────────────────────────────────────
All QC plot functions use in 01_qc_stats_report.

Each function takes an AnnData object and returns a Plotly Figure.

Usage:
    from qc_plots import (
        fig_summary_table,
        fig_violin,
        fig_scatter,
        fig_histogram_genes,
        fig_histogram_mt,
    )
────────────────────────────────────────────────────────────────
"""

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import scanpy as sc


# Shared colour palette
COLORS = px.colors.qualitative.Set2 # TO DO: Include a color function


# Summary Table 
def fig_summary_table(adata: sc.AnnData, sample_id: str, thresholds: dict) -> go.Figure:
    """
    Summary statistics table showing cell counts, medians,
    and how many cells fall outside each threshold.

    Parameters
    ----------
    adata       : AnnData with QC metrics already computed
    sample_id   : label shown in the figure title
    thresholds  : dict with keys min_genes, max_genes, max_pct_mt
    """
    obs = adata.obs
    stats = {
        "Metric": [
            "Total cells",
            "Total genes",
            "Median genes / cell",
            "Median counts / cell",
            "Median % MT",
            "Median % Ribo",
            "Cells > max_genes threshold",
            "Cells < min_genes threshold",
            "Cells > max_pct_mt threshold",
        ],
        "Value": [
            f"{adata.n_obs:,}",
            f"{adata.n_vars:,}",
            f"{obs['n_genes_by_counts'].median():.0f}",
            f"{obs['total_counts'].median():.0f}",
            f"{obs['pct_counts_mt'].median():.2f}%",
            f"{obs['pct_counts_ribo'].median():.2f}%",
            f"{(obs['n_genes_by_counts'] > thresholds['max_genes']).sum():,}",
            f"{(obs['n_genes_by_counts'] < thresholds['min_genes']).sum():,}",
            f"{(obs['pct_counts_mt'] > thresholds['max_pct_mt']).sum():,}",
        ],
    }

    fig = go.Figure(data=[go.Table(
        header=dict(
            values=["<b>Metric</b>", "<b>Value</b>"],
            fill_color="#4C72B0",
            font=dict(color="white", size=13),
            align="left",
        ),
        cells=dict(
            values=[stats["Metric"], stats["Value"]],
            fill_color=[["#f9f9f9", "white"] * 10],
            align="left",
            font=dict(size=12),
        ),
    )])
    fig.update_layout(
        title=f"Summary Statistics — {sample_id}",
        height=370,
        margin=dict(t=50, b=10),
    )
    return fig


# Violin Plots 

def fig_violin(adata: sc.AnnData, sample_id: str) -> go.Figure:
    """
    Four-panel violin plot:
        n_genes_by_counts | total_counts | pct_counts_mt | pct_counts_ribo

    Parameters
    ----------
    adata      : AnnData with QC metrics
    sample_id  : label shown in the figure title
    """
    metrics = [
        ("n_genes_by_counts", "N Genes"),
        ("total_counts",      "Total Counts (UMIs)"),
        ("pct_counts_mt",     "% MT Counts"),
        ("pct_counts_ribo",   "% Ribo Counts"),
    ]
    fig = make_subplots(rows=1, cols=4, subplot_titles=[m[1] for m in metrics])

    for i, (col, label) in enumerate(metrics, start=1):
        fig.add_trace(
            go.Violin(
                y=adata.obs[col].values,
                name=label,
                box_visible=True,
                meanline_visible=True,
                fillcolor=COLORS[i % len(COLORS)],
                opacity=0.7,
                line_color="black",
                showlegend=False,
            ),
            row=1, col=i,
        )

    fig.update_layout(
        title_text=f"QC Metrics — {sample_id}",
        height=450,
        template="plotly_white",
    )
    return fig


# Scatter Plot 

def fig_scatter(adata: sc.AnnData, sample_id: str) -> go.Figure:
    """
    Scatter plot: total_counts (x) vs n_genes_by_counts (y),
    coloured by % MT counts.

    Useful for spotting doublets (high counts + high genes)
    and dying cells (high % MT).

    Parameters
    ----------
    adata      : AnnData with QC metrics
    sample_id  : label shown in the figure title
    """
    df = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()

    fig = px.scatter(
        df,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
        color_continuous_scale="RdYlBu_r",
        labels={
            "total_counts":      "Total Counts (UMIs)",
            "n_genes_by_counts": "N Genes",
            "pct_counts_mt":     "% MT",
        },
        title=f"Counts vs Genes (coloured by % MT) — {sample_id}",
        opacity=0.5,
        template="plotly_white",
    )
    fig.update_layout(height=450)
    return fig


#  Histogram: N Genes 

def fig_histogram_genes(
    adata: sc.AnnData,
    sample_id: str,
    min_genes: int,
    max_genes: int,
) -> go.Figure:
    """
    Histogram of n_genes_by_counts with vertical threshold lines.

    Parameters
    ----------
    adata      : AnnData with QC metrics
    sample_id  : label shown in the figure title
    min_genes  : lower threshold (red dashed line)
    max_genes  : upper threshold (red dashed line)
    """
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=adata.obs["n_genes_by_counts"],
        nbinsx=50,
        marker_color=COLORS[0],
        opacity=0.8,
        name="N Genes",
    ))
    fig.add_vline(
        x=min_genes, line_dash="dash", line_color="red",
        annotation_text=f"Min: {min_genes}",
        annotation_position="top left",
    )
    fig.add_vline(
        x=max_genes, line_dash="dash", line_color="red",
        annotation_text=f"Max: {max_genes}",
        annotation_position="top right",
    )
    fig.update_layout(
        title=f"N Genes Distribution — {sample_id}",
        xaxis_title="N Genes",
        yaxis_title="Cell Count",
        template="plotly_white",
        height=400,
    )
    return fig


# Histogram: % MT 

def fig_histogram_mt(
    adata: sc.AnnData,
    sample_id: str,
    max_pct_mt: float,
) -> go.Figure:
    """
    Histogram of pct_counts_mt with a vertical threshold line.

    Parameters
    ----------
    adata      : AnnData with QC metrics
    sample_id  : label shown in the figure title
    max_pct_mt : upper threshold (red dashed line)
    """
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=adata.obs["pct_counts_mt"],
        nbinsx=50,
        marker_color=COLORS[2],
        opacity=0.8,
        name="% MT",
    ))
    fig.add_vline(
        x=max_pct_mt, line_dash="dash", line_color="red",
        annotation_text=f"Threshold: {max_pct_mt}%",
        annotation_position="top right",
    )
    fig.update_layout(
        title=f"% MT Counts Distribution — {sample_id}",
        xaxis_title="% MT Counts",
        yaxis_title="Cell Count",
        template="plotly_white",
        height=400,
    )
    return fig