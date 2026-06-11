"""
utils/qc_plots.py
────────────────────────────────────────────────────────────────
All QC plot functions — per-sample and multi-sample (cohort).

Per-sample (called by 01_qc_stats_report.py):
    fig_summary_table(adata, sample_id, thresholds)
    fig_violin(adata, sample_id)
    fig_scatter(adata, sample_id)
    fig_histogram_genes(adata, sample_id, min_genes, max_genes)
    fig_histogram_mt(adata, sample_id, max_pct_mt)
    fig_multi_mad_histograms(adata, sample_id)
    fig_mt_decay_curve(adata, sample_id)

Multi-sample / cohort (called by 02_multisample_report.py):
    fig_cohort_bar_plots(merged_df)   → list[Figure]
    fig_cohort_scatters(merged_df)    → list[Figure]
    fig_cohort_violins(merged_df)     → list[Figure]
    fig_metadata_breakdown(merged_df) → Figure

NOTE on merged_df:
    One row per CELL (not per sample). Columns include all qc metrics
    from save_qc_metrics() plus 'stage' and, if metadata was joined,
    'patient', 'timepoint', 'cells', 'response'.
    Barplot errorbar="se" therefore shows SE across cells within a sample,
    NOT across biological replicates — that is intentional for QC purposes.
────────────────────────────────────────────────────────────────
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import seaborn as sns
import scanpy as sc

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Global style ──────────────────────────────────────────────────────────────
sns.set_theme(style="whitegrid", font_scale=1.0)
plt.rcParams.update({
    "figure.titlesize": 14,
    "axes.titlesize":   12,
    "axes.labelsize":   10,
    "xtick.labelsize":   9,
    "ytick.labelsize":   9,
})

# Metadata colour maps reused across plots
RESPONSE_COLORS = {"responder": "#4C72B0", "control": "#DD8452"}
CELLS_COLORS    = {"epi_stroma": "#9FCBAD", "immune": "#4A4466"}


# ══════════════════════════════════════════════════════════════════════════════
#  PER-SAMPLE PLOTS
# ══════════════════════════════════════════════════════════════════════════════

def fig_summary_table(adata: sc.AnnData, sample_id: str, thresholds: dict) -> plt.Figure:
    """Clean summary statistics table for a single sample."""
    obs = adata.obs
    stats = {
        "Metric": [
            "Total cells", "Total genes",
            "Median genes / cell", "Median counts / cell",
            "Median % MT", "Median % Ribo",
            "Cells > max_genes", "Cells < min_genes", "Cells > max_pct_mt",
        ],
        "Value": [
            f"{adata.n_obs:,}", f"{adata.n_vars:,}",
            f"{obs['n_genes_by_counts'].median():.0f}",
            f"{obs['total_counts'].median():.0f}",
            f"{obs['pct_counts_mt'].median():.2f}%",
            f"{obs['pct_counts_ribo'].median():.2f}%",
            f"{(obs['n_genes_by_counts'] > thresholds['max_genes']).sum():,}",
            f"{(obs['n_genes_by_counts'] < thresholds['min_genes']).sum():,}",
            f"{(obs['pct_counts_mt'] > thresholds['max_pct_mt']).sum():,}",
        ],
    }
    df = pd.DataFrame(stats)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.axis("off")
    tbl = ax.table(cellText=df.values, colLabels=df.columns, loc="center", cellLoc="left")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.2, 1.4)
    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_text_props(weight="bold", color="white")
            cell.set_facecolor("#4C72B0")
        elif row % 2 == 0:
            cell.set_facecolor("#f7f7f7")
    fig.suptitle(f"Summary Statistics — {sample_id}", y=0.97)
    plt.tight_layout()
    return fig


def fig_violin(adata: sc.AnnData, sample_id: str) -> plt.Figure:
    """
    4-panel violin: linear UMIs, log10 UMIs, N genes, % MT.
    inner="quart" shows median + IQR without clutter.
    """
    palette = ["#69b3a2", "#404080", "#8fbc8f", "#e06666"]

    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    fig.suptitle(f"QC Metrics Distribution — {sample_id}", fontsize=14, fontweight="bold")

    # 1. Total UMI Counts (linear)
    sns.violinplot(y=adata.obs["total_counts"], ax=axes[0], color=palette[0], inner="quart")
    axes[0].set_title("Total UMI Counts (Linear)")
    axes[0].set_ylabel("UMIs / Cell")

    # 2. Total UMI Counts (log10)
    log10_counts = np.log10(adata.obs["total_counts"] + 1)
    sns.violinplot(y=log10_counts, ax=axes[1], color=palette[1], inner="quart")
    axes[1].set_title("Total UMI Counts (Log₁₀)")
    axes[1].set_ylabel(r"$\log_{10}(\mathrm{UMIs} + 1)$")

    # 3. Unique genes
    sns.violinplot(y=adata.obs["n_genes_by_counts"], ax=axes[2], color=palette[2], inner="quart")
    axes[2].set_title("Unique Genes Detected")
    axes[2].set_ylabel("Genes / Cell")

    # 4. MT percentage
    sns.violinplot(y=adata.obs["pct_counts_mt"], ax=axes[3], color=palette[3], inner="quart")
    axes[3].set_title("Mitochondrial Percentage")
    axes[3].set_ylabel("% MT Counts / Cell")

    plt.tight_layout()
    return fig


def fig_scatter(adata: sc.AnnData, sample_id: str) -> plt.Figure:
    """
    Counts vs genes scatter coloured by % MT.
    Colorbar capped at 99th percentile to avoid outlier saturation.
    """
    obs = adata.obs
    fig, ax = plt.subplots(figsize=(7, 5))
    scat = ax.scatter(
        x=obs["total_counts"],
        y=obs["n_genes_by_counts"],
        c=obs["pct_counts_mt"],
        cmap="RdYlBu_r",
        alpha=0.5,
        edgecolors="none",
        s=8,
        vmin=0,
        vmax=obs["pct_counts_mt"].quantile(0.99),
    )
    ax.set_xlabel("Total Counts (UMIs)")
    ax.set_ylabel("N Genes")
    ax.set_title(f"Counts vs Genes (coloured by % MT) — {sample_id}")
    cbar = fig.colorbar(scat, ax=ax, pad=0.01)
    cbar.set_label("% MT")
    plt.tight_layout()
    return fig


def fig_histogram_genes(adata: sc.AnnData, sample_id: str,
                         min_genes: int, max_genes: int) -> plt.Figure:
    """N Genes histogram with MAD-derived threshold lines."""
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.histplot(data=adata.obs, x="n_genes_by_counts", bins=50, ax=ax,
                 color=sns.color_palette("Set2")[0], kde=False)
    ax.axvline(min_genes, color="crimson",    linestyle="--", label=f"Min: {min_genes}")
    ax.axvline(max_genes, color="darkorange", linestyle="--", label=f"Max: {max_genes}")
    ax.set_xlabel("N Genes")
    ax.set_ylabel("Cell Count")
    ax.set_title(f"N Genes Distribution — {sample_id}")
    ax.legend()
    plt.tight_layout()
    return fig


def fig_histogram_mt(adata: sc.AnnData, sample_id: str, max_pct_mt: float) -> plt.Figure:
    """% MT histogram with KDE and threshold line."""
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.histplot(data=adata.obs, x="pct_counts_mt", bins=50, ax=ax,
                 color=sns.color_palette("Set2")[2], kde=True,
                 kde_kws={"bw_adjust": 1.5})
    ax.axvline(max_pct_mt, color="crimson", linestyle="--",
               linewidth=1.5, label=f"Threshold: {max_pct_mt}%")
    ax.set_xlabel("% MT Counts")
    ax.set_ylabel("Cell Count")
    ax.set_title(f"% MT Counts Distribution — {sample_id}")
    ax.legend()
    plt.tight_layout()
    return fig


def fig_mt_decay_curve(adata: sc.AnnData, sample_id: str) -> plt.Figure:
    """
    Data-retention curve: % cells and % total UMIs retained at each MT cutoff.
    X-axis shows cutoffs 1–25%; annotations are staggered to avoid overlap.
    """
    cutoffs = [1, 5, 10, 15, 20, 25]
    total_cells = adata.n_obs
    total_umis  = adata.obs["total_counts"].sum()

    rows = []
    for limit in cutoffs:
        sub      = adata.obs[adata.obs["pct_counts_mt"] <= limit]
        cell_pct = len(sub) / total_cells * 100
        umi_pct  = sub["total_counts"].sum() / total_umis * 100
        rows.append({"Cutoff": limit, "Cells (%)": cell_pct, "UMIs (%)": umi_pct})
    df = pd.DataFrame(rows)

    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax2 = ax1.twinx()

    p1, = ax1.plot(df["Cutoff"], df["Cells (%)"], marker="o",
                   color="#1f77b4", linewidth=2.5, label="Cells Retained (%)")
    p2, = ax2.plot(df["Cutoff"], df["UMIs (%)"],  marker="s",
                   color="#ff7f0e", linewidth=2.5, linestyle=":", label="UMIs Retained (%)")

    ax1.set_xlabel("Maximum Allowed MT Counts (%)")
    ax1.set_ylabel("Cells Retained (%)",      color="#1f77b4")
    ax2.set_ylabel("Total UMIs Retained (%)", color="#ff7f0e")
    ax1.set_xticks(cutoffs)
    ax1.grid(True, linestyle=":", alpha=0.6)

    # Stagger annotations to prevent overlap
    for i, row in df.iterrows():
        y_offset_cell = 6 if i % 2 == 0 else -14
        y_offset_umi  = -14 if i % 2 == 0 else 6
        ax1.annotate(f"{row['Cells (%)']:.1f}%",
                     (row["Cutoff"], row["Cells (%)"]),
                     textcoords="offset points", xytext=(0, y_offset_cell),
                     ha="center", color="#1f77b4", fontsize=8, fontweight="bold")
        ax2.annotate(f"{row['UMIs (%)']:.1f}%",
                     (row["Cutoff"], row["UMIs (%)"]),
                     textcoords="offset points", xytext=(0, y_offset_umi),
                     ha="center", color="#ff7f0e", fontsize=8)

    lines = [p1, p2]
    ax1.legend(lines, [l.get_label() for l in lines], loc="lower right", frameon=True)
    ax1.set_title(f"Data Retention vs MT Cutoff — {sample_id}", fontweight="bold")
    plt.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  MULTI-SAMPLE
# ══════════════════════════════════════════════════════════════════════════════

def fig_cohort_bar_plots(merged_df: pd.DataFrame) -> list:
    """
    Bar plots comparing per-sample distributions across the cohort.
    One figure per QC metric; bars coloured by pre/post stage when both present.

    Returns: list of plt.Figure
    """
    figures = []
    color_by = "cells" if "cells" in merged_df.columns else "sample_id"
    palette  = CELLS_COLORS if color_by == "cells" else "tab10"

    bar_metrics = [
        ("total_counts",       "Mean Total Counts (UMIs) per Sample",   "Mean Total Counts"),
        ("n_genes_by_counts",  "Mean Detected Genes per Sample",         "Mean Unique Genes"),
        ("pct_counts_mt",      "Mean Mitochondrial % per Sample",        "Mean MT %"),
        ("pct_counts_ribo",    "Mean Ribosomal % per Sample",            "Mean Ribo %"),
    ]

    sample_order = sorted(merged_df["sample_id"].unique())

    for col_name, title, ylabel in bar_metrics:
        if col_name not in merged_df.columns:
            continue
        fig, ax = plt.subplots(figsize=(12, 6))
        sns.barplot(
            data=merged_df,
            x="sample_id", y=col_name,
            order=sample_order,
            hue=color_by,
            palette= palette,
            edgecolor="0.2",
            ax=ax,
        )
        ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_xlabel("Sample ID", fontsize=11)
        ax.tick_params(axis="x", rotation=45)
        if color_by == "cells":
            handles = [mpatches.Patch(color=v, label=k) for k, v in CELLS_COLORS.items()]
            ax.legend(handles=handles, title="Cell fraction", frameon=False, fontsize=9)
        plt.tight_layout()
        figures.append(fig)

    return figures


def fig_cohort_scatters(merged_df: pd.DataFrame) -> list:
    """
    Diagnostic scatter plots across the full cohort.
    Plot A: Total counts vs N genes, coloured by % MT.
    Plot B: Total counts vs % MT, coloured by sample_id.

    Returns: list of plt.Figure
    """
    figures = []
    sample_order = sorted(merged_df["sample_id"].unique())

    # Plot A: counts vs genes, colour = %MT
    if {"total_counts", "n_genes_by_counts", "pct_counts_mt"}.issubset(merged_df.columns):
        fig, ax = plt.subplots(figsize=(10, 7))
        sc_plot = ax.scatter(
            merged_df["total_counts"],
            merged_df["n_genes_by_counts"],
            c=merged_df["pct_counts_mt"],
            cmap="magma",
            alpha=0.4,
            edgecolors="none",
            s=6,
            vmax=merged_df["pct_counts_mt"].quantile(0.99),
        )
        ax.set_xlabel("Total Counts (UMIs)", fontsize=11)
        ax.set_ylabel("Number of Unique Genes", fontsize=11)
        ax.set_title("Transcripts vs Detected Genes — All Samples", fontsize=13,
                     fontweight="bold", pad=12)
        fig.colorbar(sc_plot, ax=ax, label="Mitochondrial %", pad=0.01)
        plt.tight_layout()
        figures.append(fig)

    # Plot B: counts vs %MT, colour = sample_id
    if {"total_counts", "pct_counts_mt"}.issubset(merged_df.columns):
        fig, ax = plt.subplots(figsize=(10, 6))
        palette = sns.color_palette("tab10", n_colors=len(sample_order))
        sns.scatterplot(
            data=merged_df,
            x="total_counts", y="pct_counts_mt",
            hue="sample_id", hue_order=sample_order,
            palette=palette,
            alpha=0.35,
            linewidth=0,
            s=6,
            ax=ax,
        )
        ax.set_xlabel("Total Counts (UMIs)", fontsize=11)
        ax.set_ylabel("Mitochondrial %",     fontsize=11)
        ax.set_title("Transcripts vs MT Percentage — All Samples", fontsize=13,
                     fontweight="bold", pad=12)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left",
                  title="Sample", fontsize=8, frameon=False)
        plt.tight_layout()
        figures.append(fig)

    return figures


def fig_cohort_violins(merged_df: pd.DataFrame) -> list:
    """
    Per-metric violin plots with all samples on the x-axis.
    Coloured by 'cells' metadata column (epi_stroma / immune) when available,
    otherwise by sample_id.

    Returns: list of plt.Figure (one per metric)
    """
    metrics = [
        ("n_genes_by_counts", "N Genes / Cell"),
        ("total_counts",      "Total Counts (UMIs)"),
        ("pct_counts_mt",     "% MT Counts"),
        ("pct_counts_ribo",   "% Ribo Counts"),
    ]
    color_by = "cells" if "cells" in merged_df.columns else "sample_id"
    palette  = CELLS_COLORS if color_by == "cells" else "tab10"
    sample_order = sorted(merged_df["sample_id"].unique())
    figures = []

    for col, label in metrics:
        if col not in merged_df.columns:
            continue
        fig, ax = plt.subplots(figsize=(max(10, len(sample_order) * 0.9), 5))
        sns.violinplot(
            data=merged_df,
            x="sample_id", y=col,
            order=sample_order,
            hue=color_by,
            palette=palette,
            inner="quartile",
            linewidth=0.7,
            cut=0,
            dodge=False,
            legend=False,
            ax=ax,
        )
        ax.set_title(f"{label} — All Samples", fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(label)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)

        if color_by == "cells":
            handles = [mpatches.Patch(color=v, label=k) for k, v in CELLS_COLORS.items()]
            ax.legend(handles=handles, title="Cell fraction", frameon=False, fontsize=9)

        plt.tight_layout()
        figures.append(fig)

    return figures


def fig_metadata_breakdown(merged_df: pd.DataFrame) -> plt.Figure | None:
    """
    Two-panel bar chart using metadata columns from merged_df.
    Left:  N cells per sample, stacked by cell fraction (epi_stroma / immune).
    Right: N cells per patient × timepoint, stacked by cell fraction.

    Returns None if 'cells' column is absent.
    """
    if "cells" not in merged_df.columns:
        return None

    # One row per cell in merged_df; count occurrences
    sample_order = sorted(merged_df["sample_id"].unique())
    cell_counts  = (merged_df.groupby(["sample_id", "cells"])
                             .size().reset_index(name="n_cells"))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)

    # Left — per sample
    pivot = (cell_counts.pivot(index="sample_id", columns="cells", values="n_cells")
                        .fillna(0)
                        .reindex(sample_order))
    pivot.plot(kind="bar", stacked=True, ax=axes[0],
               color=[CELLS_COLORS.get(c, "#aaa") for c in pivot.columns],
               edgecolor="none", width=0.7)
    axes[0].set_title("N Cells per Sample (by cell fraction)", fontweight="bold")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("N Cells")
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45, ha="right", fontsize=8)
    axes[0].legend(title="Cell fraction", frameon=False)

    # Right — per patient × timepoint
    if {"timepoint", "patient"}.issubset(merged_df.columns):
        tp = (merged_df.groupby(["patient", "timepoint", "cells"])
                       .size().reset_index(name="n_cells"))
        tp["pt_tp"] = tp["patient"] + "\n" + tp["timepoint"]
        tp_piv = (tp.pivot_table(index="pt_tp", columns="cells",
                                  values="n_cells", aggfunc="sum")
                    .fillna(0))
        tp_piv.plot(kind="bar", stacked=True, ax=axes[1],
                    color=[CELLS_COLORS.get(c, "#aaa") for c in tp_piv.columns],
                    edgecolor="none", width=0.7)
        axes[1].set_title("N Cells per Patient × Timepoint", fontweight="bold")
        axes[1].set_xlabel("")
        axes[1].set_ylabel("N Cells")
        axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45, ha="right", fontsize=8)
        axes[1].legend(title="Cell fraction", frameon=False)
    else:
        axes[1].set_visible(False)

    fig.suptitle("Sample Metadata Overview", fontsize=13, fontweight="bold")
    return fig
