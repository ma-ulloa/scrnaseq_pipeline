"""
utils/mad_plots.py
────────────────────────────────────────────────────────────────
MAD figures — used only in the PRE-filter report.
These plots help decide which MAD multiplier to use for filtering.
────────────────────────────────────────────────────────────────
"""

import warnings
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
import scanpy as sc

try:
    from .qc_metrics import mad_bounds, thresholds_from_mad   # when imported as package
except ImportError:
    from qc_metrics import mad_bounds, thresholds_from_mad    # when run directly

warnings.filterwarnings("ignore", category=FutureWarning)


# ── Retention table ───────────────────────────────────────────────────────────

def compute_retention_table(adata: sc.AnnData,
                             mad_values: list[float] | None = None) -> pd.DataFrame:
    """
    Returns a DataFrame with one row per MAD value and columns:
        mad, min_genes, max_genes, min_counts, max_pct_mt,
        pass_genes, pass_counts, pass_mt, pass_combined,
        pct_genes,  pct_counts,  pct_mt,  pct_combined
    """
    if mad_values is None:
        mad_values = [1, 2, 3, 4, 5]

    rows = []
    for m in mad_values:
        t   = thresholds_from_mad(adata, m)
        obs = adata.obs

        pass_genes    = (obs["n_genes_by_counts"] >= t["min_genes"]) & \
                        (obs["n_genes_by_counts"] <= t["max_genes"])
        pass_counts   = obs["total_counts"]  >= t["min_counts"]
        pass_mt       = obs["pct_counts_mt"] <= t["max_pct_mt"]
        pass_combined = pass_genes & pass_counts & pass_mt

        rows.append({
            "mad":           m,
            "min_genes":     t["min_genes"],
            "max_genes":     t["max_genes"],
            "min_counts":    t["min_counts"],
            "max_pct_mt":    round(t["max_pct_mt"], 2),
            "pass_genes":    int(pass_genes.sum()),
            "pass_counts":   int(pass_counts.sum()),
            "pass_mt":       int(pass_mt.sum()),
            "pass_combined": int(pass_combined.sum()),
            "pct_genes":     round(pass_genes.mean()    * 100, 1),
            "pct_counts":    round(pass_counts.mean()   * 100, 1),
            "pct_mt":        round(pass_mt.mean()       * 100, 1),
            "pct_combined":  round(pass_combined.mean() * 100, 1),
        })

    return pd.DataFrame(rows)


#  Figure 1: Retention curve 

def fig_mad_retention_curve(adata: sc.AnnData,
                             sample_id: str,
                             mad_values: list[float] | None = None,
                             highlight_mad: Optional[float] = 3.0) -> plt.Figure:
    """
    Line plot: MAD multiplier (x-axis) vs % cells retained (y-axis).

    One line per filter metric (genes, counts, MT %) plus the combined filter.
    A vertical dotted line marks the chosen MAD and annotates the cell count
    that would be retained at that setting.
    """
    if mad_values is None:
        mad_values = [1, 2, 3, 4, 5]

    df = compute_retention_table(adata, mad_values)

    palette = {
        "Genes (low + high)": "#4C72B0",
        "Counts (low)":       "#55A868",
        "MT % (high)":        "#C44E52",
        "Combined":           "#1a1a1a",
    }
    styles  = {"Genes (low + high)": "-",  "Counts (low)": "-",
               "MT % (high)": "-",         "Combined": "--"}
    markers = {"Genes (low + high)": "o",  "Counts (low)": "s",
               "MT % (high)": "^",         "Combined": "D"}
    cols    = {"Genes (low + high)": "pct_genes",  "Counts (low)": "pct_counts",
               "MT % (high)": "pct_mt",            "Combined": "pct_combined"}

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.suptitle(f"MAD Threshold Sensitivity — {sample_id}",
                 fontsize=13, fontweight="bold", y=1.01)

    for label, col in cols.items():
        ax.plot(df["mad"], df[col],
                label=label,
                color=palette[label],
                linestyle=styles[label],
                marker=markers[label],
                linewidth=2.5 if label == "Combined" else 2,
                markersize=7)
        ax.annotate(f"{df[col].iloc[-1]:.1f}%",
                    xy=(df["mad"].iloc[-1], df[col].iloc[-1]),
                    xytext=(6, 0), textcoords="offset points",
                    fontsize=8, color=palette[label], va="center")

    if highlight_mad is not None and highlight_mad in df["mad"].values:
        row = df[df["mad"] == highlight_mad].iloc[0]
        ax.axvline(highlight_mad, color="grey", linestyle=":", linewidth=1.5,
                   label=f"Selected MAD = {int(highlight_mad)}×")
        ax.annotate(
            f"  {row['pass_combined']:,} cells\n  ({row['pct_combined']:.1f}%)",
            xy=(highlight_mad, row["pct_combined"]),
            xytext=(8, -12), textcoords="offset points",
            fontsize=8.5, color="black",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="grey", alpha=0.8),
        )

    ax.set_xlabel("MAD multiplier", fontsize=11)
    ax.set_ylabel("Cells retained (%)", fontsize=11)
    ax.set_ylim(0, 105)
    ax.set_xticks(df["mad"].tolist())
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.legend(frameon=True, fontsize=9)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    sns.despine(ax=ax)
    plt.tight_layout()
    return fig


# ── Figure 2: Sensitivity table ───────────────────────────────────────────────

def fig_mad_sensitivity_table(adata: sc.AnnData,
                               sample_id: str,
                               mad_values: list[float] | None = None,
                               highlight_mad: Optional[float] = 3.0) -> plt.Figure:
    """
    Formatted table: thresholds + per-metric and combined retention for each MAD.
    The row matching highlight_mad is shaded blue to indicate the chosen setting.
    """
    if mad_values is None:
        mad_values = [1, 2, 3, 4, 5]

    df = compute_retention_table(adata, mad_values)

    col_labels = [
        "MAD",
        "min\ngenes", "max\ngenes", "min\ncounts", "max\nMT %",
        "pass genes\nn  (%)",
        "pass counts\nn  (%)",
        "pass MT\nn  (%)",
        "COMBINED\nn  (%)",
    ]
    table_data = []
    for _, r in df.iterrows():
        table_data.append([
            f"{int(r['mad'])}×",
            f"{r['min_genes']:,}",
            f"{r['max_genes']:,}",
            f"{r['min_counts']:,}",
            f"{r['max_pct_mt']:.1f}",
            f"{r['pass_genes']:,}  ({r['pct_genes']:.1f}%)",
            f"{r['pass_counts']:,}  ({r['pct_counts']:.1f}%)",
            f"{r['pass_mt']:,}  ({r['pct_mt']:.1f}%)",
            f"{r['pass_combined']:,}  ({r['pct_combined']:.1f}%)",
        ])

    fig_h = 0.55 + 0.42 * len(mad_values)
    fig, ax = plt.subplots(figsize=(14, fig_h))
    ax.axis("off")
    fig.suptitle(
        f"MAD Sensitivity Summary — {sample_id}  "
        f"(total cells before filter: {adata.n_obs:,})",
        fontsize=12, fontweight="bold", y=1.02,
    )

    tbl = ax.table(cellText=table_data, colLabels=col_labels,
                   cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.6)

    HEADER_BG    = "#2c3e50"
    HIGHLIGHT_BG = "#d6eaf8"
    ALT_BG       = "#f8f9fa"
    COMBINED_COL = len(col_labels) - 1   # last column

    for col_idx in range(len(col_labels)):
        cell = tbl[0, col_idx]
        cell.set_facecolor(HEADER_BG)
        cell.set_text_props(color="white", fontweight="bold")

    for row_idx in range(1, len(table_data) + 1):
        mad_val      = float(df.iloc[row_idx - 1]["mad"])
        is_highlight = highlight_mad is not None and mad_val == highlight_mad
        for col_idx in range(len(col_labels)):
            cell = tbl[row_idx, col_idx]
            if is_highlight:
                cell.set_facecolor(HIGHLIGHT_BG)
            elif row_idx % 2 == 0:
                cell.set_facecolor(ALT_BG)
            else:
                cell.set_facecolor("white")
            if col_idx == COMBINED_COL:
                cell.set_text_props(fontweight="bold")

    plt.tight_layout()
    return fig