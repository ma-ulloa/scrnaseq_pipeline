"""
utils/qc_metrics.py
────────────────────────────────────────────────────────────────
Shared QC computation logic used by both pre- and post-filter
report scripts.

MT_PREFIXES, RIBO_PREFIXES      — gene name prefix constants
compute_qc_metrics(adata)       — flag MT/ribo genes, run sc.pp.calculate_qc_metrics
mad_bounds(series, n_mads)      — raw (unscaled) MAD lower/upper bounds
thresholds_from_mad(adata, n)   — derive the four filter thresholds for a given MAD
add_outlier_flags(adata, mads)  — write per-MAD outlier columns to adata.obs
save_qc_metrics(adata, sid)     — extract per-cell QC metrics as a DataFrame
────────────────────────────────────────────────────────────────
"""

import numpy as np
import pandas as pd
import scanpy as sc

# ── Gene-name prefixes ────────────────────────────────────────────────────────
# Human: MT-  |  Mouse: mt-  |  other conventions: MT_
HUMAN_MT      = ("MT-", "MT_")
MOUSE_MT      = ("mt-",)   
HUMAN_RIBO    = ("RPS", "RPL")
MOUSE_RIBO    = ("Rps", "Rpl")

# ── Core computation ──────────────────────────────────────────────────────────

def compute_qc_metrics(adata: sc.AnnData, species: str) -> sc.AnnData:
    """ Get metrics"""

    if species == "human":
        adata.var["mt"]   = adata.var_names.str.startswith(HUMAN_MT)
        adata.var["ribo"] = adata.var_names.str.startswith(HUMAN_RIBO)
    else:
        adata.var["mt"]   = adata.var_names.str.startswith(MOUSE_MT)
        adata.var["ribo"] = adata.var_names.str.startswith(MOUSE_RIBO)

    n_mt   = adata.var["mt"].sum()
    n_ribo = adata.var["ribo"].sum()
    print(f"[INFO] MT genes: {n_mt} | Ribo genes: {n_ribo}")
    if n_mt == 0:
        print("[WARN] No mitochondrial genes detected — check gene naming convention.")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"],
        percent_top=None, log1p=True, inplace=True,
    )
    print(f"[INFO] QC metrics computed. Cells: {adata.n_obs} | Genes: {adata.n_vars}")
    return adata


#  MAD helpers 
def mad_bounds(series: pd.Series, n_mads: float) -> tuple[float, float]:
    """
    MAD calculation
    lower = max(0, median − n_mads × MAD)
    upper =        median + n_mads × MAD
    """
    median = series.median()
    mad    = np.median(np.abs(series - median))
    return max(0.0, median - n_mads * mad), median + n_mads * mad


def thresholds_from_mad(adata: sc.AnnData, n_mads: float) -> dict:
    """
    Derive the four filter thresholds for a given MAD multiplier.
    Genes and counts: thresholds computed on log1p-transformed values,
    then converted back to linear space via expm1.
    MT %: threshold computed directly in linear space.
    """
    g_lo, g_hi = mad_bounds(adata.obs["log1p_n_genes_by_counts"], n_mads)
    c_lo, _    = mad_bounds(adata.obs["log1p_total_counts"],       n_mads)
    _, mt_hi   = mad_bounds(adata.obs["pct_counts_mt"],            n_mads)

    return {
        "min_genes":  int(np.expm1(g_lo)),
        "max_genes":  int(np.expm1(g_hi)),
        "min_counts": int(np.expm1(c_lo)),
        "max_pct_mt": float(mt_hi),
    }


# ── Outlier flags ─────────────────────────────────────────────────────────────

def add_outlier_flags(adata: sc.AnnData,
                      mad_values: list[float] | None = None) -> sc.AnnData:
    """
    Write per-MAD outlier columns to adata.obs.

    For each value m in mad_values, adds four columns:
        outlier_genes_{m}mad   — True if n_genes outside [min_genes, max_genes]
        outlier_counts_{m}mad  — True if total_counts < min_counts
        outlier_mt_{m}mad      — True if pct_counts_mt > max_pct_mt
        outlier_{m}mad         — True if flagged by ANY of the three above

    The downstream filtering script (03_filtering.py) reads the combined
    column directly, with no need to recompute thresholds.
    """
    if mad_values is None:
        mad_values = [3, 5]

    for m in mad_values:
        t   = thresholds_from_mad(adata, m)
        obs = adata.obs
        tag = f"{int(m)}mad"

        adata.obs[f"outlier_genes_{tag}"] = (
            (obs["n_genes_by_counts"] < t["min_genes"]) |
            (obs["n_genes_by_counts"] > t["max_genes"])
        )
        adata.obs[f"outlier_counts_{tag}"] = obs["total_counts"] < t["min_counts"]
        adata.obs[f"outlier_mt_{tag}"]     = obs["pct_counts_mt"] > t["max_pct_mt"]
        adata.obs[f"outlier_{tag}"]        = (
            adata.obs[f"outlier_genes_{tag}"]  |
            adata.obs[f"outlier_counts_{tag}"] |
            adata.obs[f"outlier_mt_{tag}"]
        )

    return adata


#  CSV export 

def save_qc_metrics(adata: sc.AnnData, sample_id: str) -> pd.DataFrame:
    """
    Extract per-cell QC metrics from adata.obs into a standalone DataFrame.
    Column order: sample_id, cell_id, <qc metrics>, <outlier flags>
    """
    base_cols = [
        "n_genes_by_counts",
        "log1p_n_genes_by_counts",
        "total_counts",
        "log1p_total_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
    ]
    outlier_cols = sorted(
        [c for c in adata.obs.columns if c.startswith("outlier_")]
    )

    all_cols = [c for c in base_cols + outlier_cols if c in adata.obs.columns]
    df = adata.obs[all_cols].copy()

    df.insert(0, "cell_id",   adata.obs.index)
    df.insert(0, "sample_id", sample_id)

    return df