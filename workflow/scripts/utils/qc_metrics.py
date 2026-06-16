"""
utils/qc_metrics.py
────────────────────────────────────────────────────────────────
Shared QC computation logic used by both pre- and post-filter
report scripts.

MT_PREFIXES, RIBO_PREFIXES          — gene name prefix constants
compute_qc_metrics(adata, species)  — flag MT/ribo genes, run sc.pp.calculate_qc_metrics
run_scrublet(adata, rate)           — doublet detection via sc.pp.scrublet; adds doublet_score / predicted_doublet
mad_bounds(series, n_mads)          — raw (unscaled) MAD lower/upper bounds
thresholds_from_mad(adata, n)       — derive display thresholds at a given MAD (kept for sensitivity plots)
add_outlier_flags(adata, mads)      — write per-MAD outlier columns to adata.obs (kept for MAD sensitivity plots)
THRESHOLD_COLS                      — ordered list of per-sample threshold column names
get_sample_thresholds(path, sid)    — read per-sample threshold dict from samples.csv
compute_outlier_flags(adata, thr)   — compute 'outlier' flag from per-sample thresholds
resolve_display_thresholds(adata, thr) — effective threshold values for histogram annotation
save_qc_metrics(adata, sid)         — extract per-cell QC metrics as a DataFrame
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

# ── Per-sample threshold column names ────────────────────────────────────────
THRESHOLD_COLS = [
    "mad_mt",
    "mad_counts_lower",
    "mad_counts_upper",
    "mad_genes_lower",
    "mad_genes_upper",
    "threshold_mt_upper",
    "threshold_counts_upper",
    "threshold_counts_lower",
    "threshold_genes_lower",
    "threshold_genes_upper",
    "min_cells_per_gene",
]

# ── Core computation ──────────────────────────────────────────────────────────

def run_scrublet(adata: sc.AnnData, expected_doublet_rate: float = 0.05) -> sc.AnnData:
    """
    Run Scrublet doublet detection on raw counts via sc.pp.scrublet.
    Adds 'doublet_score' and 'predicted_doublet' to adata.obs.
    Must be called before any normalisation.
    """
    sc.pp.scrublet(adata, expected_doublet_rate=expected_doublet_rate)
    n_doublets = int(adata.obs["predicted_doublet"].sum())
    pct = n_doublets / adata.n_obs * 100
    print(
        f"[INFO] Scrublet: {n_doublets:,} predicted doublets "
        f"out of {adata.n_obs:,} cells ({pct:.1f}%)"
    )
    return adata


def compute_qc_metrics(adata: sc.AnnData, species: str) -> sc.AnnData:
    """ Get metrics"""

    if species == "human":
        adata.var["mt"]   = adata.var_names.str.startswith(HUMAN_MT)
        adata.var["ribo"] = adata.var_names.str.startswith(HUMAN_RIBO)
    else: # To do mouse
        adata.var["mt"]   = adata.var_names.str.startswith(MOUSE_MT)
        adata.var["ribo"] = adata.var_names.str.startswith(MOUSE_RIBO)

    n_mt   = adata.var["mt"].sum()
    n_ribo = adata.var["ribo"].sum()
    print(f"[INFO] MT genes: {n_mt} | Ribo genes: {n_ribo}")
    if n_mt == 0:
        print("[WARN] No mitochondrial genes detected — MT genes may have been removed during filtering.")

    # If MT genes were removed upstream, pct_counts_mt in .obs already holds the
    # correct pre-removal values. calculate_qc_metrics would overwrite them with 0,
    # so we save and restore them.
    saved_mt = (
        adata.obs["pct_counts_mt"].copy()
        if n_mt == 0 and "pct_counts_mt" in adata.obs.columns
        else None
    )

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"],
        percent_top=None, log1p=True, inplace=True,
    )

    if saved_mt is not None:
        adata.obs["pct_counts_mt"] = saved_mt.values

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
    Kept for MAD sensitivity plots.
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


# ── Outlier flags (MAD sensitivity — kept for diagnostic plots) ───────────────

def add_outlier_flags(adata: sc.AnnData,
                      mad_values: list[float] | None = None) -> sc.AnnData:
    """
    Write per-MAD outlier columns to adata.obs for diagnostic sensitivity plots.

    For each value m in mad_values, adds four columns:
        outlier_genes_{m}mad   — True if n_genes outside [min_genes, max_genes]
        outlier_counts_{m}mad  — True if total_counts < min_counts
        outlier_mt_{m}mad      — True if pct_counts_mt > max_pct_mt
        outlier_{m}mad         — True if flagged by ANY of the three above
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


# ── Per-sample threshold helpers ──────────────────────────────────────────────

def get_sample_thresholds(samples_path: str, sample_id: str) -> dict:
    """
    Read per-sample QC threshold columns from samples CSV.
    Returns a dict keyed by THRESHOLD_COLS; values are float or None.
    Empty/NaN cells in the CSV become None (= skip that filter).
    """
    df  = pd.read_csv(samples_path)
    row = df[df.iloc[:, 0].astype(str) == sample_id]
    if row.empty:
        raise ValueError(f"sample_id '{sample_id}' not found in {samples_path}")
    row = row.iloc[0]

    result = {}
    for col in THRESHOLD_COLS:
        if col not in row.index or pd.isna(row[col]) or str(row[col]).strip() == "":
            result[col] = None
        else:
            result[col] = float(row[col])
    return result


def compute_outlier_flags(adata: sc.AnnData, thresholds: dict) -> sc.AnnData:
    """
    Compute per-cell outlier flag from a per-sample threshold dict.
    Adds 'outlier' column to adata.obs (True = exclude).

    MAD thresholds use log-space for genes/counts (consistent with thresholds_from_mad).
    Any key set to None is skipped.
    """
    obs   = adata.obs
    flags = pd.Series(False, index=adata.obs_names)

    def _lo(series, m): return np.expm1(mad_bounds(series, m)[0])
    def _hi(series, m): return np.expm1(mad_bounds(series, m)[1])

    # MAD-based filters
    if thresholds.get("mad_genes_lower") is not None:
        flags |= obs["n_genes_by_counts"] < _lo(obs["log1p_n_genes_by_counts"],
                                                 thresholds["mad_genes_lower"])
    if thresholds.get("mad_genes_upper") is not None:
        flags |= obs["n_genes_by_counts"] > _hi(obs["log1p_n_genes_by_counts"],
                                                 thresholds["mad_genes_upper"])
    if thresholds.get("mad_counts_lower") is not None:
        flags |= obs["total_counts"] < _lo(obs["log1p_total_counts"],
                                           thresholds["mad_counts_lower"])
    if thresholds.get("mad_counts_upper") is not None:
        flags |= obs["total_counts"] > _hi(obs["log1p_total_counts"],
                                           thresholds["mad_counts_upper"])
    if thresholds.get("mad_mt") is not None:
        _, hi = mad_bounds(obs["pct_counts_mt"], thresholds["mad_mt"])
        flags |= obs["pct_counts_mt"] > hi

    # Absolute threshold filters
    if thresholds.get("threshold_genes_lower") is not None:
        flags |= obs["n_genes_by_counts"] < thresholds["threshold_genes_lower"]
    if thresholds.get("threshold_genes_upper") is not None:
        flags |= obs["n_genes_by_counts"] > thresholds["threshold_genes_upper"]
    if thresholds.get("threshold_counts_lower") is not None:
        flags |= obs["total_counts"] < thresholds["threshold_counts_lower"]
    if thresholds.get("threshold_counts_upper") is not None:
        flags |= obs["total_counts"] > thresholds["threshold_counts_upper"]
    if thresholds.get("threshold_mt_upper") is not None:
        flags |= obs["pct_counts_mt"] > thresholds["threshold_mt_upper"]

    adata.obs["outlier"] = flags.values
    return adata


def resolve_display_thresholds(adata: sc.AnnData, thresholds: dict) -> dict:
    """
    Derive effective threshold values for histogram/plot annotation.

    When both a MAD-based and an absolute threshold are active for the same
    bound, the more restrictive value is returned (higher lower-bound or lower
    upper-bound).  Returns None for any bound that has no active filter.
    """
    obs = adata.obs

    def _lo(series, m): return float(np.expm1(mad_bounds(series, m)[0]))
    def _hi(series, m): return float(np.expm1(mad_bounds(series, m)[1]))

    def _most_restrictive_lower(*vals):
        vals = [v for v in vals if v is not None]
        return max(vals) if vals else None

    def _most_restrictive_upper(*vals):
        vals = [v for v in vals if v is not None]
        return min(vals) if vals else None

    min_genes = _most_restrictive_lower(
        _lo(obs["log1p_n_genes_by_counts"], thresholds["mad_genes_lower"])
            if thresholds.get("mad_genes_lower") is not None else None,
        thresholds.get("threshold_genes_lower"),
    )
    max_genes = _most_restrictive_upper(
        _hi(obs["log1p_n_genes_by_counts"], thresholds["mad_genes_upper"])
            if thresholds.get("mad_genes_upper") is not None else None,
        thresholds.get("threshold_genes_upper"),
    )
    min_counts = _most_restrictive_lower(
        _lo(obs["log1p_total_counts"], thresholds["mad_counts_lower"])
            if thresholds.get("mad_counts_lower") is not None else None,
        thresholds.get("threshold_counts_lower"),
    )
    max_pct_mt_vals = []
    if thresholds.get("mad_mt") is not None:
        _, hi = mad_bounds(obs["pct_counts_mt"], thresholds["mad_mt"])
        max_pct_mt_vals.append(float(hi))
    if thresholds.get("threshold_mt_upper") is not None:
        max_pct_mt_vals.append(thresholds["threshold_mt_upper"])
    max_pct_mt = _most_restrictive_upper(*max_pct_mt_vals) if max_pct_mt_vals else None

    return {
        "min_genes":  int(min_genes)  if min_genes  is not None else None,
        "max_genes":  int(max_genes)  if max_genes  is not None else None,
        "min_counts": int(min_counts) if min_counts is not None else None,
        "max_pct_mt": max_pct_mt,
    }


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
        "doublet_score",
        "predicted_doublet",
    ]
    outlier_cols = sorted(
        [c for c in adata.obs.columns if c.startswith("outlier")]
    )

    all_cols = [c for c in base_cols + outlier_cols if c in adata.obs.columns]
    df = adata.obs[all_cols].copy()

    df.insert(0, "cell_id",   adata.obs.index)
    df.insert(0, "sample_id", sample_id)

    return df
