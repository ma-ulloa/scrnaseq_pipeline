#!/usr/bin/env python3

"""
QC report script for single-cell RNA-seq data

Generates an interactive HTML report and pdf with QC plots .
To run BEFORE and AFTER filtering (controlled by --stage).

- Number of genes and transcripts per cell
- Proportion of mitochondrial genes

This script accepts:
Count matrices in the following formats:
- .h5
- .hd5a

Metadata in the following formats:
- .txt, .csv
- First column is always treated as the sample ID index

Usage:
    python 01_qc_report.py \
        --input data/raw/sample.h5 \
        --metadata config/metadata.csv \
        --config config/config.yaml \
        --sample_id sample_id \
        --technology 10x \
        --stage pre \
        --output results/qc/pre_qc.html
"""
import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import yaml

sc.settings.verbosity = 1

# Parsing arguments - to be run in terminal
def parse_args():
    parser = argparse.ArgumentParser(description="Generate QC report for scRNA-seq data")
    parser.add_argument("--input",       required=True,  help="Path to .h5 count matrix")
    parser.add_argument("--metadata",    required=True,  help="Path to metadata CSV")
    parser.add_argument("--config",      required=True,  help="Path to config.yaml")
    parser.add_argument("--sample_id",   required=True,  help="Sample ID (must match metadata)")
    parser.add_argument("--technology",  required=True,  choices=["10x", "other"], help="Input technology")
    parser.add_argument("--stage",       required=True,  choices=["pre", "post"], help="pre or post filtering")
    parser.add_argument("--output",      required=True,  help="Path to output HTML file")
    return parser.parse_args()
 
# Load data

# ── Data Loading ─────────────────────────────────────────────

def load_data(input_path: str, technology: str) -> sc.AnnData:
    """
    Load count matrix depending on technology.
    The input-technology is decided on format file (.h5 or .h5ad)
    """
    if not os.path.exists(input_path):
        sys.exit(f"[ERROR] Input file not found: {input_path}")
 
    if technology == "10x":
        print(f"[INFO] Loading 10x .h5 file: {input_path}")
        adata = sc.read_10x_h5(input_path)
        adata.var_names_make_unique()
 
    elif technology == "other":
        print(f"[INFO] Loading .h5 file: {input_path}")
        adata = sc.read_h5ad(input_path)
        adata.var_names_make_unique()
 
    return adata

    # TO DO: check if necessary

def attach_metadata(adata: sc.AnnData, metadata_path: str, sample_id: str) -> sc.AnnData:
    """
    Attach all metadata columns to adata.obs.
 
    Rules:
      - First column is always the sample ID index 
      - All remaining columns are attached automatically
      - sample_id must match a value in the first column
    """
    if not os.path.exists(metadata_path):
        sys.exit(f"[ERROR] Metadata file not found: {metadata_path}")
 
    meta = pd.read_csv(metadata_path)
 
    if meta.shape[1] < 1:
        sys.exit("[ERROR] Metadata file appears to be empty.")
 
    # First column is the sample ID index regardless of its name
    id_col = meta.columns[0]
    meta[id_col] = meta[id_col].astype(str)
 
    row = meta[meta[id_col] == str(sample_id)]
    if row.empty:
        sys.exit(
            f"[ERROR] sample_id '{sample_id}' not found in metadata column '{id_col}'.\n"
            f"        Available IDs: {meta[id_col].tolist()}"
        )
 
    # Attach all non-index columns dynamically
    extra_cols = [c for c in meta.columns if c != id_col]
    for col in extra_cols:
        adata.obs[col] = str(row[col].values[0])
 
    # Always store sample_id under a consistent key
    adata.obs["sample_id"] = sample_id
 
    print(f"[INFO] Metadata attached for '{sample_id}': {extra_cols}")
    return adata

# QC Metrics 

# Mitochondrial gene prefixes for human and mouse
MT_PREFIXES = (
    "MT-",   # human (Homo sapiens)
    "mt-",   # mouse (Mus musculus) — lowercase
    "MT_",   # some annotation conventions
)
 
# Ribosomal gene prefixes (human + mouse share these)
RIBO_PREFIXES = ("RPS", "RPL", "Rps", "Rpl")
 

def compute_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    """
    Compute standard QC metrics.
 
    Mitochondrial detection is flexible:
      - Human: MT- prefix (e.g. MT-CO1, MT-ND1)
      - Mouse:  mt- prefix (e.g. mt-Co1, mt-Nd1)
    Both are detected automatically.
    """
    # Mitochondrial — check human AND mouse prefixes
    adata.var["mt"] = adata.var_names.str.startswith(MT_PREFIXES)
    n_mt = adata.var["mt"].sum()
    print(f"[INFO] Mitochondrial genes detected: {n_mt}")
    if n_mt == 0:
        print("[WARN] No mitochondrial genes found. Check gene naming convention.")
 
    # Ribosomal
    adata.var["ribo"] = adata.var_names.str.startswith(RIBO_PREFIXES)
    n_ribo = adata.var["ribo"].sum()
    print(f"[INFO] Ribosomal genes detected: {n_ribo}")
 
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    print(f"[INFO] QC metrics computed. Cells: {adata.n_obs}, Genes: {adata.n_vars}")
    return adata

#  Plotly Figures

COLORS = px.colors.qualitative.Set2
 
 
def fig_violin(adata: sc.AnnData, sample_id: str) -> go.Figure:
    """Violin plots for n_genes, total_counts, pct_counts_mt, pct_counts_ribo."""
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
            row=1, col=i
        )
 
    fig.update_layout(
        title_text=f"QC Metrics — {sample_id}",
        height=450,
        template="plotly_white"
    )
    return fig
 
 
def fig_scatter(adata: sc.AnnData, sample_id: str) -> go.Figure:
    """Scatter: total_counts vs n_genes, colored by pct_counts_mt."""
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
            "pct_counts_mt":     "% MT"
        },
        title=f"Counts vs Genes (colored by % MT) — {sample_id}",
        opacity=0.5,
        template="plotly_white"
    )
    fig.update_layout(height=450)
    return fig
 
 
def fig_histogram_mt(adata: sc.AnnData, sample_id: str, max_pct_mt: float) -> go.Figure:
    """Histogram of % MT counts with threshold line."""
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=adata.obs["pct_counts_mt"],
        nbinsx=50,
        marker_color=COLORS[2],
        opacity=0.8,
        name="% MT"
    ))
    fig.add_vline(
        x=max_pct_mt,
        line_dash="dash", line_color="red",
        annotation_text=f"Threshold: {max_pct_mt}%",
        annotation_position="top right"
    )
    fig.update_layout(
        title=f"% MT Counts Distribution — {sample_id}",
        xaxis_title="% MT Counts",
        yaxis_title="Cell Count",
        template="plotly_white",
        height=400
    )
    return fig
 
 
def fig_histogram_genes(adata: sc.AnnData, sample_id: str,
                         min_genes: int, max_genes: int) -> go.Figure:
    """Histogram of n_genes with threshold lines."""
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=adata.obs["n_genes_by_counts"],
        nbinsx=50,
        marker_color=COLORS[0],
        opacity=0.8,
        name="N Genes"
    ))
    fig.add_vline(x=min_genes, line_dash="dash", line_color="red",
                  annotation_text=f"Min: {min_genes}", annotation_position="top left")
    fig.add_vline(x=max_genes, line_dash="dash", line_color="red",
                  annotation_text=f"Max: {max_genes}", annotation_position="top right")
    fig.update_layout(
        title=f"N Genes Distribution — {sample_id}",
        xaxis_title="N Genes",
        yaxis_title="Cell Count",
        template="plotly_white",
        height=400
    )
    return fig
 
 
def fig_summary_table(adata: sc.AnnData, sample_id: str, thresholds: dict) -> go.Figure:
    """Summary statistics table."""
    obs = adata.obs
    stats = {
        "Metric": [
            "Total cells",
            "Median genes/cell",
            "Median counts/cell",
            "Median % MT",
            "Median % Ribo",
            "Cells > max_genes threshold",
            "Cells < min_genes threshold",
            "Cells > max_pct_mt threshold",
        ],
        "Value": [
            f"{adata.n_obs:,}",
            f"{obs['n_genes_by_counts'].median():.0f}",
            f"{obs['total_counts'].median():.0f}",
            f"{obs['pct_counts_mt'].median():.2f}%",
            f"{obs['pct_counts_ribo'].median():.2f}%",
            f"{(obs['n_genes_by_counts'] > thresholds['max_genes']).sum():,}",
            f"{(obs['n_genes_by_counts'] < thresholds['min_genes']).sum():,}",
            f"{(obs['pct_counts_mt'] > thresholds['max_pct_mt']).sum():,}",
        ]
    }
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=["<b>Metric</b>", "<b>Value</b>"],
            fill_color="#4C72B0",
            font=dict(color="white", size=13),
            align="left"
        ),
        cells=dict(
            values=[stats["Metric"], stats["Value"]],
            fill_color=[["#f9f9f9", "white"] * 10],
            align="left",
            font=dict(size=12)
        )
    )])
    fig.update_layout(
        title=f"Summary Statistics — {sample_id}",
        height=350,
        margin=dict(t=50, b=10)
    )
    return fig
 
 
# ── PDF Export ───────────────────────────────────────────────
 
def export_pdf(figures: list, pdf_path: str, sample_id: str, stage: str):
    """
    Export all figures to a single PDF (one figure per page).
    Requires kaleido: pip install kaleido
    """
    try:
        import kaleido  # noqa: F401 — just to check availability
    except ImportError:
        print("[WARN] kaleido not installed — skipping PDF export. Install with: pip install kaleido")
        return
 
    import tempfile
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.platypus import SimpleDocTemplate, Image, Spacer, Paragraph
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.units import cm
 
    stage_label = "Pre-Filtering" if stage == "pre" else "Post-Filtering"
    page_w, page_h = landscape(A4)
 
    doc = SimpleDocTemplate(
        pdf_path,
        pagesize=landscape(A4),
        leftMargin=1.5*cm, rightMargin=1.5*cm,
        topMargin=1.5*cm, bottomMargin=1.5*cm
    )
 
    styles = getSampleStyleSheet()
    story = []
 
    # Title page element
    story.append(Paragraph(
        f"<b>QC Report — {sample_id} — {stage_label}</b>",
        styles["Title"]
    ))
    story.append(Spacer(1, 0.5*cm))
 
    with tempfile.TemporaryDirectory() as tmpdir:
        for i, fig in enumerate(figures):
            img_path = os.path.join(tmpdir, f"fig_{i}.png")
            fig.write_image(img_path, width=1100, height=500, scale=2)
            img_w  = page_w - 3*cm
            img_h  = img_w * 500 / 1100
            story.append(Image(img_path, width=img_w, height=img_h))
            story.append(Spacer(1, 0.4*cm))
 
        doc.build(story)
 
    print(f"[INFO] PDF report saved to: {pdf_path}")
 
 
# ── HTML Report ──────────────────────────────────────────────
 
def build_html_report(figures: list, sample_id: str, stage: str,
                      thresholds: dict, technology: str,
                      n_cells: int, metadata_cols: list) -> str:
    """Assemble all figures into a single self-contained HTML string."""
 
    stage_label = "Pre-Filtering" if stage == "pre" else "Post-Filtering"
    color       = "#2196F3" if stage == "pre" else "#4CAF50"
 
    plots_html = "\n".join([
        f'<div class="plot-container">{fig.to_html(full_html=False, include_plotlyjs=False)}</div>'
        for fig in figures
    ])
 
    threshold_rows = "".join([
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in thresholds.items()
    ])
 
    # Info cards: sample_id + all metadata columns
    info_cards = f"""
        <div class="info-card"><div class="label">Sample ID</div><div class="value">{sample_id}</div></div>
        <div class="info-card"><div class="label">Stage</div><div class="value">{stage_label}</div></div>
        <div class="info-card"><div class="label">Total Cells</div><div class="value">{n_cells:,}</div></div>
        <div class="info-card"><div class="label">Technology</div><div class="value">{technology.upper()}</div></div>
    """
 
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>QC Report — {sample_id}</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
            background: #f5f5f5; margin: 0; padding: 0; color: #333;
        }}
        .header {{
            background: {color}; color: white; padding: 30px 40px;
        }}
        .header h1 {{ margin: 0 0 8px 0; font-size: 1.8em; }}
        .header p  {{ margin: 0; opacity: 0.9; font-size: 1em; }}
        .badge {{
            display: inline-block; background: rgba(255,255,255,0.25);
            border-radius: 20px; padding: 4px 14px; font-size: 0.85em; margin-top: 10px;
        }}
        .content {{ max-width: 1400px; margin: 30px auto; padding: 0 20px; }}
        .info-grid {{
            display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 16px; margin-bottom: 30px;
        }}
        .info-card {{
            background: white; border-radius: 10px; padding: 18px 22px;
            box-shadow: 0 1px 4px rgba(0,0,0,0.08);
        }}
        .info-card .label {{ font-size: 0.75em; color: #888; text-transform: uppercase; }}
        .info-card .value {{ font-size: 1.4em; font-weight: 600; color: #333; margin-top: 4px; }}
        .plot-container {{
            background: white; border-radius: 10px; padding: 20px;
            margin-bottom: 24px; box-shadow: 0 1px 4px rgba(0,0,0,0.08);
        }}
        .section-title {{
            font-size: 1.1em; font-weight: 600; color: #555;
            margin: 30px 0 16px 0; padding-bottom: 8px; border-bottom: 2px solid #eee;
        }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ padding: 8px 14px; text-align: left; border-bottom: 1px solid #eee; }}
        th {{ background: #f0f0f0; font-weight: 600; }}
        .threshold-table {{
            background: white; border-radius: 10px; padding: 20px;
            margin-bottom: 24px; box-shadow: 0 1px 4px rgba(0,0,0,0.08);
        }}
        footer {{ text-align: center; padding: 30px; color: #aaa; font-size: 0.85em; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🔬 Single-Cell QC Report</h1>
        <p>Sample: <strong>{sample_id}</strong> &nbsp;|&nbsp; Technology: <strong>{technology.upper()}</strong></p>
        <span class="badge">{stage_label}</span>
    </div>
 
    <div class="content">
        <div class="info-grid">{info_cards}</div>
 
        <div class="section-title">📋 QC Thresholds Applied</div>
        <div class="threshold-table">
            <table>
                <tr><th>Parameter</th><th>Value</th></tr>
                {threshold_rows}
            </table>
        </div>
 
        <div class="section-title">📊 QC Plots</div>
        {plots_html}
    </div>
 
    <footer>Generated by scRNA-seq QC Pipeline &nbsp;|&nbsp; Scanpy + Plotly</footer>
</body>
</html>"""
    return html
 
 
# ── Main ─────────────────────────────────────────────────────
 
def main():
    args = parse_args()
 
    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)
    thresholds = config["qc_thresholds"]
 
    # Load data & attach metadata
    adata = load_data(args.input, args.technology)
    adata = attach_metadata(adata, args.metadata, args.sample_id)
    adata = compute_qc_metrics(adata)
 
    # Collect metadata column names for display
    meta      = pd.read_csv(args.metadata)
    id_col    = meta.columns[0]
    meta_cols = [c for c in meta.columns if c != id_col]
 
    n_cells = adata.n_obs
 
    # Build figures
    figures = [
        fig_summary_table(adata, args.sample_id, thresholds),
        fig_violin(adata, args.sample_id),
        fig_scatter(adata, args.sample_id),
        fig_histogram_genes(adata, args.sample_id,
                            thresholds["min_genes"], thresholds["max_genes"]),
        fig_histogram_mt(adata, args.sample_id, thresholds["max_pct_mt"]),
    ]
 
    # ── HTML output
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    html = build_html_report(
        figures=figures,
        sample_id=args.sample_id,
        stage=args.stage,
        thresholds=thresholds,
        technology=args.technology,
        n_cells=n_cells,
        metadata_cols=meta_cols
    )
    with open(args.output, "w") as f:
        f.write(html)
    print(f"[INFO] HTML report saved to: {args.output}")
 
    # ── PDF output (same path, .pdf extension)
    pdf_path = os.path.splitext(args.output)[0] + ".pdf"
    export_pdf(figures, pdf_path, args.sample_id, args.stage)
 
 
if __name__ == "__main__":
    main()