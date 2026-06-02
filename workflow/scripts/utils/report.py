"""
utils/report.py

HTML and PDF report for the QC step.

Functions:
    build_html_report(...)  — assembles Matplotlib figures into a self-contained HTML
    export_pdf(...)         — renders Matplotlib figures to a multi-page PDF (no kaleido needed!)
"""

import os
import base64
from io import BytesIO

#  HTML 
def build_html_report(
    figures: list,
    sample_id: str,
    stage: str,
    thresholds: dict,
    n_cells: int,
) -> str:
    """
    Assemble Matplotlib figures into a single self-contained HTML string
    by converting them into static embedded base64 PNG tags.

    Parameters
    ----------
    figures     : list of Matplotlib Figure objects
    sample_id   : sample identifier shown in the header
    stage       : 'pre' or 'post' (controls header colour and label)
    thresholds  : dict of QC threshold values shown in the table
    n_cells     : total cell count shown in the info cards
    """
    stage_label = "Pre-Filtering" if stage == "pre" else "Post-Filtering"
    color       = "#2196F3"        if stage == "pre" else "#4CAF50"

    # Convert each Matplotlib figure to a static embedded base64 image tag
    plots_html_list = []
    for fig in figures:
        buf = BytesIO()
        # Save figure to memory buffer
        fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        buf.close()

        # Append structured image block inside your original container styling
        plots_html_list.append(
            f'<div class="plot-container" style="text-align: center;">\n'
            f'    <img src="data:image/png;base64,{img_base64}" style="max-width: 100%; height: auto;" />\n'
            f'</div>'
        )
    
    plots_html = "\n".join(plots_html_list)

    threshold_rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in thresholds.items()
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>QC Report — {sample_id}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
            background: #f5f5f5; margin: 0; padding: 0; color: #333;
        }}
        .header {{ background: {color}; color: white; padding: 30px 40px; }}
        .header h1 {{ margin: 0 0 8px 0; font-size: 1.8em; }}
        .header p  {{ margin: 0; opacity: 0.9; font-size: 1em; }}
        .badge {{
            display: inline-block; background: rgba(255,255,255,0.25);
            border-radius: 20px; padding: 4px 14px;
            font-size: 0.85em; margin-top: 10px;
        }}
        .content {{ max-width: 1400px; margin: 30px auto; padding: 0 20px; }}
        .info-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
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
            margin: 30px 0 16px 0; padding-bottom: 8px;
            border-bottom: 2px solid #eee;
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
        <p>Sample: <strong>{sample_id}</strong>
           &nbsp;|&nbsp;
        <span class="badge">{stage_label}</span>
    </div>

    <div class="content">

        <div class="info-grid">
            <div class="info-card">
                <div class="label">Sample ID</div>
                <div class="value">{sample_id}</div>
            </div>
            <div class="info-card">
                <div class="label">Stage</div>
                <div class="value">{stage_label}</div>
            </div>
            <div class="info-card">
                <div class="label">Total Cells</div>
                <div class="value">{n_cells:,}</div>
            </div>
        </div>

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

    <footer>Generated by scRNA-seq QC Pipeline &nbsp;|&nbsp; Scanpy + Matplotlib</footer>
</body>
</html>"""


#  PDF 

def export_pdf(figures: list, pdf_path: str, sample_id: str, stage: str):
    """
    Export all Matplotlib figures to a single multi-page landscape PDF using ReportLab.
     Bypasses Chromium and Kaleido completely.

    Parameters
    ----------
    figures    : list of Matplotlib Figure objects
    pdf_path   : output path for the PDF file
    sample_id  : used in the document title
    stage      : 'pre' or 'post'
    """
    try:
        from reportlab.lib.pagesizes import A4, landscape
        from reportlab.platypus import SimpleDocTemplate, Image, Spacer, Paragraph
        from reportlab.lib.styles import getSampleStyleSheet
        from reportlab.lib.units import cm
    except ImportError:
        print("[WARN] reportlab not installed — skipping PDF. Install with: pip install reportlab")
        return

    import tempfile

    stage_label = "Pre-Filtering" if stage == "pre" else "Post-Filtering"
    page_w, page_h = landscape(A4)

    doc = SimpleDocTemplate(
        pdf_path,
        pagesize=landscape(A4),
        leftMargin=1.5*cm, rightMargin=1.5*cm,
        topMargin=1.5*cm,  bottomMargin=1.5*cm,
    )

    styles = getSampleStyleSheet()
    story  = [
        Paragraph(f"<b>QC Report — {sample_id} — {stage_label}</b>", styles["Title"]),
        Spacer(1, 0.5*cm),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        # Calculate maximum allowed width inside the margins safely
        # page_w is the full A4 landscape width; 3.5*cm gives a safe buffer zone
        max_img_w = page_w - 3.5 * cm 

        for i, fig in enumerate(figures):
            img_path = os.path.join(tmpdir, f"fig_{i}.png")
            
            # Save the Matplotlib figure to a temporary file locally
            fig.savefig(img_path, format='png', bbox_inches='tight', dpi=150)
            
            # Dynamically grab aspect ratio from the Matplotlib figure size bounds
            fig_w, fig_h = fig.get_size_inches()
            aspect_ratio = fig_h / fig_w
            
            # Set target width to our safe boundary limit
            img_w = max_img_w
            img_h = img_w * aspect_ratio
            
            # OPTIONAL GUARD: If the image height ends up too tall for the page, scale down based on height instead
            max_img_h = page_h - 4.5 * cm
            if img_h > max_img_h:
                img_h = max_img_h
                img_w = img_h / aspect_ratio

            story.append(Image(img_path, width=img_w, height=img_h))
            story.append(Spacer(1, 0.4 * cm))

        doc.build(story)

    print(f"[INFO] PDF report saved to: {pdf_path}")