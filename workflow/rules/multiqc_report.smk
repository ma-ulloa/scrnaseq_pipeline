# rules/multiqc_report.smk
# ══════════════════════════════════════════════════════════════════════════════
# Rule: multiqc_summary
# ──────────────────────────────────────────────────────────────────────────────
# Combines per-sample per-cell QC CSVs into a single cohort-level report.
# Runs once for each stage (pre / post) via the {stage} wildcard.
#
# Inputs:
#   csvs     — all per-sample *_cell_qc_metrics.csv for the requested stage
#   metadata — config/metadata.csv (optional, joined for metadata-aware plots)
#
# Outputs:
#   csv      — master merged CSV (one row per cell across all samples)
#   plots    — directory of PNG figures
#   html     — self-contained HTML report
#   pdf      — multi-page PDF report
# ══════════════════════════════════════════════════════════════════════════════

rule multiqc_summary:
    input:
        csvs = lambda wildcards: expand(
            os.path.join(QC_DIR, f"{wildcards.stage}_{{sample}}_cell_qc_metrics.csv"),
            sample=SAMPLES
        ),
        metadata = config["input"]["metadata"],
    output:
        csv   = os.path.join(QC_DIR, "{stage}_cohort_qc_metrics.csv"),
        plots = directory(os.path.join(QC_DIR, "cohort_plots_{stage}")),
        html  = os.path.join(QC_DIR, "{stage}_cohort_qc_report.html"),
        pdf   = os.path.join(QC_DIR, "{stage}_cohort_qc_report.pdf"),
    wildcard_constraints:
        # only fire for valid stage values — prevents accidental wildcard matches
        stage = "pre|post"
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/scanpy-env"
    log:
        "logs/qc/cohort_summary_{stage}.log"
    shell:
        """
        python workflow/scripts/02_multisample_report.py \
            --input_files {input.csvs}        \
            --metadata    {input.metadata}    \
            --output_csv  "{output.csv}"      \
            --plot_dir    "{output.plots}"    \
            --output_html "{output.html}"     \
            --stage       "{wildcards.stage}" \
        > {log} 2>&1
        """
