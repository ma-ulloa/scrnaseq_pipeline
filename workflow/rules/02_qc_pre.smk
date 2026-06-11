# ============================================================
# Rule: qc_report_pre
# Generates a QC HTML + PDF report and per-cell metrics CSV BEFORE filtering.
# Input:  per-sample count matrix — SoupX-corrected .h5ad if
#         config["soupx"]["use_soupx"] is true, otherwise the
#         original raw/filtered .h5 matrix
# Output: results/02_qc_pre/files/pre_<sample>_qc.html
#         results/02_qc_pre/files/pre_<sample>_qc.pdf
#         results/02_qc_pre/files/pre_<sample>_cell_qc_metrics.csv  ← includes outlier_* cols
# ============================================================
rule qc_report_pre:
    input:
        h5       = qc_input,                    # helper function in Snakefile
        samples = config["input"]["samples"],
        cfg      = "config/config.yaml",
    output:
        html = os.path.join(QC_PRE_DIR, "pre_{sample}_qc.html"),
        pdf  = os.path.join(QC_PRE_DIR, "pre_{sample}_qc.pdf"),
        csv  = os.path.join(QC_PRE_DIR, "pre_{sample}_cell_qc_metrics.csv"),
    params:
        sample_id = "{sample}",
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/scanpy-env"
    log:
        "logs/qc/pre_{sample}.log"
    shell:
        """
        python workflow/scripts/02_qc_stats_report.py \
            --input     "{input.h5}"           \
            --samples  {input.samples}       \
            --cfg       "{input.cfg}"          \
            --sample_id {params.sample_id}     \
            --output    "{output.html}"        \
        > {log} 2>&1
        """
