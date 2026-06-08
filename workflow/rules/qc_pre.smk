# Rule: qc_report_pre
# Generates a QC HTML + PDF report BEFORE filtering.
# Input:  per-sample count matrix — SoupX-corrected .h5ad if
#         config["soupx"]["use_soupx"] is true, otherwise the
#         original raw/filtered .h5 matrix
# Output: results/qc/pre_<sample>_qc.html
# ============================================================

rule qc_report_pre:
    input:
        h5       = qc_input, # helper function in Snakemake file
        metadata = config["input"]["metadata"],
        cfg      = "config/config.yaml",
    output:
        html = os.path.join(QC_DIR, "pre_{sample}_qc.html"),
        pdf  = os.path.join(QC_DIR, "pre_{sample}_qc.pdf"),
        csv  = os.path.join(QC_DIR, "pre_{sample}_cell_qc_metrics.csv")
    params:
        sample_id = "{sample}",
        stage     = "pre",
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/scanpy-env"
    log:
        "logs/qc/pre_{sample}.log"
    shell:
        """
        python workflow/scripts/01_qc_stats_report.py \
            --input     "{input.h5}"         \
            --metadata  {input.metadata}   \
            --cfg    "{input.cfg}"        \
            --sample_id {params.sample_id} \
            --stage     "{params.stage}"     \
            --output    "{output.html}"      \
        > {log} 2>&1
        """
