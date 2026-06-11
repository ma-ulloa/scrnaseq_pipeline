# ============================================================
# Rule: qc_report_post
# Generates a QC HTML, PDF, and cell-metrics CSV report AFTER filtering.
# Input is the filtered .h5ad (not raw counts).
# No MAD sensitivity section — just confirmation of what survived.
# ============================================================
rule qc_report_post:
    input:
        h5ad     = os.path.join(FILT_DIR, "{sample}_filtered.h5ad"),
        samples = config["input"]["samples"],
        cfg      = "config/config.yaml",
    output:
        html = os.path.join(QC_POST_DIR, "post_{sample}_qc.html"),
        pdf  = os.path.join(QC_POST_DIR, "post_{sample}_qc.pdf"),
        csv  = os.path.join(QC_POST_DIR, "post_{sample}_cell_qc_metrics.csv"),
    params:
        sample_id = "{sample}",
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/scanpy-env"
    log:
        "logs/qc/post_{sample}.log"
    shell:
        """
        python workflow/scripts/05_qc_post_report.py \
            --input     "{input.h5ad}"         \
            --samples  {input.samples}       \
            --cfg       "{input.cfg}"          \
            --sample_id {params.sample_id}     \
            --output    "{output.html}"        \
        > {log} 2>&1
        """
