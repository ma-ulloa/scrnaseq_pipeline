# ============================================================
# Rule: qc_report_post
# Generates a QC HTML + PDF report AFTER filtering.
# Input:  filtered .h5ad per sample (output of filter_cells)
# Output: results/qc/post_<sample>_qc.html
# ============================================================

rule qc_report_post:
    input:
        h5ad     = os.path.join(FILT_DIR, "{sample}_filtered.h5ad"),
        metadata = config["input"]["metadata"],
        cfg      = "config/config.yaml",
    output:
        html = os.path.join(QC_DIR, "post_{sample}_qc.html"),
    params:
        sample_id = "{sample}",
        stage     = "post",
    conda:
        "../envs/scanpy.yaml"
    log:
        "logs/qc/post_{sample}.log"
    shell:
        """
        python workflow/scripts/01_qc_stats_report.py \
            --input     "{input.h5ad}"       \
            --metadata  "{input.metadata}"   \
            --config    "{input.cfg}"        \
            --sample_id "{params.sample_id}" \
            --stage     "{params.stage}"     \
            --output    "{output.html}"      \
        > {log} 2>&1
        """
