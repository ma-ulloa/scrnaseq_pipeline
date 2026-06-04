# Rule: qc_report_pre
# Generates a QC HTML + PDF report BEFORE filtering.
# Input:  raw .h5 count matrix per sample
# Output: results/qc/pre_<sample>_qc.html
# ============================================================

rule qc_report_pre:
    input:
        h5       = sample_h5, # helper function in Snakemake file
        metadata = config["input"]["metadata"],
        cfg      = "config/config.yaml",
    output:
        html = os.path.join(QC_DIR, "pre_{sample}_qc.html"),
    params:
        sample_id = "{sample}",
        stage     = "pre",
    conda:
        "../envs/scanpy.yaml"
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
