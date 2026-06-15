rule multiqc_pre:
    input:
        csvs     = expand(os.path.join(QC_PRE_DIR, "files", "pre_{sample}_cell_qc_metrics.csv"),
                          sample=SAMPLES),
        samples  = config["input"]["samples"],
    output:
        csv      = os.path.join(MULTIQC_PRE_DIR, "files", "pre_cohort_qc_metrics.csv"),
        html     = os.path.join(MULTIQC_PRE_DIR, "files", "pre_cohort_qc_report.html"),
        pdf      = os.path.join(MULTIQC_PRE_DIR, "files", "pre_cohort_qc_report.pdf"),
        plot_dir = directory(os.path.join(MULTIQC_PRE_DIR, "figures")),
    params:
        color_by = lambda wc: " ".join(config["qc_plots"]["color_by"]),
    conda:
        "../envs/scanpy_env.yml"
    log:
        "logs/qc/cohort_summary_pre.log"
    shell:
        """
        python workflow/scripts/03_multisample_report.py \
            --input_files {input.csvs}         \
            --samples     {input.samples}      \
            --output_csv  {output.csv}         \
            --plot_dir    {output.plot_dir}    \
            --output_html {output.html}        \
            --stage       pre                  \
            --color_by    {params.color_by}    \
        > {log} 2>&1
        """
