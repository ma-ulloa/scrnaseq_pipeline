rule multiqc_pre:
    input:
        csvs     = expand(os.path.join(QC_PRE_DIR, "pre_{sample}_cell_qc_metrics.csv"),
                          sample=SAMPLES),
        samples  = config["input"]["samples"],
    output:
        csv      = os.path.join(MULTIQC_PRE_DIR, "pre_cohort_qc_metrics.csv"),
        html     = os.path.join(MULTIQC_PRE_DIR, "pre_cohort_qc_report.html"),
        pdf      = os.path.join(MULTIQC_PRE_DIR, "pre_cohort_qc_report.pdf"),
        plot_dir = directory(os.path.join(os.path.dirname(MULTIQC_PRE_DIR), "figures")),
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/scanpy-env"
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
        > {log} 2>&1
        """
