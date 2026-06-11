rule multiqc_post:
    input:
        csvs     = expand(os.path.join(QC_POST_DIR, "post_{sample}_cell_qc_metrics.csv"),
                          sample=SAMPLES),
        samples  = config["input"]["samples"],
    output:
        csv      = os.path.join(MULTIQC_POST_DIR, "post_cohort_qc_metrics.csv"),
        html     = os.path.join(MULTIQC_POST_DIR, "post_cohort_qc_report.html"),
        pdf      = os.path.join(MULTIQC_POST_DIR, "post_cohort_qc_report.pdf"),
        plot_dir = directory(os.path.join(os.path.dirname(MULTIQC_POST_DIR), "figures")),
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/scanpy-env"
    log:
        "logs/qc/cohort_summary_post.log"
    shell:
        """
        python workflow/scripts/03_multisample_report.py \
            --input_files {input.csvs}         \
            --samples     {input.samples}      \
            --output_csv  {output.csv}         \
            --plot_dir    {output.plot_dir}    \
            --output_html {output.html}        \
            --stage       post                 \
        > {log} 2>&1
        """
