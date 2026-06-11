# ============================================================
# Rule: filtering
# Filters a single sample using the outlier flags written by qc_report_pre.
# Reads outlier_{n_mads}mad from the pre-filter metrics CSV —
# no MAD recomputation inside this script.
# ============================================================
rule filtering:
    input:
        h5       = qc_input,                    # same input as qc_report_pre
        metrics  = os.path.join(QC_DIR, "pre_{sample}_cell_qc_metrics.csv"),
        samples = config["input"]["samples"],
        cfg      = "config/config.yaml",
    output:
        h5ad = os.path.join(FILT_DIR, "{sample}_filtered.h5ad"),
    params:
        sample_id = "{sample}",
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/scanpy-env"
    log:
        "logs/filtering/{sample}.log"
    shell:
        """
        python workflow/scripts/04_filtering.py \
            --input     "{input.h5}"           \
            --metrics   "{input.metrics}"      \
            --samples  {input.samples}       \
            --cfg       "{input.cfg}"          \
            --sample_id {params.sample_id}     \
            --output    "{output.h5ad}"        \
        > {log} 2>&1
        """
