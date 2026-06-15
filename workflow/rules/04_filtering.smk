# ============================================================
# Rule: filtering
# Filters a single sample using the outlier flags written by qc_report_pre.
# Reads outlier_{n_mads}mad from the pre-filter metrics CSV —
# no MAD recomputation inside this script.
# ============================================================
rule filtering:
    input:
        h5      = qc_input,                    # same input as qc_report_pre
        metrics = os.path.join(QC_PRE_DIR, "files", "pre_{sample}_cell_qc_metrics.csv"),
        samples = config["input"]["samples"],
    output:
        h5ad = os.path.join(FILT_DIR, "files", "{sample}_filtered.h5ad"),
    params:
        sample_id       = "{sample}",
        remove_doublets = config["doublets"]["remove_doublets"],
        remove_mt       = config["filter_genes"]["remove_mt"],
        remove_ribo     = config["filter_genes"]["remove_ribo"],
        remove_malat1   = config["filter_genes"]["remove_malat1"],
        species         = config["species"],
    conda:
        "../envs/scanpy_env.yml"
    log:
        "logs/filtering/{sample}.log"
    shell:
        """
        python workflow/scripts/04_filtering.py \
            --input           "{input.h5}"                \
            --metrics         "{input.metrics}"           \
            --samples         "{input.samples}"           \
            --sample_id       {params.sample_id}          \
            --remove_doublets {params.remove_doublets}    \
            --remove_mt       {params.remove_mt}          \
            --remove_ribo     {params.remove_ribo}        \
            --remove_malat1   {params.remove_malat1}      \
            --species         {params.species}            \
            --output          "{output.h5ad}"             \
        > {log} 2>&1
        """
