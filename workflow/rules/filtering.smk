rule filter_cells:
    input:
        h5     = qc_input,                                                    # Raw data matrix
        csv    = os.path.join(QC_DIR, "pre_{sample}_cell_qc_metrics.csv"),   # MAD info
        meta   = config["input"]["metadata"],
        cfg    = "config/config.yaml"
    output:
        h5ad   = os.path.join(FILT_DIR, "{sample}_filtered.h5ad")
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/scanpy-env"
    log:
        "logs/filtering/{sample}.log"
    shell:
        """
        python workflow/scripts/02_filtering.py \
            --input     "{input.h5}" \
            --metadata  "{input.meta}" \
            --qc_csv    "{input.csv}" \
            --config    "{input.cfg}" \
            --sample_id "{wildcards.sample}" \
            --output    "{output.h5ad}" \
        > {log} 2>&1
        """