# Rule: filter_cells
# Applies QC thresholds from config.yaml and saves filtered AnnData.
# Input:  per-sample count matrix — SoupX-corrected .h5ad if
#         config["soupx"]["use_soupx"] is true, otherwise the
#         original raw/filtered .h5 matrix
# Output: results/filtered/<sample>_filtered.h5ad
# ============================================================

rule filter_cells:
    input:
        h5       = qc_input,
        metadata = config["input"]["metadata"],
        cfg      = "config/config.yaml",
    output:
        h5ad = os.path.join(FILT_DIR, "{sample}_filtered.h5ad"),
    params:
        sample_id = "{sample}",
    conda:
        "../envs/scanpy.yaml"
    log:
        "logs/filtering/{sample}.log"
    shell:
        """
        python workflow/scripts/02_filtering.py \
            --input     "{input.h5}"         \
            --metadata  "{input.metadata}"   \
            --config    "{input.cfg}"        \
            --sample_id "{params.sample_id}" \
            --output    "{output.h5ad}"      \
        > {log} 2>&1
        """
