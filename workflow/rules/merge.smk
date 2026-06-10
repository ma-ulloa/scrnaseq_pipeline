
# ============================================================
# Rule: merge_adatas
# Ambient Merging adatas after QC
# Input:  .h5ad files
# Output: One .h5ad file containing all samples
# ============================================================

rule merge_adatas:
    input:
        h5ads = expand(os.path.join(FILT_DIR, "{sample}_filtered.h5ad"),sample=SAMPLES)
    output:
        h5ad = os.path.join(MERGE_DIR, "sc.h5ad")
    params:
        join = config["merge"]["join"]
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/scanpy-env/"
    log:
        "logs/merge_adatas/log.log"
    shell:
        """
        python workflow/scripts/06_merge_adatas.py \
            --inputs            {input.h5ads}          \
            --join             "{params.join}" \
            --output           "{output.h5ad}"
        > {log} 2>&1
        """
