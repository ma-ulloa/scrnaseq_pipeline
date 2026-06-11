
def sample_raw_h5(wildcards):
    filtered = SAMPLE_PATHS[wildcards.sample]
    return filtered.replace("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5")



# ============================================================
# Rule: soupx
# Ambient RNA correction with SoupX — per sample.
# Input:  filtered + raw .h5 count matrices
# Output: results/soupx/<sample>_soupx.h5ad
# ============================================================

rule soupx:
    input:
        filtered  = sample_h5,
        raw       = sample_raw_h5,
        soupgenes = config["soupx"]["soupgenes"],
        samples  = config["input"]["samples"]
    output:
        h5ad = os.path.join(SOUPX_DIR, "{sample}.h5ad"),
    params:
        sample_id          = "{sample}",
        apply_to_fractions = lambda wc: " ".join(config["soupx"]["apply_to_fractions"]),
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/pipelineR"
    log:
        "logs/soupx/{sample}.log"
    shell:
        """
        Rscript workflow/scripts/01_SoupX.R \
            --filtered           "{input.filtered}"          \
            --raw                "{input.raw}"               \
            --soupgenes          "{input.soupgenes}"         \
            --samples            "{input.samples}"          \
            --sample_id          "{params.sample_id}"        \
            --output             "{output.h5ad}"             \
            --apply_to_fractions {params.apply_to_fractions} \
        > {log} 2>&1
        """
