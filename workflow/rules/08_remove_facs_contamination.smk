
# ============================================================
# Rule: remove_facs_contamination
# Remove cross-fraction FACS contamination from merged object
# Input:  merged .h5ad (all samples concatenated after QC)
# Output: cleaned .h5ad with contaminant cells removed
# ============================================================

rule remove_facs_contamination:
    input:
        h5ad = os.path.join(MERGE_DIR, "files", "sc.h5ad")
    output:
        h5ad = os.path.join(FACS_DIR, "files", "sc.h5ad")
    params:
        species          = config["species"],
        non_sorted_flag  = config["non_sorted_flag"],
        immune_label     = config["immune_label"],
        epidn_label      = config["epidn_label"],
        seed             = config["clustering"]["seed"],
        remove_clusters  = lambda wc: " ".join(str(c) for c in (config.get("remove_clusters") or []))
    conda:
        "../envs/scanpy_env.yml"
    log:
        "logs/remove_facs_contamination/log.log"
    shell:
        """
        python workflow/scripts/07_remove_facs_contamination.py \
            --input           "{input.h5ad}"             \
            --output          "{output.h5ad}"            \
            --species         "{params.species}"         \
            --non_sorted_flag "{params.non_sorted_flag}" \
            --immune_label    "{params.immune_label}"    \
            --epidn_label     "{params.epidn_label}"     \
            --seed            {params.seed}              \
            --remove_clusters {params.remove_clusters}   \
        > {log} 2>&1
        """
