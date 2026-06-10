
# ============================================================
# Rule: clustering
# Leiden clustering at multiple resolutions + marker gene scoring
# Input:  Harmony-integrated sc.h5ad
# Output: annotated sc.h5ad + figures/
# ============================================================

rule clustering:
    input:
        h5ad        = clustering_input,
        marker_file = config["clustering"]["marker_file"]
    output:
        h5ad = os.path.join(CLUSTERING_DIR, "files", "sc.h5ad")
    params:
        resolutions = lambda wc: " ".join(str(r) for r in config["clustering"]["resolution_test"]),
        seed        = config["clustering"]["seed"],
        basis       = "umap_harmony" if config["clustering"]["use_harmony"] else "umap_pca"
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/scanpy-env/"
    log:
        "workflow/logs/clustering/log.log"
    shell:
        """
        python workflow/scripts/07_clustering.py \
            --input       "{input.h5ad}"          \
            --markers     "{input.marker_file}"   \
            --output      "{output.h5ad}"         \
            --resolutions {params.resolutions}    \
            --seed        {params.seed}           \
            --basis       {params.basis}          \
        > {log} 2>&1
        """
