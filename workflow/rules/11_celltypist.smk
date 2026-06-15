
# ============================================================
# Rule: celltypist
# Automated cell type annotation with CellTypist majority voting
# Input:  clustered sc.h5ad from step 09
# Output: annotated sc.h5ad + UMAP figures
# ============================================================


rule celltypist:
    input:
        h5ad = os.path.join(CLUSTERING_DIR, "files", "sc.h5ad")
    output:
        h5ad = os.path.join(CELLTYPIST_DIR, "files", "sc.h5ad")
    params:
        model      = CELLTYPIST_MODEL,
        leiden_key = config["celltypist"].get("leiden_key") or "",
        basis      = "umap_harmony" if config["clustering"]["use_harmony"] else "umap_pca"
    conda:
        "../envs/scanpy_env.yml"
    log:
        "logs/celltypist/log.log"
    shell:
        """
        python workflow/scripts/10_celltypist.py \
            --input      "{input.h5ad}"      \
            --output     "{output.h5ad}"     \
            --model      "{params.model}"    \
            --leiden_key "{params.leiden_key}" \
            --basis      {params.basis}      \
        > {log} 2>&1
        """
