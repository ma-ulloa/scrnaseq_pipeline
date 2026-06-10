
# ============================================================
# Rule: harmony
# Sample integration with harmonypy
# Input:  merged object
# Output: integrated .h5ad + plots PDF
# ============================================================

rule harmony:
    input:
        h5ad = os.path.join(MERGE_DIR, "sc.h5ad")
    output:
        h5ad  = os.path.join(HARMONY_DIR, "sc.h5ad"),
        plots = os.path.join(HARMONY_DIR, "sc_harmony.pdf")
    params:
        integration_keys = lambda wc: " ".join(config["harmony"]["integration_keys"]),
        n_pcs            = config["harmony"]["n_pcs"]
    conda:
        "/srv/data/users/shared_conda_alejandra_martin/scanpy-env/"
    log:
        "logs/harmony/log.log"
    shell:
        """
        python workflow/scripts/06_harmony.py \
            --input            "{input.h5ad}"             \
            --config           "config/config.yaml"       \
            --integration_keys {params.integration_keys}  \
            --n_pcs            {params.n_pcs}             \
            --output           "{output.h5ad}"            \
            --plots            "{output.plots}"           \
        > {log} 2>&1
        """
