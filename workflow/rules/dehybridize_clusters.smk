rule dehybridize_clusters:
    input:
        merged = RESULTS_DIR + "/merged/clusters_with_taxonomy.tsv"
    output:
        final_assign = RESULTS_DIR + "/merged/clusters_non-hybrid.tsv",
        progress     = RESULTS_DIR + "/merged/dehybridization_summary.tsv"
    log:
        logO = "logs/dehybridize_clusters/dehybridize_clusters.log",
        logE = "logs/dehybridize_clusters/dehybridize_clusters.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/dehybridize_clusters.py \
          --input {input.merged} \
          --final_assign {output.final_assign} \
          --progress {output.progress} \
          > {log.logO} 2> {log.logE}
        """
