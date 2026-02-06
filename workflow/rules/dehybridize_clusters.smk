rule dehybridize_clusters:
    input:
        merged = RESULTS_DIR + "/merged/clusters_with_taxonomy_noApr24_ANI90_AF85.tsv"
    output:
        final_assign = RESULTS_DIR + "/merged/final_cluster_assignments_non-hybrid_noApr24_ANI90_AF85.tsv",
        progress     = RESULTS_DIR + "/merged/dehybridization_progress_summary_noApr24_ANI90_AF85.tsv"
    log:
        logO = "logs/dehybridize_clusters/merged.log",
        logE = "logs/dehybridize_clusters/merged.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/dehybridize_clusters_script.py \
          --input {input.merged} \
          --final_assign {output.final_assign} \
          --progress {output.progress} \
          > {log.logO} 2> {log.logE}
        """
