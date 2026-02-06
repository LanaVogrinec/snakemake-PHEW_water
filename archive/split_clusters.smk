rule split_clusters:
    input:
        clusters_with_tax = RESULTS_DIR + "/merged/clusters_with_taxonomy.tsv"
    output:
        wo_singletons = RESULTS_DIR + "/analysis/clusters_with_taxonomy_wo_singletons.tsv",
        singletons    = RESULTS_DIR + "/analysis/clusters_with_taxonomy_singletons.tsv"
    log:
        logO = "logs/split_clusters/split_clusters.log",
        logE = "logs/split_clusters/split_clusters.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/split_clusters.py \
            --input {input.clusters_with_tax} \
            --wo_singletons {output.wo_singletons} \
            --singletons {output.singletons} \
            > {log.logO} 2> {log.logE}
        """
