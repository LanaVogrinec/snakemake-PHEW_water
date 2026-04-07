rule dehybridize_clusters:
    input:
        merged = RESULTS_DIR + "/merged/test_Apr24_DS1/clusters_with_taxonomy_no_DS1.tsv"
    output:
        final_assign = RESULTS_DIR + "/merged/test_Apr24_DS1/clusters_non-hybrid_no_DS1.tsv",
        progress     = RESULTS_DIR + "/merged/test_Apr24_DS1/dehybridization_summary_no_DS1.tsv"
    log:
        logO = "logs/dehybridize_clusters/test_Apr24_DS1/dehybridize_clusters.log",
        logE = "logs/dehybridize_clusters/test_Apr24_DS1/dehybridize_clusters.err.log"
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
