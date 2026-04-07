rule merge_taxonomy_clusters:
    input:
        clusters = RESULTS_DIR + "/merged/test_Apr24_DS1/clusters_no_DS1.tsv",
        taxonomy = RESULTS_DIR + "/merged/merged_taxonomy.tsv"
    output:
        merged = RESULTS_DIR + "/merged/test_Apr24_DS1/clusters_with_taxonomy_no_DS1.tsv"
    log:
        logO = "logs/merge_taxonomy_clusters/test_Apr24_DS1/merge_taxonomy_clusters.log",
        logE = "logs/merge_taxonomy_clusters/test_Apr24_DS1/merge_taxonomy_clusters.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/merge_taxonomy_clusters.py \
          --clusters {input.clusters} \
          --taxonomy {input.taxonomy} \
          --output {output.merged} \
          > {log.logO} 2> {log.logE}
        """
