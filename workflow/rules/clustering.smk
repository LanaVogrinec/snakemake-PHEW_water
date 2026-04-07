rule clustering:
    input:
        fasta = RESULTS_DIR + "/merged/viral_contigs.fasta",
        edges = RESULTS_DIR + "/merged/test_Apr24_DS1/edges_filtered_no_DS1.tsv"
    output:
        clusters = RESULTS_DIR + "/merged/test_Apr24_DS1/clusters_no_DS1.tsv",
        stats    = RESULTS_DIR + "/merged/test_Apr24_DS1/clusters_stats_no_DS1.tsv"
    params:
        resolution = 1.0
    log:
        logO = "logs/clustering/test_Apr24_DS1/clustering.log",
        logE = "logs/clustering/test_Apr24_DS1/clustering.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/clustering.py \
            --fasta {input.fasta} \
            --edges {input.edges} \
            --clusters {output.clusters} \
            --stats {output.stats} \
            --resolution {params.resolution} \
            > {log.logO} 2> {log.logE}
        """
