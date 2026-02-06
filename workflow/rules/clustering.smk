rule clustering:
    input:
        fasta = RESULTS_DIR + "/merged/viral_contigs.fasta",
        edges = RESULTS_DIR + "/merged/edges_noApr24_ANI90_AF85.tsv"
    output:
        clusters = RESULTS_DIR + "/merged/leiden_clusters_noApr24_ANI90_AF85.tsv",
        stats    = RESULTS_DIR + "/merged/leiden_stats_noApr24_ANI90_AF85.tsv"
    params:
        resolution = 1.0
    log:
        logO = "logs/leiden_clustering/merged.log",
        logE = "logs/leiden_clustering/merged.err.log"
    conda:
        "../envs/clustering_env.yaml"
    threads: 1
    shell:
        """
        python workflow/scripts/clustering_script.py \
            --fasta {input.fasta} \
            --edges {input.edges} \
            --clusters {output.clusters} \
            --stats {output.stats} \
            --resolution {params.resolution} \
            > {log.logO} 2> {log.logE}
        """
