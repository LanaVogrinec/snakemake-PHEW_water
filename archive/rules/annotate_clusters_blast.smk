rule annotate_clusters_blast:
    input:
        blast    = RESULTS_DIR + "/merged/filtered_blastn_pairs.tsv",
        clusters = RESULTS_DIR + "/merged/clusters.tsv"
    output:
        intra = RESULTS_DIR + "/analysis/blast_pairs_intra_clusters.tsv"
    log:
        logO = "logs/annotate_clusters_blast/annotate_clusters_blast.log",
        logE = "logs/annotate_clusters_blast/annotate_clusters_blast.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/annotate_clusters_blast.py \
            --blast {input.blast} \
            --clusters {input.clusters} \
            --out {output.intra} \
            > {log.logO} 2> {log.logE}
        """
