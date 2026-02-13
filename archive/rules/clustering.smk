rule clustering:
    input:
        f = RESULTS_DIR + "/merged/viral_contigs.fasta",
        b = RESULTS_DIR + "/merged/blastn_pairs.tsv"
    output:
        edges    = RESULTS_DIR + "/merged/edges.txt",
        filtered = RESULTS_DIR + "/merged/filtered_blastn_pairs.tsv",
        clusters = RESULTS_DIR + "/merged/clusters.tsv"
    params:
        min_identity = 90.0,
        min_length = 150
    log:
        logO = "logs/clustering/merged.log",
        logE = "logs/clustering/merged.err.log"
    conda:
        "../envs/clustering_env.yaml"
    threads: 30
    shell:
        """
        python workflow/scripts/clustering_script.py \
            --fasta {input.f} \
            --blast {input.b} \
            --edges {output.edges} \
            --filtered {output.filtered} \
            --clusters {output.clusters} \
            --min_identity {params.min_identity} \
            --min_length {params.min_length} \
            > {log.logO} 2> {log.logE}
        """
