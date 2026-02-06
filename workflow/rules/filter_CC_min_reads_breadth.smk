rule filter_CC_min_reads_breadth:
    input:
        samples = RESULTS_DIR + "/merged/rep_contigs_minlen500_read_counts_taxonomy_samples_noApr24_ANI90_AF85.tsv",
        carryovers = RESULTS_DIR + "/merged/rep_contigs_minlen500_read_counts_taxonomy_carryovers_noApr24_ANI90_AF85.tsv",
        breadth = RESULTS_DIR + "/merged/rep_contigs_minlen500_breadth_noApr24_ANI90_AF85.tsv"
    output:
        filtered = RESULTS_DIR + "/merged/rep_contigs_minlen500_read_counts_noApr24_ANI90_AF85_filtered.tsv",
        stats    = RESULTS_DIR + "/merged/rep_contigs_minlen500_read_counts_noApr24_ANI90_AF85_filtered_stats.tsv",
        check    = RESULTS_DIR + "/merged/rep_contigs_minlen500_read_counts_noApr24_ANI90_AF85_filtered.done"
    params:
        min_reads = 2,
        min_breadth = 0.25
    log:
        logO = "logs/filter_CC_min_reads_breadth.log",
        logE = "logs/filter_CC_min_reads_breadth.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/filter_CC_min_reads_breadth.py \
            --samples {input.samples} \
            --carryovers {input.carryovers} \
            --breadth {input.breadth} \
            --out-filtered {output.filtered} \
            --out-stats {output.stats} \
            --min-reads {params.min_reads} \
            --min-breadth {params.min_breadth} \
            > {log.logO} 2> {log.logE}

        touch {output.check}
        """
