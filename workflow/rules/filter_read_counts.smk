rule filter_read_counts:
    input:
        samples = RESULTS_DIR + "/merged/reps_read_counts_samples.tsv",
        carryover = RESULTS_DIR + "/merged/reps_read_counts_carry-over.tsv",
        breadth = RESULTS_DIR + "/merged/merged_breadth.tsv",
        mapping = RESULTS_DIR + "/merged/carry-over_matches.tsv"
    output:
        filtered = RESULTS_DIR + "/merged/reps_read_counts_samples_filtered.tsv",
        stats    = RESULTS_DIR + "/merged/reps_read_counts_samples_filtered_stats.tsv"
    params:
        min_reads = 2,
        min_breadth = 0.25
    log:
        logO = "logs/filter_read_counts/filter_read_counts.log",
        logE = "logs/filter_read_counts/filter_read_counts.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/filter_read_counts.py \
            --samples {input.samples} \
            --carryover {input.carryover} \
            --breadth {input.breadth} \
            --mapping {input.mapping} \
            --out-filtered {output.filtered} \
            --out-stats {output.stats} \
            --min-reads {params.min_reads} \
            --min-breadth {params.min_breadth} \
            > {log.logO} 2> {log.logE}
        """
