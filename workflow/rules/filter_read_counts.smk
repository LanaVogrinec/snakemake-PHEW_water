rule filter_read_counts:
    input:
        samples = RESULTS_DIR + "/merged/test_Apr24_DS1/reps_read_counts_samples_no_carryover_no_DS1.tsv",
        breadth = RESULTS_DIR + "/merged/test_Apr24_DS1/merged_breadth_no_DS1.tsv"
    output:
        filtered = RESULTS_DIR + "/merged/test_Apr24_DS1/reps_read_counts_samples_filtered_no_DS1.tsv",
        stats    = RESULTS_DIR + "/merged/test_Apr24_DS1/reps_read_counts_samples_filtered_stats_no_DS1.tsv"
    params:
        min_reads = 2,
        min_breadth = 0.25
    log:
        logO = "logs/filter_read_counts/test_Apr24_DS1/filter_breadth.log",
        logE = "logs/filter_read_counts/test_Apr24_DS1/filter_breadth.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/filter_read_counts.py \
            --samples {input.samples} \
            --breadth {input.breadth} \
            --out-filtered {output.filtered} \
            --out-stats {output.stats} \
            --min-reads {params.min_reads} \
            --min-breadth {params.min_breadth} \
            > {log.logO} 2> {log.logE}
        """