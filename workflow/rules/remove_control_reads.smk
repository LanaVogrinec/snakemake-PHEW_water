rule remove_control_reads:
    input:
        samples = RESULTS_DIR + "/merged/reps_read_counts_samples.tsv",
        carryover = RESULTS_DIR + "/merged/reps_read_counts_carry-over.tsv",
        mapping = RESULTS_DIR + "/merged/carry-over_matches.tsv"
    output:
        cleaned = RESULTS_DIR + "/merged/reps_read_counts_samples_no_carryover.tsv"
    log:
        logO = "logs/filter_read_counts/filter_carryover.log",
        logE = "logs/filter_read_counts/filter_carryover.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/remove_control_reads.py \
            --samples {input.samples} \
            --carryovers {input.carryover} \
            --mapping {input.mapping} \
            --out-filtered {output.cleaned} \
            > {log.logO} 2> {log.logE}
        """