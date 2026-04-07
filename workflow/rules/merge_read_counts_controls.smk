rule merge_read_counts_controls:
    input:
        reps = RESULTS_DIR + "/merged/test_Apr24_DS1/cluster_representatives_list_no_DS1.tsv",
        taxonomy = RESULTS_DIR + "/merged/merged_taxonomy.tsv"
    output:
        nki = RESULTS_DIR + "/merged/test_Apr24_DS1/reps_read_counts_nki_no_DS1.tsv",
        carryover = RESULTS_DIR + "/merged/test_Apr24_DS1/reps_read_counts_carry-over_no_DS1.tsv"
    log:
        logO = "logs/merge_read_counts/test_Apr24_DS1/merge_read_counts_controls.log",
        logE = "logs/merge_read_counts/test_Apr24_DS1/merge_read_counts_controls.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/merge_read_counts.py \
            --viral-reps {input.reps} \
            --taxonomy {input.taxonomy} \
            --mapping-root {RESULTS_DIR} \
            --out-NKIs {output.nki} \
            --out-carryovers {output.carryover} \
            > {log.logO} 2> {log.logE}
        """