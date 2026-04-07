rule merge_read_counts_samples:
    input:
        reps = RESULTS_DIR + "/merged/test_Apr24_DS1/cluster_representatives_list_no_DS1.tsv",
        taxonomy = RESULTS_DIR + "/merged/merged_taxonomy.tsv"
    output:
        reps_tax = RESULTS_DIR + "/merged/test_Apr24_DS1/cluster_representatives_taxonomy_no_DS1.tsv",
        samples = RESULTS_DIR + "/merged/test_Apr24_DS1/reps_read_counts_samples_no_DS1.tsv"
    log:
        logO = "logs/merge_read_counts/test_Apr24_DS1/merge_read_counts_samples.log",
        logE = "logs/merge_read_counts/test_Apr24_DS1/merge_read_counts_samples.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/merge_read_counts.py \
            --viral-reps {input.reps} \
            --taxonomy {input.taxonomy} \
            --mapping-root {RESULTS_DIR} \
            --out-reps-tax {output.reps_tax} \
            --out-samples {output.samples} \
            > {log.logO} 2> {log.logE}
        """