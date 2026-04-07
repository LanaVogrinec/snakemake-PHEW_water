rule calculate_rpkm:
    input:
        counts = RESULTS_DIR + "/merged/test_Apr24_DS1/reps_read_counts_samples_filtered_no_DS1.tsv"
    output:
        rpkm_long = RESULTS_DIR + "/merged/test_Apr24_DS1/rpkm_long_no_DS1.tsv",
        pa_long   = RESULTS_DIR + "/merged/test_Apr24_DS1/presence_absence_long_no_DS1.tsv"
    log:
        logO = "logs/calculate_rpkm/test_Apr24_DS1/calculate_rpkm.log",
        logE = "logs/calculate_rpkm/test_Apr24_DS1/calculate_rpkm.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/calculate_rpkm.py \
            --input {input.counts} \
            --output_rpkm {output.rpkm_long} \
            --output_pa {output.pa_long} \
            > {log.logO} 2> {log.logE}
        """
