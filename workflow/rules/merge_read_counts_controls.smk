rule merge_read_counts_controls:
    input:
        reps = RESULTS_DIR + "/merged/cluster_representatives_list.tsv",
        taxonomy = RESULTS_DIR + "/merged/merged_taxonomy.tsv"
    output:
        nki = RESULTS_DIR + "/merged/reps_read_counts_nki.tsv",
        carryover = RESULTS_DIR + "/merged/reps_read_counts_carry-over.tsv"
    log:
        logO = "logs/merge_read_counts/merge_read_counts_controls.log",
        logE = "logs/merge_read_counts/merge_read_counts_controls.err.log"
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