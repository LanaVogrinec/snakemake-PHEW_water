rule merge_read_counts_controls:
    input:
        reps = RESULTS_DIR + "/merged/cluster_representatives_list.tsv",
        nkis = expand(
            RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_filtered_reps.tsv",
            sample_ID=nki_ids
        ),
        carry = expand(
            RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_filtered_reps.tsv",
            sample_ID=carry_ids
        )
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
            --mode controls \
            --nki-files {input.nkis} \
            --carry-files {input.carry} \
            --out-nki {output.nki} \
            --out-carry {output.carryover} \
            > {log.logO} 2> {log.logE}
        """