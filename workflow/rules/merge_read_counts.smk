rule merge_read_counts:
    input:
        reps = RESULTS_DIR + "/merged/cluster_representatives_list.tsv",
        taxonomy = RESULTS_DIR + "/merged/merged_taxonomy.tsv",
        mapping_files = expand(
            RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_filtered_reps.tsv",
            sample_ID=samples["samples"]
        )
    output:
        reps_tax = RESULTS_DIR + "/merged/cluster_representatives_taxonomy.tsv",
        samples = RESULTS_DIR + "/merged/reps_read_counts_samples.tsv",
        nki = RESULTS_DIR + "/merged/reps_read_counts_nki.tsv",
        carryover = RESULTS_DIR + "/merged/reps_read_counts_carry-over.tsv"
    log:
        logO = "logs/merge_read_counts/merge_read_counts.log",
        logE = "logs/merge_read_counts/merge_read_counts.err.log"
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
            --out-nki {output.nki} \
            --out-carryover {output.carryover} \
            > {log.logO} 2> {log.logE}
        """
