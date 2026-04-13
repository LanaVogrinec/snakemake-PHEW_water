rule merge_read_counts_samples:
    input:
        reps = RESULTS_DIR + "/merged/cluster_representatives_list.tsv",
        taxonomy = RESULTS_DIR + "/merged/merged_taxonomy.tsv",
        coverm = expand(
            RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_filtered_reps.tsv",
            sample_ID=sample_ids
        )
    output:
        reps_tax = RESULTS_DIR + "/merged/cluster_representatives_taxonomy.tsv",
        samples = RESULTS_DIR + "/merged/reps_read_counts_samples.tsv"
    log:
        logO = "logs/merge_read_counts/merge_read_counts_samples.log",
        logE = "logs/merge_read_counts/merge_read_counts_samples.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/merge_read_counts.py \
            --viral-reps {input.reps} \
            --taxonomy {input.taxonomy} \
            --mode samples \
            --sample-list {input.coverm} \
            --out {output.samples} \
            > {log.logO} 2> {log.logE}
        """