rule combine_read_counts_taxonomy:
    input:
        reps = RESULTS_DIR + "/merged/viral_representatives_list_noApr24_ANI90_AF85_min_500nt.tsv",
        taxonomy = RESULTS_DIR + "/merged/merged_taxonomy.tsv",
        mapping_files = expand(
            RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_bbmap_rep_all_minlen_500_noApr24_ANI90_AF85.tsv",
            sample_ID = samples["samples"]
        )
    output:
        reps_tax = RESULTS_DIR + "/merged/viral_representatives_minlen500_taxonomy_noApr24_ANI90_AF85.tsv",
        samples = RESULTS_DIR + "/merged/rep_contigs_minlen500_read_counts_taxonomy_samples_noApr24_ANI90_AF85.tsv",
        NKIs = RESULTS_DIR + "/merged/rep_contigs_minlen500_read_counts_taxonomy_NKIs_noApr24_ANI90_AF85.tsv",
        carryovers = RESULTS_DIR + "/merged/rep_contigs_minlen500_read_counts_taxonomy_carryovers_noApr24_ANI90_AF85.tsv",
        check   = RESULTS_DIR + "/merged/rep_contigs_minlen500_taxonomy_read_counts_noApr24_ANI90_AF85.done"
    log:
        logO = "logs/combine_taxonomy_and_read_counts/merged.log",
        logE = "logs/combine_taxonomy_and_read_counts/merged.err.log"
    conda:
        "../envs/calculate_breadth_env.yaml"
    shell:
        """
        python workflow/scripts/combine_read_counts_taxonomy.py \
            --viral-reps {input.reps} \
            --taxonomy {input.taxonomy} \
            --mapping-root {RESULTS_DIR} \
            --out-reps-tax {output.reps_tax} \
            --out-samples {output.samples} \
            --out-NKIs {output.NKIs} \
            --out-carryovers {output.carryovers} \
            > {log.logO} 2> {log.logE}

        touch {output.check}
        """
