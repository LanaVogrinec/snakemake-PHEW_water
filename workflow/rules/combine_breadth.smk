rule combine_breadth:
    input:
        breadth_files = expand(
            RESULTS_DIR + "/{sample_ID}/09_{sample_ID}_breadth_rep_all_minlen_500_noApr24_ANI90_AF85.tsv",
            sample_ID = samples["samples"]
        )
    output:
        breadth_matrix = RESULTS_DIR + "/merged/rep_contigs_minlen500_breadth_noApr24_ANI90_AF85.tsv",
        check = RESULTS_DIR + "/merged/rep_contigs_minlen500_breadth_noApr24_ANI90_AF85.done"
    log:
        logO = "logs/combine_breadth.log",
        logE = "logs/combine_breadth.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/combine_breadth.py \
            --breadth-files {input.breadth_files} \
            --out {output.breadth_matrix} \
            > {log.logO} 2> {log.logE}

        touch {output.check}
        """
