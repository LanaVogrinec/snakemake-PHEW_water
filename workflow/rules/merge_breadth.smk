rule merge_breadth:
    input:
        breadth_files = expand(
            RESULTS_DIR + "/{sample_ID}/09_{sample_ID}_breadth.tsv",
            sample_ID = samples["samples"])
    output:
        breadth_matrix = RESULTS_DIR + "/merged/merged_breadth.tsv"
    log:
        logO = "logs/merge_breadth/merge_breadth.log",
        logE = "logs/merge_breadth/merge_breadth.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/combine_breadth.py \
            --breadth-files {input.breadth_files} \
            --out {output.breadth_matrix} \
            > {log.logO} 2> {log.logE}
        """
