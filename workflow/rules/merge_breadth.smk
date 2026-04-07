rule merge_breadth:
    input:
        breadth_files = expand(
            RESULTS_DIR + "/{sample_ID}/test_Apr24_DS1/09_{sample_ID}_breadth_no_DS1.tsv",
            sample_ID = samples["samples"])
    output:
        breadth_matrix = RESULTS_DIR + "/merged/test_Apr24_DS1/merged_breadth_no_DS1.tsv"
    log:
        logO = "logs/merge_breadth/test_Apr24_DS1/merge_breadth.log",
        logE = "logs/merge_breadth/test_Apr24_DS1/merge_breadth.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/merge_breadth.py \
            --breadth-files {input.breadth_files} \
            --out {output.breadth_matrix} \
            > {log.logO} 2> {log.logE}
        """
