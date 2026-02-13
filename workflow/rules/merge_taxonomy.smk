rule merge_taxonomy:
    input:
        taxonomy = expand(
            RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_mmseqs2_taxonomy.tsv",
            sample_ID=samples["samples"])
    output:
        merged = RESULTS_DIR + "/merged/merged_taxonomy.tsv"
    log:
        logO = "logs/merge_taxonomy/merge_taxonomy.out.log",
        logE = "logs/merge_taxonomy/merged_taxonomy.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/merged
        cat {input.taxonomy} > {output.merged} 2> {log.logE} > {log.logO}
        """
