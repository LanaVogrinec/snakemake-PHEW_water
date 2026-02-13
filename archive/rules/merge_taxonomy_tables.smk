rule merge_taxonomy_tables:
    input:
        taxonomy = expand(RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_mmseqs2_taxonomy.tsv",
                          sample_ID=samples["samples"])
    output:
        merged = RESULTS_DIR + "/merged/merged_taxonomy.tsv"
    shell:
        """
        mkdir -p {RESULTS_DIR}/merged
        cat {input.taxonomy} > {output.merged}
        """
