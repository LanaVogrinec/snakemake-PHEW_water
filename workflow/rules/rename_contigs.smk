rule rename_contigs:
    input:
        fasta = RESULTS_DIR + "/{sample_ID}/03_{sample_ID}_spades_contigs.fasta"
    output:
        renamed = RESULTS_DIR + "/{sample_ID}/03_{sample_ID}_spades_contigs_renamed.fasta"
    log:
        log = "logs/rename_contigs/{sample_ID}.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/renaming_contigs_script.py \
            --input {input.fasta} \
            --output {output.renamed} \
            --prefix {wildcards.sample_ID} \
            > {log.log} 2>&1
        """
