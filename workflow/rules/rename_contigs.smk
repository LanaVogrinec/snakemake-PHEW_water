rule rename_contigs:
    input:
        fasta = RESULTS_DIR + "/{sample_ID}/03_{sample_ID}_spades_contigs.fasta"
    output:
        renamed = RESULTS_DIR + "/{sample_ID}/03_{sample_ID}_spades_contigs_renamed.fasta"
    log:
        logO = "logs/rename_contigs/{sample_ID}.out.log",
        logE = "logs/rename_contigs/{sample_ID}.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/rename_contigs.py \
            --input {input.fasta} \
            --output {output.renamed} \
            --prefix {wildcards.sample_ID} \
            > {log.logO} 2> {log.logE}
        """
