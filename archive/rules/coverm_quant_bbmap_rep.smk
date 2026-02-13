rule coverm_quant_bbmap_rep:
    input:
        bam = RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_rep_sorted.bam",
        bai = RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_rep_sorted.bam.bai"
    output:
        tsv   = RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_bbmap_rep.tsv",
        check = RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_bbmap_rep.done"
    log:
        logO = "logs/coverm_bbmap_rep/{sample_ID}.log",
        logE = "logs/coverm_bbmap_rep/{sample_ID}.err.log"
    conda:
        "../envs/mapping_env.yaml"
    threads: 12
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --threads {threads} \
            --min-read-percent-identity 95 \
            --min-read-aligned-percent 95 \
            --methods count tpm rpkm \
            --output-file {output.tsv} \
            > {log.logO} 2> {log.logE}

        touch {output.check}
        """
