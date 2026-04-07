rule coverm:
    input:
        bam = RESULTS_DIR + "/{sample_ID}/test_Apr24_DS1/07_{sample_ID}_bbmap_reps_no_DS1.bam",
        bai = RESULTS_DIR + "/{sample_ID}/test_Apr24_DS1/07_{sample_ID}_bbmap_reps_no_DS1.bam.bai"
    output:
        tsv     = RESULTS_DIR + "/{sample_ID}/test_Apr24_DS1/08_{sample_ID}_coverm_filtered_reps_no_DS1.tsv",
        filtered_bam = RESULTS_DIR + "/{sample_ID}/test_Apr24_DS1/08_{sample_ID}_coverm_filtered_reps_no_DS1.bam",
        check   = RESULTS_DIR + "/{sample_ID}/test_Apr24_DS1/08_{sample_ID}_coverm_filtered_reps_no_DS1.done"
    log:
        logO = "logs/coverm/test_Apr24_DS1/{sample_ID}.log",
        logE = "logs/coverm/test_Apr24_DS1/{sample_ID}.err.log"
    conda:
        "../envs/mapping_env.yaml"
    threads: 12
    shell:
        """
        coverm filter \
            --bam-files {input.bam} \
            --threads {threads} \
            --min-read-percent-identity 95 \
            --min-read-aligned-percent 95 \
            --output-bam-files {output.filtered_bam} \
            > {log.logO} 2> {log.logE}

        samtools index {output.filtered_bam} 2>> {log.logE}

        coverm contig \
            --bam-files {output.filtered_bam} \
            --threads {threads} \
            --methods count tpm rpkm \
            --output-file {output.tsv} \
            >> {log.logO} 2>> {log.logE}

        touch {output.check}
        """
