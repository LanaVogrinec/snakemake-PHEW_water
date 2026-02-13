rule coverm_quant_bbmap_luciferase:
    input:
        bam = RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_luciferase.bam",
        bai = RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_luciferase.bam.bai"
    output:
        tsv     = RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_bbmap_luciferase.tsv",
        filtered_bam = RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_bbmap_luciferase.bam",
        check   = RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_bbmap_luciferase.done"
    log:
        logO = "logs/coverm_bbmap_luciferase/{sample_ID}.log",
        logE = "logs/coverm_bbmap_luciferase/{sample_ID}.err.log"
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
