rule coverm_quant_bbmap_all_to_rep:
    input:
        bam = RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_rep_all_sorted_minlen_500_noApr24_ANI90_AF85.bam",
        bai = RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_rep_all_sorted_minlen_500_noApr24_ANI90_AF85.bam.bai"
    output:
        tsv     = RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_bbmap_rep_all_minlen_500_noApr24_ANI90_AF85.tsv",
        filtered_bam = RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_bbmap_rep_all_minlen_500_noApr24_ANI90_AF85.filtered.bam",
        check   = RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_bbmap_rep_all_minlen_500_noApr24_ANI90_AF85.done"
    log:
        logO = "logs/coverm_bbmap_rep_all_minlen_500_noApr24_ANI90_AF85/{sample_ID}.log",
        logE = "logs/coverm_bbmap_rep_all_minlen_500_noApr24_ANI90_AF85/{sample_ID}.err.log"
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
