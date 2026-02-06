rule coverm_quant_bt2:
    input:
        bam = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bt2_sorted.bam",
        bai = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bt2_sorted.bam.bai"
    output:
        tsv   = RESULTS_DIR + "/{sample_ID}/06_{sample_ID}_coverm_bt2.tsv",
        check = RESULTS_DIR + "/{sample_ID}/06_{sample_ID}_coverm_bt2.done"
    log:
        logO = "logs/coverm/{sample_ID}_bt2.log",
        logE = "logs/coverm/{sample_ID}_bt2.err.log"
    conda:
        "../envs/mapping_env.yaml"
    threads: 12
    shell:
        r"""
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
