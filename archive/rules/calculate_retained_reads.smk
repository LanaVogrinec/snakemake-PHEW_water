rule calculate_retained_reads:
    input:
        cutadapt_log = "logs/cutadapt/{sample_ID}.log",
        original_bam = RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_rep_all_sorted.bam",
        filtered_bam = RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_bbmap_rep_all.filtered.bam"
    output:
        report = RESULTS_DIR + "/{sample_ID}/10_{sample_ID}_retained_reads.txt",
        check  = RESULTS_DIR + "/{sample_ID}/10_{sample_ID}_retained_reads.done"
    log:
        logO = "logs/retention/{sample_ID}.log",
        logE = "logs/retention/{sample_ID}.err.log"
    conda:
        "../envs/mapping_env.yaml"
    shell:
        """
        bash workflow/scripts/calculating_retained_reads.sh \
            {input.cutadapt_log} \
            {input.original_bam} \
            {input.filtered_bam} \
            {output.report} \
            > {log.logO} 2> {log.logE}

        touch {output.check}
        """
