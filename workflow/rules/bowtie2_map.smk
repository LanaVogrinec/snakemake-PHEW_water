rule bowtie2_map:
    input:
        r1  = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_1_paired.fq.gz",
        r2  = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_2_paired.fq.gz",
        idx = rules.bowtie2_index.output.idx
    output:
        bam   = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bt2_sorted.bam",
        bai   = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bt2_sorted.bam.bai",
        check = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bt2_mapping.done"
    log:
        logO = "logs/bowtie2_map/{sample_ID}.log",
        logE = "logs/bowtie2_map/{sample_ID}.err.log"
    conda:
        "../envs/mapping_env.yaml"
    threads: 13
    shell:
        r"""
        bowtie2 \
          --sensitive-local \
          -x {RESULTS_DIR}/merged/bt2_index/merged_contigs \
          -1 {input.r1} -2 {input.r2} \
          -p {threads} \
          --no-unal \
          2> {log.logE} \
        | samtools view -b -F 4 -@ {threads} - \
        | samtools sort -@ {threads} -o {output.bam} -

        samtools index -@ {threads} {output.bam} > {log.logO} 2>> {log.logE}
        touch {output.check}
        """
