rule bbmap_map:
    input:
        r1  = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_1_paired.fq.gz",
        r2  = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_2_paired.fq.gz",
        ref = RESULTS_DIR + "/merged/merged_contigs.fasta"
    output:
        bam   = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bbmap_sorted.bam",
        bai   = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bbmap_sorted.bam.bai",
        check = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bbmap_mapping.done"
    log:
        logO = "logs/bbmap/{sample_ID}.log",
        logE = "logs/bbmap/{sample_ID}.err.log"
    conda:
        "../envs/mapping_env.yaml"
    threads: 20
    shell:
        """
        bbmap.sh \
          in={input.r1} in2={input.r2} \
          ref={input.ref} \
          ambig=random \
          minid=0.90 \
          maxindel=3 \
          threads={threads} \
          vslow=t \
          outm=stdout \
        | samtools view -b -@ {threads} - \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index -@ {threads} {output.bam} > {log.logO} 2> {log.logE}

        touch {output.check}
        """
