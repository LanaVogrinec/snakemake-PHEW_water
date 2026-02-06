rule bbmap_map_rep:
    input:
        r1  = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bbmap_viral_R1.fastq.gz",
        r2  = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bbmap_viral_R2.fastq.gz",
        ref = RESULTS_DIR + "/merged/viral_representatives.fasta"
    output:
        bam   = RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_rep_sorted.bam",
        bai   = RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_rep_sorted.bam.bai",
        check = RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_rep_mapping.done"
    log:
        logO = "logs/bbmap_rep/{sample_ID}.log",
        logE = "logs/bbmap_rep/{sample_ID}.err.log"
    conda:
        "../envs/mapping_env.yaml"
    threads: 20
    shell:
        """
        bbmap.sh \
            in={input.r1} in2={input.r2} \
            ref={input.ref} \
            ambig=best \
            minid=0.90 \
            maxindel=3 \
            threads={threads} \
            vslow=t \
            outm=stdout \
            2> {log.logE} \
        | samtools view -b -@ {threads} - 2>> {log.logE} \
        | samtools sort -@ {threads} -o {output.bam} \
            > {log.logO} 2>> {log.logE}

        samtools index -@ {threads} {output.bam} 2>> {log.logE}

        touch {output.check}
        """
