rule bbmap_to_reps:
    input:
        r1  = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_1_paired.fq.gz",
        r2  = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_2_paired.fq.gz",
        ref = RESULTS_DIR + "/merged/test_Apr24_DS1/cluster_representatives_filtered_length_no_DS1.fasta"
    output:
        bam   = temp(RESULTS_DIR + "/{sample_ID}/test_Apr24_DS1/07_{sample_ID}_bbmap_reps_no_DS1.bam"),
        bai   = temp(RESULTS_DIR + "/{sample_ID}/test_Apr24_DS1/07_{sample_ID}_bbmap_reps_no_DS1.bam.bai"),
        check = RESULTS_DIR + "/{sample_ID}/test_Apr24_DS1/07_{sample_ID}_bbmap_reps_no_DS1.done"
    log:
        logE = "logs/bbmap_to_reps/test_Apr24_DS1/{sample_ID}.err.log"
    conda:
        "../envs/mapping_env.yaml"
    threads: 20
    shell:
        r"""
        bbmap.sh \
            nodisk=t \
            in={input.r1} in2={input.r2} \
            ref={input.ref} \
            ambig=random \
            minid=0.90 \
            maxindel=3 \
            threads={threads} \
            out={output.bam}.sam \
            2> {log.logE}

        samtools view -b -@ {threads} {output.bam}.sam \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index -@ {threads} {output.bam} 2>> {log.logE}

        rm {output.bam}.sam
        touch {output.check}
        """
