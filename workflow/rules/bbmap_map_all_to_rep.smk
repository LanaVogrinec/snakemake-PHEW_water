rule bbmap_map_all_to_rep:
    input:
        r1  = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_1_paired.fq.gz",
        r2  = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_2_paired.fq.gz",
        ref = RESULTS_DIR + "/merged/viral_representatives_noApr24_ANI90_AF85_min_500nt.fasta"
    output:
        bam   = temp(RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_rep_all_sorted_minlen_500_noApr24_ANI90_AF85.bam"),
        bai   = temp(RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_rep_all_sorted_minlen_500_noApr24_ANI90_AF85.bam.bai"),
        check = RESULTS_DIR + "/{sample_ID}/07_{sample_ID}_bbmap_rep_all_mapping_minlen_500_noApr24_ANI90_AF85.done"
    log:
        logO = "logs/bbmap_rep_all_minlen_500_noApr24_ANI90_AF85/{sample_ID}.log",
        logE = "logs/bbmap_rep_all_minlen_500_noApr24_ANI90_AF85/{sample_ID}.err.log"
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
