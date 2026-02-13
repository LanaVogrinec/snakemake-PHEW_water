rule trimmomatic:
    input:
        r1 = "resources/{sample_ID}/{sample_ID}_1.fq.gz",
        r2 = "resources/{sample_ID}/{sample_ID}_2.fq.gz",
        iclip = "resources/TruSeq3-PE.fa"
    output:
        tp1 = RESULTS_DIR + "/{sample_ID}/01_{sample_ID}_trim_truseq_1_paired.fq.gz",
        tp2 = RESULTS_DIR + "/{sample_ID}/01_{sample_ID}_trim_truseq_2_paired.fq.gz",
        tup1 = RESULTS_DIR + "/{sample_ID}/01_{sample_ID}_trim_truseq_1_unpaired.fq.gz",
        tup2 = RESULTS_DIR + "/{sample_ID}/01_{sample_ID}_trim_truseq_2_unpaired.fq.gz",
        check = RESULTS_DIR + "/{sample_ID}/01_{sample_ID}_trim_truseq.done"
    log:
        logO = "logs/trimmomatic/{sample_ID}.log",
        logE = "logs/trimmomatic/{sample_ID}.err.log"
    conda:
        "../envs/trimmomatic_env.yaml"
    threads: 4
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 \
            {input.r1} {input.r2} \
            {output.tp1} {output.tup1} {output.tp2} {output.tup2} \
            ILLUMINACLIP:{input.iclip}:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 \
            > {log.logO} 2> {log.logE}

        touch {output.check}
        """
