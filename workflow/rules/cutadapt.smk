rule cutadapt:
    input:
        r1 = RESULTS_DIR + "/{sample_ID}/01_{sample_ID}_trim_truseq_1_paired.fq.gz",
        r2 = RESULTS_DIR + "/{sample_ID}/01_{sample_ID}_trim_truseq_2_paired.fq.gz"
    output:
        tp1 = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_1_paired.fq.gz",
        tp2 = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_2_paired.fq.gz",
        check = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera.done"
    log:
        logO = "logs/cutadapt/{sample_ID}.log",
        logE = "logs/cutadapt/{sample_ID}.err.log"
    conda:
        "../envs/cutadapt_env.yaml"
    threads: 4
    shell:
        """
        cutadapt -j {threads} \
          -g GTTTCCCAGTCACGATA -G GTTTCCCAGTCACGATA \
          -a TATCGTGACTGGGAAAC -A TATCGTGACTGGGAAAC \
          -o {output.tp1} -p {output.tp2} \
          {input.r1} {input.r2} \
          > {log.logO} 2> {log.logE}
        touch {output.check}
        """
