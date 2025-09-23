rule spades:
    input:
        r1 = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_1_paired.fq.gz",
        r2 = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_2_paired.fq.gz"
    output:
        d = temp(directory(RESULTS_DIR + "/{sample_ID}/{sample_ID}_spades")),
        f = RESULTS_DIR + "/{sample_ID}/03_{sample_ID}_spades_contigs.fasta"
    log:
        logO = "logs/spades/{sample_ID}.log",
        logE = "logs/spades/{sample_ID}.err.log"
    conda:
        "../envs/spades_env.yaml"
    threads: 30
    shell:
        """
        metaspades.py -t {threads} \
            -1 {input.r1} -2 {input.r2} \
            -o {output.d} \
            > {log.logO} 2> {log.logE}

        cp {output.d}/contigs.fasta {output.f}
        """
