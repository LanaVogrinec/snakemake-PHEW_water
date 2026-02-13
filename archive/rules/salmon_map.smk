rule salmon_map:
    input:
        index = RESULTS_DIR + "/merged/salmon_index",
        r1 = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_1_paired.fq.gz",
        r2 = RESULTS_DIR + "/{sample_ID}/02_{sample_ID}_trim_primera_2_paired.fq.gz"
    output:
        quant_dir = RESULTS_DIR + "/{sample_ID}/salmon_map"
    log:
        logO = "logs/salmon/{sample_ID}.log",
        logE = "logs/salmon/{sample_ID}.err.log"
    conda:
        "../envs/salmon_env.yaml"
    threads: 40
    shell:
        """
        salmon quant -i {input.index} \
                     -l A \
                     -1 {input.r1} -2 {input.r2} \
                     -p {threads} \
                     --validateMappings \
                     -o {output.quant_dir} \
                     > {log.logO} 2> {log.logE}
        """
