rule bowtie2_index:
    input:
        merged = RESULTS_DIR + "/merged/viral_contigs.fasta"
    output:
        idx = expand(RESULTS_DIR + "/merged/bt2_viral_contigs_index/viral_contigs.{ext}",
                     ext=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"])
    log:
        logO = "logs/bowtie2_index/viral_contigs.log",
        logE = "logs/bowtie2_index/viral_contigs.err.log"
    conda:
        "../envs/mapping_env.yaml"
    threads: 32
    shell:
        r"""
        mkdir -p {RESULTS_DIR}/merged/bt2_viral_contigs_index
        bowtie2-build --threads {threads} {input.merged} {RESULTS_DIR}/merged/bt2_viral_contigs_index/viral_contigs \
            > {log.logO} 2> {log.logE}
        """
