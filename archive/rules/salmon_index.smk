rule salmon_index:
    input:
        merged = RESULTS_DIR + "/merged/merged_contigs.fasta"
    output:
        index_dir = RESULTS_DIR + "/merged/salmon_index"
    log:
        "logs/salmon/index.log"
    conda:
        "../envs/salmon_env.yaml"
    threads: 20
    shell:
        """
        salmon index -t {input.merged} -i {output.index_dir} -p {threads} > {log} 2>&1
        """
