rule merge_contigs:
    input:
        expand(
            RESULTS_DIR + "/{sample_ID}/03_{sample_ID}_spades_contigs_renamed.fasta",
            sample_ID=samples["samples"])
    output:
        merged = RESULTS_DIR + "/merged/merged_contigs.fasta"
    log:
        logO = "logs/merge_contigs/merge_contigs.out.log",
        logE = "logs/merge_contigs/merge_contigs.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        cat {input} > {output.merged} 2> {log.logE} > {log.logO}
        """
