rule extract_viral_contigs:
    input:
        fasta = RESULTS_DIR + "/merged/merged_contigs.fasta",
        taxonomy = RESULTS_DIR + "/merged/merged_taxonomy.tsv"
    output:
        viral_fasta = RESULTS_DIR + "/merged/viral_contigs.fasta"
    log:
        logO = "logs/extract_viral_contigs/extract_viral_contigs.log",
        logE = "logs/extract_viral_contigs/extract_viral_contigs.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/extract_viral_contigs.py \
            --fasta {input.fasta} \
            --taxonomy {input.taxonomy} \
            --output {output.viral_fasta} \
            > {log.logO} 2> {log.logE}
        """
