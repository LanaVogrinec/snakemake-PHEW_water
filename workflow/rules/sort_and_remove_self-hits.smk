rule sort_and_remove_self_hits:
    input:
        blast = RESULTS_DIR + "/merged/blastn_pairs.tsv"
    output:
        sorted = RESULTS_DIR + "/merged/blastn_pairs_sorted_no-self.tsv.gz"
    log:
        logO = "logs/sort_blast/merged.log",
        logE = "logs/sort_blast/merged.err.log"
    conda:
        "../envs/clustering_env.yaml"
    threads: 1
    shell:
        """
        python workflow/scripts/sorting_and_removing_self-hits.py \
            --input {input.blast} \
            --output {output.sorted} \
            > {log.logO} 2> {log.logE}
        """
