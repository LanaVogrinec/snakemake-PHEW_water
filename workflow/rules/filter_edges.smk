rule filter_edges:
    input:
        fasta = RESULTS_DIR + "/merged/viral_contigs.fasta",
        blast = RESULTS_DIR + "/merged/blastn_pairs_sorted_no-self_no_april2024.tsv.gz"
    output:
        pairs = RESULTS_DIR + "/merged/filtered_pairs_noApr24_ANI90_AF85.tsv",
        edges = RESULTS_DIR + "/merged/edges_noApr24_ANI90_AF85.tsv"
    params:
        min_ani = 90.0,
        filter_mode = "AF",
        min_hsp_len = 150,
        min_af = 0.85
    log:
        logO = "logs/filter_edges/merged.log",
        logE = "logs/filter_edges/merged.err.log"
    conda:
        "../envs/clustering_env.yaml"
    threads: 1
    shell:
        """
        python workflow/scripts/calculating_ANI_AF_longest_HSP_and_filtering.py \
            --blast_sorted_gz {input.blast} \
            --fasta {input.fasta} \
            --out_pairs {output.pairs} \
            --out_edges {output.edges} \
            --min_ani {params.min_ani} \
            --filter_mode {params.filter_mode} \
            --min_hsp_len {params.min_hsp_len} \
            --min_af {params.min_af} \
            > {log.logO} 2> {log.logE}
        """
