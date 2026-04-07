rule filter_edges:
    input:
        fasta = RESULTS_DIR + "/merged/viral_contigs.fasta",
        blast = RESULTS_DIR + "/merged/test_Apr24_DS1/blastn_pairs_sorted_no_self_no_DS1.tsv.gz"
    output:
        pairs = RESULTS_DIR + "/merged/test_Apr24_DS1/blastn_pairs_sorted_no_self_filtered_no_DS1.tsv",
        edges = RESULTS_DIR + "/merged/test_Apr24_DS1/edges_filtered_no_DS1.tsv"
    params:
        min_ani = 90.0,
        filter_mode = "AF",
        min_hsp_len = 150,
        min_af = 0.85
    log:
        logO = "logs/filter_edges/test_Apr24_DS1/filter_edges.log",
        logE = "logs/filter_edges/test_Apr24_DS1/filter_edges.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/filter_edges.py \
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
