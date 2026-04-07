rule select_reps:
    input:
        fasta    = RESULTS_DIR + "/merged/viral_contigs.fasta",
        clusters = RESULTS_DIR + "/merged/test_Apr24_DS1/clusters_non-hybrid_no_DS1.tsv"
    output:
        reps_fasta = RESULTS_DIR + "/merged/test_Apr24_DS1/cluster_representatives_no_DS1.fasta",
        reps_tsv   = RESULTS_DIR + "/merged/test_Apr24_DS1/cluster_representatives_list_no_DS1.tsv",
        stats      = RESULTS_DIR + "/merged/test_Apr24_DS1/cluster_representatives_stats_no_DS1.tsv"
    log:
        logO = "logs/select_reps/test_Apr24_DS1/select_reps.log",
        logE = "logs/select_reps/test_Apr24_DS1/select_reps.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/select_reps.py \
          --fasta {input.fasta} \
          --clusters {input.clusters} \
          --reps-fasta {output.reps_fasta} \
          --reps-tsv {output.reps_tsv} \
          --stats {output.stats} \
          > {log.logO} 2> {log.logE}
        """
