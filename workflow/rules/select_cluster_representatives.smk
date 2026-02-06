rule select_cluster_representatives:
    input:
        fasta    = RESULTS_DIR + "/merged/viral_contigs.fasta",
        clusters = RESULTS_DIR + "/merged/final_cluster_assignments_non-hybrid_noApr24_ANI90_AF85.tsv"
    output:
        reps_fasta = RESULTS_DIR + "/merged/viral_representatives_noApr24_ANI90_AF85.fasta",
        reps_tsv   = RESULTS_DIR + "/merged/viral_representatives_list_noApr24_ANI90_AF85.tsv",
        stats      = RESULTS_DIR + "/merged/viral_representatives_noApr24_ANI90_AF85.stats.tsv",
        check      = RESULTS_DIR + "/merged/viral_representatives_noApr24_ANI90_AF85.done"
    log:
        logO = "logs/select_cluster_representatives/merged.log",
        logE = "logs/select_cluster_representatives/merged.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/selecting_cluster_representatives.py \
          --fasta {input.fasta} \
          --clusters {input.clusters} \
          --reps-fasta {output.reps_fasta} \
          --reps-tsv {output.reps_tsv} \
          --stats {output.stats} \
          > {log.logO} 2> {log.logE}

        touch {output.check}
        """
