rule filter_representatives_minlen_500:
    input:
        reps_fasta = RESULTS_DIR + "/merged/viral_representatives_noApr24_ANI90_AF85.fasta",
        reps_tsv   = RESULTS_DIR + "/merged/viral_representatives_list_noApr24_ANI90_AF85.tsv"
    output:
        reps_fasta = RESULTS_DIR + "/merged/viral_representatives_noApr24_ANI90_AF85_min_500nt.fasta",
        reps_tsv   = RESULTS_DIR + "/merged/viral_representatives_list_noApr24_ANI90_AF85_min_500nt.tsv",
        stats      = RESULTS_DIR + "/merged/viral_representatives_noApr24_ANI90_AF85_min_500nt.stats.tsv",
        check      = RESULTS_DIR + "/merged/viral_representatives_noApr24_ANI90_AF85_min_500nt.done"
    log:
        logO = "logs/filter_representatives_minlen_500/merged.log",
        logE = "logs/filter_representatives_minlen_500/merged.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/filter_representatives_minlen.py \
          --reps-fasta {input.reps_fasta} \
          --reps-tsv {input.reps_tsv} \
          --out-fasta {output.reps_fasta} \
          --out-tsv {output.reps_tsv} \
          --stats {output.stats} \
          --min-len 500 \
          > {log.logO} 2> {log.logE}

        touch {output.check}
        """
