rule mmseqs2_tsv:
    input:
        querydb_idx  = RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_queryDB.index",
        resultdb_idx = RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_resultDB.index"
    params:
        qprefix = RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_queryDB",
        rprefix = RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_resultDB"
    output:
        tsv = RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_mmseqs2_taxonomy.tsv"
    log:
        logO = "logs/mmseqs2_tsv/{sample_ID}.tsv.log",
        logE = "logs/mmseqs2_tsv/{sample_ID}.tsv.err.log"
    conda:
        "../envs/mmseqs2_env.yaml"
    threads: 4
    shell:
        r"""
        mmseqs createtsv \
          {params.qprefix} \
          {params.rprefix} \
          {output.tsv} \
          > {log.logO} 2> {log.logE}
        """
