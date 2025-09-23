rule mmseqs2:
    input:
        f = RESULTS_DIR + "/{sample_ID}/03_{sample_ID}_spades_contigs_renamed.fasta"
    output:
        querydb_idx  = RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_queryDB.index",
        resultdb_idx = RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_resultDB.index",
        tmpdir       = temp(directory(RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_tmp"))
    params:
        queryprefix  = RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_queryDB",
        resultprefix = RESULTS_DIR + "/{sample_ID}/04_{sample_ID}_resultDB"
    log:
        logO = "logs/mmseqs2/{sample_ID}.db.log",
        logE = "logs/mmseqs2/{sample_ID}.db.err.log"
    conda:
        "../envs/mmseqs2_env.yaml"
    threads: 30
    shell:
        r"""
        mmseqs createdb {input.f} {params.queryprefix} \
          > {log.logO} 2> {log.logE}

        mmseqs taxonomy \
          {params.queryprefix} \
          /biodbs/mmseqs2/nr_database/nr.fnaDB \
          {params.resultprefix} \
          {output.tmpdir} \
          --threads {threads} --tax-lineage 1 \
          >> {log.logO} 2>> {log.logE}
        """
