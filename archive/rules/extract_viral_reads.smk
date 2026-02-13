rule extract_viral_reads:
    input:
        bam = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bbmap_sorted.bam",
        viral_fasta = RESULTS_DIR + "/merged/viral_contigs.fasta"
    output:
        r1 = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bbmap_viral_R1.fastq.gz",
        r2 = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bbmap_viral_R2.fastq.gz",
        check = RESULTS_DIR + "/{sample_ID}/05_{sample_ID}_bbmap_viral_reads.done"
    conda:
        "../envs/mapping_env.yaml"
    threads: 12
    shell:
        """
        python workflow/scripts/extracting_viral_reads.py \
            --bam {input.bam} \
            --viral-fasta {input.viral_fasta} \
            --out-r1 {output.r1} \
            --out-r2 {output.r2}

        touch {output.check}
        """
