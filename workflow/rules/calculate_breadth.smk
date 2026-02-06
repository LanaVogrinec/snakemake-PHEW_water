rule calculate_breadth:
    input:
        bam = RESULTS_DIR + "/{sample_ID}/08_{sample_ID}_coverm_bbmap_rep_all_minlen_500_noApr24_ANI90_AF85.filtered.bam"
    output:
        breadth = RESULTS_DIR + "/{sample_ID}/09_{sample_ID}_breadth_rep_all_minlen_500_noApr24_ANI90_AF85.tsv",
        check   = RESULTS_DIR + "/{sample_ID}/09_{sample_ID}_breadth_rep_all_minlen_500_noApr24_ANI90_AF85.done"
    log:
        logO = "logs/calculate_breadth_minlen_500/{sample_ID}.log",
        logE = "logs/calculate_breadth_minlen_500/{sample_ID}.err.log"
    conda:
        "../envs/calculate_breadth_env.yaml"
    threads: 20
    shell:
        """
        mkdir -p $(dirname {output.breadth})

        cov_file=$(mktemp $(dirname {output.breadth})/tmp_genomecov_XXXXXX.txt)

        bedtools genomecov -d -ibam {input.bam} > "$cov_file" 2>> {log.logE}

        python workflow/scripts/summarize_breadth.py "$cov_file" {output.breadth} >> {log.logO} 2>> {log.logE}

        rm -f "$cov_file"

        touch {output.check}
        """
