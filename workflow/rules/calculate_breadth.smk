rule calculate_breadth:
    input:
        bam = RESULTS_DIR + "/{sample_ID}/test_Apr24_DS1/08_{sample_ID}_coverm_filtered_reps_no_DS1.bam"
    output:
        breadth = RESULTS_DIR + "/{sample_ID}/test_Apr24_DS1/09_{sample_ID}_breadth_no_DS1.tsv"
    log:
        logO = "logs/calculate_breadth/test_Apr24_DS1/{sample_ID}.log",
        logE = "logs/calculate_breadth/test_Apr24_DS1/{sample_ID}.err.log"
    conda:
        "../envs/calculate_breadth_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.breadth})

        cov_file=$(mktemp $(dirname {output.breadth})/tmp_genomecov_XXXXXX.txt)

        bedtools genomecov -d -ibam {input.bam} > "$cov_file" 2>> {log.logE}

        python workflow/scripts/calculate_breadth.py "$cov_file" {output.breadth} >> {log.logO} 2>> {log.logE}

        rm -f "$cov_file"
        """  
