rule isoform_count:
    input:
        isoform_quants=expand("results/step1/rsem_results/{Sample}/{Sample}.isoforms.results", Sample=Samples),
    output:
        "results/step1/isoform_count_matrix.txt"
    log:
        "logs/isoform_count.log"
    script:
        "../scripts/rsem_isoform_countMatrix.py"
