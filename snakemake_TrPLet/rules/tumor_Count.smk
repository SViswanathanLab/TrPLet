rule gene_count:
    input:
        isoform_quants=expand("results/step1/rsem_results/{Sample}/{Sample}.genes.results", Sample=Samples),
    output:
        "results/step1/gene_count_matrix.txt"
    log:
        "logs/gene_count.log"
    script:
        "../scripts/rsem_gene_countMatrix.py"
