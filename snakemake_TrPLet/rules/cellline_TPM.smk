rule isoform_count:
    input:
        isoform_quants=expand("results/step1/rsem_results/{Sample}/{Sample}.genes.results", Sample=Samples),
    output:
        "results/step1/gene_log2_TPM+1_matrix.txt"
    log:
        "logs/isoform_count.log"
    script:
        "../scripts/rsem_gene_TPM_Matrix.py"
