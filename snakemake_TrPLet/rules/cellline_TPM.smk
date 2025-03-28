rule gene_TPM:
    input:
        isoform_quants=expand("results/step1/rsem_results/{Sample}/{Sample}.genes.results", Sample=Samples),
    output:
        "results/step1/gene_TPM_matrix.txt"
    log:
        "logs/gene_TPM.log"
    script:
        "../scripts/rsem_gene_TPM_Matrix.py"
