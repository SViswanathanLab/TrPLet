rule merge_TCGA:
    input:
        TCGA="data/merged_TCGA_geneLevel_CountMatrix.tsv.gz",
        sample="data/gene_count_matrix.txt",

    output:
        "results/step2/merged_with_TCGA11373_counts.tsv.gz",
    log:
        "logs/merge_TCGA.log"
    resources:
        mem_mb=64000
    script:
        "../scripts/step2_merge.py"
