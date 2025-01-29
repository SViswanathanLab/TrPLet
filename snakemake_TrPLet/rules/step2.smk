rule merge_TCGA:
    input:
        TCGA="data/merged_TCGA_counts.tsv.gz",
        sample="results/step1/isoform_count_matrix.txt",

    output:
        "results/step2/merged_with_TCGA11373_counts.tsv.gz",
    log:
        "logs/merge_TCGA.log"
    threads: 30
    params:
        threads=30
    resources:
        mem_mb=64000
    script:
        "../scripts/step2_merge.py"
