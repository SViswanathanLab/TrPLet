rule merge_DepMap:
    input:
        DepMap="data/Expression_Internal_23Q2.csv",
        sample="results/step1/gene_log2_TPM+1_matrix.txt",

    output:
        "results/step2/genelog2TPMp1_Zscored_with_DepMap_Celllines.tsv.gz",
    log:
        "logs/merge_DepMap.log"
    threads: 30
    params:
        threads=30
    resources:
        mem_mb=64000
    script:
        "../scripts/cellline_step2_zscore.py"
