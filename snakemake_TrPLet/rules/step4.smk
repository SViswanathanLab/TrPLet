rule normalize:
    input:
        "results/step3/batch_corrected_TCGA_counts_with_lineageCOV.tsv.gz",
    output:
        "results/step4/genelog2TPMp1_Zscored_batchCorr_with_TCGA_lineageCOV_final_no_NORMALS.tsv.gz",
    threads: 30
    params:
        threads=30,
        num=num_sample,
    resources:
        mem_mb=64000
    log:
        "logs/normalize.log"
    script:
        "../scripts/step4_normalize.py"
