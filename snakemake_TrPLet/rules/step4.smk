rule normalize:
    input:
        "results/step3/batch_corrected_TCGA_counts_with_lineageCOV.tsv.gz",
    output:
        "results/step4/genelog2TPMp1_Zscored_batchCorr_with_TCGA_lineageCOV_final_no_NORMALS.tsv.gz",
    params:
        num=num_sample,
    resources:
        mem_mb=256000,
    log:
        "logs/normalize.log"
    script:
        "../scripts/step4_normalize.py"
