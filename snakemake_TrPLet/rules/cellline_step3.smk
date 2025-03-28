rule dep_prediction:
    input:
        "results/step2/genelog2TPMp1_Zscored_with_DepMap_Celllines.tsv.gz",
        "data/Expression_Internal_23Q2.csv",
        "data/CRISPR_DepMap_Internal_23Q2_Score_Chronos.csv",
        "data/summary_sheet_bestmodel_bestMpermodel_noKNN.csv",
    output:
        "results/step3/dep_predictions_Zscoredlog2TPMp1_bestMODEL_all_genes_geq0.2.csv",
    resources:
        mem_mb=256000
    log:
        "logs/dep_prediction.log"
    script:
        "../scripts/cellline_step3_dep_predict.py"
