rule dep_prediction:
    input:
        "results/step4/genelog2TPMp1_Zscored_batchCorr_with_TCGA_lineageCOV_final_no_NORMALS.tsv.gz",
        "data/Expression_Internal_23Q2.csv",
        "data/CRISPR_DepMap_Internal_23Q2_Score_Chronos.csv",
        "data/summary_sheet_bestmodel_bestMpermodel.csv",
    output:
        "results/step5_to_7/dep_predictions_Zscoredlog2TPMp1_bestMODEL_all_genes_geq0.2.csv",
    threads: 30
    params:
        threads=30,
    resources:
        mem_mb=256000
    log:
        "logs/dep_prediction.log"
    script:
        "../scripts/step5_to_7_dep_predict.py"

# rule dep_prediction:
#     input:
#         "results/step4/Wang_genelog2TPMp1_Zscored_batchCorr_with_TCGA_lineageCOV_final_no_NORMALS.tsv.gz",
#         "data/Expression_Internal_23Q2.csv",
#         "data/CRISPR_DepMap_Internal_23Q2_Score_Chronos.csv",
#         "data/summary_sheet_bestmodel_bestMpermodel.csv",
#     output:
#         "results/step5_to_7/dep_predictions_Zscoredlog2TPMp1_bestMODEL_all_genes_geq0.2.csv",
#     threads: 30
#     params:
#         threads=30,
#     resources:
#         mem_mb=256000
#     log:
#         "logs/dep_prediction.log"
#     script:
#         "../scripts/step5_to_7_dep_predict.py"
