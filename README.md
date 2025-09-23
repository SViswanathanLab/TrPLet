# TrPLet: Cancer Dependency Prediction from RNA-Seq

**Transcriptional Prediction of Lethality (TrPLet)**

If used, please cite both:

_"Assembling a landscape of vulnerabilities across rare kidney cancers"_ (Li* & Sadagopan* et al.), bioRxiv 2024

_"Cancer target discovery enabled by transcriptome-based virtual CRISPR screening"_ (Sadagopan*, B. Li*, J. Li* et al.), bioRxiv 2025

## Summary

This repo provides the scripts and workflow to predict cancer dependency scores from tumor or cell-line RNA-seq data for a subset of highly predictable genes (N=657). Although you can predict dependency scores for all genes, the accuracy will be substantially lower since most genetic dependencies are not predictable from RNA-seq data alone. The most general workflow involves:

1. Generate/download isoform-level* RNA count data (e.g. RNA fastq -> bam -> counts using STAR/RSEM)
2. Merge your data with a large RNA-seq dataset (cell lines: DepMap/CCLE, tumors: TCGA)
3. Batch correct your data, read batch correction section if you are considering this
4. Normalize RNA-seq counts
- Calculate transcripts per kilobase million (TPM) for each isoform
- Calculate gene TPM by summing isoform TPM per gene
- Convert TPM to log<sub>2</sub>(TPM+1) to generate log-normal distributions
- Z-score each feature (i.e. expression of each gene)
5. Reduce dimensionality (subset the train+test data to the top M features with the highest |Pearson correlation coefficient| to the dependency being predicted in the train data; by default M=5000)
6. Predict dependencies on your sample's normalized RNA-seq data

*gene level count data can also be used

The model is trained on the entirety of DepMap and tested on your dataset. When assessing model performance, we used 5-fold cross-validation across DepMap.

## Workflow Considerations

The workflow differs slightly depending on input data. There are three major considerations:

### Cell Line or Tumor RNA-Seq Data

- If predicting on __TCGA tumor RNA-seq__, you can use gene log<sub>2</sub>(TPM+1), which can be calculated from count data here: https://osf.io/gqrz9/files/osfstorage. Z-score the expression of each gene and continue at step 5.
- If predicting on __Non-TCGA tumor RNA-seq__, you should merge your data with TCGA (batch correction is likely necessary); if batch correcting, you will merge the external dataset with TCGA isoform-level count data (available here: https://osf.io/gqrz9/files/osfstorage) and continue at step 2 of the workflow.
- If predicting on __Cell Line RNA-seq__, you can use gene log<sub>2</sub>(TPM+1) available from here: https://depmap.org/portal/data_page/?tab=currentRelease). Z-score the expression of each gene using the mean/standard deviation from DepMap, and continue at step 5. For this use case, batch correction usually isn't necessary. If it is required, merge with CCLE isoform-level counts (available here: https://osf.io/gqrz9/files/osfstorage) and continue at step 2 of the workflow.

### Batch Correction

In general, we noticed that cell line RNA-seq often does not require batch correction with CCLE, while tumor RNA-seq with TCGA almost always does (based on tSNE analysis). We recommend batch correcting using ComBat-seq using lineage as a covariate. Choose a lineage of your sample that most closely matches with TCGA or CCLE lineages. TCGA lineages available here: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations; CCLE lineages are indicated by the string following the underscore "_" in the cell line name as indicated in the sample_info.csv file (https://depmap.org/portal/data_page/?tab=allData).

### Computing Isoform-Level Counts

We recommend using STAR/RSEM, we have not tested other methods of quantification, though they may also work.


