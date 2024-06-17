# Triplet: Cancer Dependency Prediction from RNA-Seq

**Transcriptional Prediction of Lethality (Triplet)**

If used, please cite: _"A landscape and prediction of dependencies in TFE3 fusion-driven cancers"_ (Li* & Sadagopan* et al.)

## Summary

This repo provides the scripts and workflow to accurately predict cancer dependency scores from tumor or cell-line RNA-seq data for a subset of highly predictable genes (N=648). Although you can predict dependency scores for all genes, the accuracy will be substantially lower since most genetic dependencies are not predictable from RNA-seq data alone. The most general workflow involves:

1. Generate isoform-level RNA count data (e.g. RNA fastq -> bam -> counts using STAR/RSEM)
2. (Batch correct and) merge your data with a large RNA-seq dataset (cell lines: DepMap/CCLE, tumors: TCGA), read batch correction section if you plan to batch correct
3. Normalize RNA-seq counts
- Calculate transcripts per kilobase million (TPM) for each isoform
- Calculate gene TPM by summing isoform TPM per gene
- Convert TPM to log<sub>2</sub>(TPM+1) to generate log-normal distributions
- Z-score each feature (i.e. each gene)
4. Reduce dimensionality (subset the train+test data to the top M features with the highest |Pearson correlation coefficient| to the dependency being predicted in the train data; by default M=5000)
5. Train the support vector regression model
6. Predict dependencies on your sample's normalized RNA-seq data

The model is trained on the entirety of DepMap/CCLE. When assessing model performance, we used 5-fold cross-validation.

## Workflow Considerations

The workflow differs slightly depending on input data. There are four major considerations:

### Cell Line or Tumor RNA-Seq Data

- If predicting on __CCLE RNA-seq__, you can directly use normalized CCLE RNA-seq data (gene TPM is fine for this purpose and is available here: https://depmap.org/portal/data_page/?tab=currentRelease; Note: convert this to Z-scored log<sub>2</sub>(TPM+1) then continue at step 4 of the main workflow)
- If predicting on __non-CCLE RNA-seq__, you should merge your data with CCLE (batch correcting, if needed); if you are not batch correcting, merge your sample's gene TPMs with CCLE gene TPMs, compute Z-scored log<sub>2</sub>(TPM+1), then continue at step 4 of the main workflow. If you are batch correcting, you will merge with CCLE isoform-level count data (available here: https://osf.io/gqrz9/files/osfstorage, and continue at step 2 of the workflow)
- If predicting on __TCGA tumor RNA-seq__, you can directly use normalized TCGA RNA-seq data (gene FPKM-UQ is fine for this purpose and is available here: https://gdc.cancer.gov/about-data/publications/pancanatlas (file: EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv); Note: convert this to Z-scored log<sub>2</sub>(FPKM-UQ+1) then continue at step 4 of the main workflow)
- If predicting on __non-TCGA tumor RNA-seq__, you should merge your data with TCGA (batch correction is likely necessary); if batch correcting, you will merge with TCGA isoform-level count data (available here: https://osf.io/gqrz9/files/osfstorage, and continue at step 2 of the workflow)

### Batch Correction

In general, we noticed that cell line RNA-seq often does not require batch correction with CCLE, while tumor RNA-seq with TCGA almost always does (based on two-component PCA). We recommend batch correcting using ComBat-seq using lineage as a covariate. Choose a lineage of your sample that most closely matches with TCGA or CCLE lineages. TCGA lineages available here: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations; CCLE lineages are indicated by the string following the underscore "_" in the cell line name as indicated in the sample_info.csv file (https://depmap.org/portal/data_page/?tab=allData).

### Inclusion of Mutation Calls

This is completely optional; including mutation calls may increase the accuracy of predicting dependency scores for a few genes. If you choose to inclue mutation data, it should be binarized (0: no mutation, 1: mutation present) then Z-scored. You can add mutation status of a gene manually during the RNA-seq normalization stage. CCLE has mutation calls available here: https://depmap.org/portal/data_page. TCGA has mutation calls available here: https://gdc.cancer.gov/about-data/publications/pancanatlas. One should filter to high impact / deleterious / hotspot mutations.

If you wish to call mutations from RNA-seq data, consider using an established workflow such as GATK's: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels.

### Computing Isoform-Level Counts

We recommend using STAR/RSEM, we have not tested other methods of quantification, though they may also work.


