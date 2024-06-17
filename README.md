# CancerDepPrediction

**Cancer dependency prediction from RNA (CaDeRNA)**

If used, please cite: _"A landscape and prediction of dependencies in TFE3 fusion-driven cancers"_ (Li* & Sadagopan* et al.)

## Summary

This repo provides the scripts and workflow to accurately predict cancer dependency scores from tumor or cell-line RNA-seq data for a subset of highly predictable genes (N=648). Although you can predict dependency scores for all genes, the accuracy will be substantially lower since most genetic dependencies are not predictable from RNA-seq data alone. The overall workflow involves taking isoform-level RNA count data, batch correcting it (if needed), merging it with a large RNA-sequencing dataset (cell lines: DepMap / CCLE, tumors: TCGA), normalizing RNA-seq counts, reducing dimensionality (by subsetting features), then predicting dependencies with support vector regression.

The calculations you will perform on a merged, batch-corrected isoform-level count matrix during normalization are the following:

1) Calculate transcripts per kilobase million (TPM) for each isoform
2) Calculate gene TPM by summing isoform TPM per gene
3) Convert TPM to log<sub>2</sub>(TPM+1) to generate log-normal distributions
4) Z-score each feature (i.e. each gene)

## Workflow Considerations

The workflow differs slightly depending on input data. There are three major considerations:

### Cell Line or Tumor RNA-Seq Data

If using cell line RNA-seq, you will merging your data with CCLE count or TPM . If using tumor RNA-seq data, you will merge your data with TCGA count data (available here: https://osf.io/gqrz9/files/osfstorage). 

### Batch Correction

ComBat-seq

### Inclusion of Mutation Calls

Mutation data is binarized then the features are Z-scored. 
