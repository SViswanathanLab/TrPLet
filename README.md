# CancerDepPrediction

**Cancer dependency prediction from RNA (CaDeRNA)**

If used, please cite: _"A landscape and prediction of dependencies in TFE3 fusion-driven cancers"_ (Li* & Sadagopan* et al.)

## Summary

This repo provides the scripts and workflow to accurately predict cancer dependency scores from tumor or cell-line RNA-seq data for a subset of highly predictable genes (N=648). Although you can predict dependency scores for all genes, the accuracy will be substantially lower since most genetic dependencies are not predictable from RNA-seq data alone. The overall workflow involves taking isoform-level RNA count data, batch correcting it (if needed), merging it with a large RNA-sequencing dataset (cell lines: DepMap / CCLE, tumors: TCGA), normalizing RNA-seq counts, reducing dimensionality (by subsetting the train+test data to the top M features with the highest |Pearson correlation coefficient| to the dependency being predicted using the train data, by default M=5000), training the support vector regression model, and lastly predicting dependencies on your Z-scored normalized merged batch-corrected RNA-seq data. The data the model is trained on is the entirety of DepMap. When assessing model performance, we used 5-fold cross-validation, training on 80% of DepMap and assessing performance using the remaining 20% (then repeating).

The calculations you will perform on a merged, batch-corrected isoform-level count matrix during normalization are the following:

1) Calculate transcripts per kilobase million (TPM) for each isoform
2) Calculate gene TPM by summing isoform TPM per gene
3) Convert TPM to log<sub>2</sub>(TPM+1) to generate log-normal distributions
4) Z-score each feature (i.e. each gene)

## Workflow Considerations

The workflow differs slightly depending on input data. There are four major considerations:

### Cell Line or Tumor RNA-Seq Data

If using cell line RNA-seq, you will merging your data with CCLE count or TPM ???
. If using tumor RNA-seq data, you will merge your data with TCGA count data (available here: https://osf.io/gqrz9/files/osfstorage). 

### Batch Correction

ComBat-seq with lineage co-variate. Make sure you use isoform-level counts

### Inclusion of Mutation Calls

This is completely optional; including mutation calls may increase the accuracy of predicting dependency scores for a few genes. If you choose to inclue mutation data, it should be binarized (0: no mutation, 1: mutation present) then Z-scored. You can add mutation status of a gene manually during the RNA-seq normalization stage. CCLE has mutation calls available here: https://depmap.org/portal/data_page. TCGA has mutation calls available here: https://gdc.cancer.gov/about-data/publications/pancanatlas. One should filter to high impact / deleterious / hotspot mutations if possible.

If you wish to call mutations from RNA-seq data, please use an established workflow. One possible method is using HaplotypeCaller after generating a processed RNA-seq .bam file, then some variant filtering method. One could consider the GATK workflow: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels. 

### Computing Isoform-Level Counts

We recommend using STAR/RSEM, we have not tested with other methods.



