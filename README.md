# CancerDepPrediction

**Cancer dependency prediction from RNA (CaDeRNA)**

If used, please cite: _"A landscape and prediction of dependencies in TFE3 fusion-driven cancers"_ (Li* & Sadagopan* et al.)

## Summary

This repo provides the scripts and workflow to accurately predict cancer dependency scores from tumor or cell-line RNA-seq data for a subset of highly predictable genes (N=648). Although you can predict dependency scores for all genes, the accuracy will be substantially lower since most genetic dependencies are not predictable from RNA-seq data alone. The overall workflow is to take RNA count data, batch correct it (if needed), merge it with a large RNA-sequencing dataset (cell lines: DepMap, tumors: TCGA), Z-score your data with the RNA-seq dataset, reduce dimensionality (by subsetting features), then predict dependencies with support vector regression.

## Specific Workflow

workflow differs slightly depending on input data. There are a few things to consider prior to implementing the general procedure. Your goal is to get your RNA-seq data batch corrected with either TCGA or CCLE 

- Batch Correction Considerations

- Mutation Data (optional)

- Cell Line RNA-Seq (no batch correction)

- Cell Line RNA-Seq (with batch correction)
