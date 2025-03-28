import sys
# log file
sys.stderr = open(snakemake.log[0], "w")
import pandas as pd
import os
import numpy as np


# Initialize an empty DataFrame
merged_count = pd.DataFrame(columns=["gene_id", "expected_count", "Sample"])

# load the genes.results for each sample and merge them together
for quant_path in snakemake.input:
    print(f"Processing file: {quant_path}")
    quant = pd.read_table(quant_path)
    effLen_df = quant[["gene_id", "effective_length"]] # randomly pick effective_length from 1 sample
    quant = quant[["gene_id", "expected_count"]]
    sample_name = os.path.basename(quant_path).replace(".genes.results", "") # extract sample name from path
    quant["Sample"] = sample_name # Add a column of sample names

    merged_count = pd.concat([merged_count, quant], ignore_index=True) # merge quant matrices of all samples

# Pivot the dataframe to create a wide format
gene_count_wide = merged_count.pivot(index=["gene_id"], columns="Sample", values="expected_count")
# Reset the index to make transcript_id a regular column
gene_count_wide.reset_index(inplace=True)

# Add effective_length column to the final count matrix based on the transcript_id
gene_count_wide_withEffLen = pd.merge(gene_count_wide, effLen_df, on='gene_id',how='left')
# move effective_length forward for easier later manipulation
columns = gene_count_wide_withEffLen.columns.tolist()
columns.insert(1, columns.pop(columns.index("effective_length")))
gene_count_wide_withEffLen = gene_count_wide_withEffLen[columns]

# save into txt file
gene_count_wide_withEffLen.to_csv(snakemake.output[0], sep="\t", index=False)


