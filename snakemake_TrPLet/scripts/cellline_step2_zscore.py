import os
import pandas as pd
import gzip
import numpy as np

# log file
sys.stderr = open(snakemake.log[0], "w")


# User input RSEM TPM dataset - gene_log2_TPM+1_matrix.txt
df2 = pd.read_csv(snakemake.input[1], sep="\t") # load

# log2(TPM+1) transformation
numeric_cols = df2.columns[1:]
df2[numeric_cols] = np.log2(df2[numeric_cols] + 1)

df2_ids = df2['gene_id'].tolist()
#ID format: ENSG00000000003.15_TSPAN6
df2_genes = [] # extract the list of gene names
for a in df2_ids:
        df2_genes.append(a.split("_")[1])
df2.insert(0,'genes',df2_genes)
del df2['gene_id']


# DepMap dataset - data/Expression_Internal_23Q2.csv
df = pd.read_csv(snakemake.input[0]) # load DepMap Cellline TPM matrix
subset_df = df.iloc[:, [0]].join(df.iloc[:, 8:]) # keep only columns of depmap_id & TPM for genes
subset_df.set_index('depmap_id', inplace=True) # set 'depmap_id' as the rowname
transposed_df = subset_df.transpose() # swap rows and columns
transposed_df = transposed_df.reset_index() # reset the index so that gene names become the first column
transposed_df.rename(columns={'index': 'gene_id'}, inplace=True) # rename the first column to 'gene_id'
# Get mean&sd for each gene
transposed_df['mean'] = transposed_df.iloc[:, 1:].mean(axis=1)
transposed_df['sd'] = transposed_df.iloc[:, 1:].std(axis=1)

# merge
df_merge = df2.merge(transposed_df,left_on='genes', right_on='gene_id',how='inner')
del df_merge['gene_id']

# z-score transformation
subset_df_merge = df_merge.iloc[:, 0:1]
cols = df2.iloc[:,1:].columns.tolist()
for col in cols:
    subset_df_merge[col] = (df_merge[col] - df_merge['mean']) / df_merge['sd']
subset_df_merge.set_index('genes', inplace=True) # set 'depmap_id' as the rowname
t_subset_df_merge = subset_df_merge.transpose() # swap rows and columns
t_subset_df_merge = t_subset_df_merge.reset_index() # reset the index so that gene names become the first column
t_subset_df_merge.rename(columns={'index': 'Cell_lines'}, inplace=True)


# save
os.makedirs("results/step2", exist_ok=True) # Create the directory results/step2
t_subset_df_merge.to_csv(snakemake.output[0], index=False, sep="\t", compression="gzip")
