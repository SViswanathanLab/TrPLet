import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.stats import zscore

df = pd.read_csv("batch_corrected_Wang_TCGA_counts_with_lineageCOV.tsv.gz", compression='gzip', header=0, sep='\t', quotechar='"')
print(df)

sample_columns = df.columns.tolist()[3:]
temp_df = df[sample_columns]
temp_df.replace([np.inf, -np.inf], np.nan, inplace=True)
temp_df = temp_df.fillna(0)
print(temp_df)

#temp_row = pd.Series(df['effective_length'].tolist())
temp_row = df['effective_length']
temp_df = temp_df.div(temp_row, axis=0)
print(temp_df)

cols = temp_df.columns.tolist()

for a in cols:
    sub_val = temp_df[a].tolist()
    temp_val = []
    for b in sub_val:
        if b != b or b == float("inf"):
            temp_val.append(0)
        else:
            temp_val.append(float(b))
    sum_val = sum(temp_val)
    new_vals = [x*1000000/sum_val for x in temp_val]
    temp_df[str(a)] = new_vals

#temp_df = temp_df.div(temp_df.sum())*10**6
print(temp_df)

genes_unp = df['genes'].tolist()
new_genes = []
for a in genes_unp:
    new_genes.append(a.split("-")[0])

temp_df['gene'] = new_genes
temp_df = temp_df.groupby('gene').mean().reset_index()
print(temp_df)

gene_list = temp_df['gene'].tolist()
del temp_df['gene']

temp_df = temp_df + 1
temp_df = np.log2(temp_df)
print(temp_df)

vals = temp_df.values.tolist()
new_vals = []
for a in vals:
    try:
        new_val = zscore(a)
    except:
        new_val = [0]*len(a)
    new_vals.append(new_val)

new_df = pd.DataFrame(new_vals)
new_df.columns = temp_df.columns.tolist()

print(new_df)

new_df.index = gene_list
new_df = new_df.T
new_df = new_df.head(7) #taking only the first 7 samples (ASPS from Wang) that we are going to predict on

print(new_df)

new_df.to_csv("Wang_genelog2TPMp1_Zscored_batchCorr_with_TCGA_lineageCOV_final_no_NORMALS.tsv.gz", sep="\t", compression="gzip")
