# install packages if needed
import sys
import subprocess
import os

# load required packages
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.stats import zscore

df = pd.read_csv(snakemake.input[0], compression='gzip', header=0, sep='\t', quotechar='"')
print(df)

sample_columns = df.columns.tolist()[3:]
temp_df = df[sample_columns]
temp_df = temp_df.apply(pd.to_numeric, errors="coerce")
temp_df.replace([np.inf, -np.inf], np.nan, inplace=True)
temp_df = temp_df.fillna(0)
print(temp_df)

# temp_row = df['effective_length']
temp_row = pd.to_numeric(df['effective_length'], errors="coerce")

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

print(temp_df)

gene_list = df['genes'].tolist()

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
new_df = new_df.head(snakemake.params['num'])

print(new_df)

os.makedirs("results/step4", exist_ok=True) # Create the directory results/step2
new_df.to_csv(snakemake.output[0], sep="\t", compression="gzip")
