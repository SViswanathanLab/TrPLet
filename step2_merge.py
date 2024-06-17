import os
import pandas as pd
import gzip

def round_numeric_columns(df, decimals):
    numeric_cols = df.select_dtypes(include=[float, int]).columns
    df[numeric_cols] = df[numeric_cols].round(decimals)
    return df

path_to_your_data = "ASPS_data.txt" #replace with path to your isoform-level count data
path_to_TCGA_count_matrices = "/home/au409/TCGA_count_matrices/" #replace to path with TCGA count matrices downloaded from here: https://osf.io/gqrz9/files/osfstorage

files = sorted([x for x in os.listdir(path_to_TCGA_count_matrices) if "html" in x])
file0 = files[0]
files = files[1:]
out_df = pd.read_csv(path_to_TCGA_count_matrices + file0, compression='gzip', header=0, sep='\t', quotechar='"')

for a in files:
        df = pd.read_csv(path_to_TCGA_count_matrices + a, compression='gzip', header=0, sep='\t', quotechar='"')
        out_df = out_df.merge(df,how='outer',on="Unnamed: 0")

out_df.to_csv(path_to_TCGA_count_matrices + "merged_TCGA_counts.tsv.gz", index=False, sep="\t", compression="gzip")

#df = pd.read_csv(path_to_TCGA_count_matrices + "merged_TCGA_counts.tsv.gz", compression='gzip', header=0, sep='\t', quotechar='"')
df = out_df

df_ids = df['Unnamed: 0'].tolist()
#ID format: ENST00000000233.9|ENSG00000004059.10|OTTHUMG00000023246.6|OTTHUMT00000059567.2|ARF5-001|ARF5|1103|protein_coding|

df2 = pd.read_csv(path_to_your_data, sep="\t")
df2_ids = df2['transcript_id'].tolist()
#ID format: ENST00000373020.9_TSPAN6-201

df2_genes = []
df2_ids_proc = []
for a in df2_ids:
        df2_ids_proc.append(a.split(".")[0])
        df2_genes.append(a.split("_")[-1])

df_ids_proc = []
for a in df_ids:
        df_ids_proc.append(a.split(".")[0])

df['revised_IDs'] = df_ids_proc
del df['Unnamed: 0']

#df2['revised_IDs'] = df2_ids_proc
#df2['genes'] = df2_genes

df2.insert(0,'genes',df2_genes)
df2.insert(0,'revised_IDs',df2_ids_proc)
del df2['transcript_id']
#del df2['length']

df_merge = df2.merge(df,on='revised_IDs',how='inner')
df_merge = round_numeric_columns(df_merge,0)

df_merge.to_csv("ASPS7_TCGA11373_merged_counts.tsv.gz", index=False, sep="\t", compression="gzip")
