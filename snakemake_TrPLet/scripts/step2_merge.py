import os
import pandas as pd
import gzip

# log file
sys.stderr = open(snakemake.log[0], "w")

def round_numeric_columns(df, decimals):
    numeric_cols = df.select_dtypes(include=[float, int]).columns
    df[numeric_cols] = df[numeric_cols].round(decimals)
    return df

df = pd.read_csv(snakemake.input[0], compression='gzip', header=0, sep='\t', quotechar='"')

df['revised_IDs'] = df['gene_name']
del df['gene_name'] #ID format: ARF5
#ID format: ENST00000000233.9|ENSG00000004059.10|OTTHUMG00000023246.6|OTTHUMT00000059567.2|ARF5-001|ARF5|1103|protein_coding|

df2 = pd.read_csv(snakemake.input[1], sep="\t")
df2_ids = df2['gene_id'].tolist()
#ID format: 'ENSG00000000005.6_TNMD'

df2_genes = []
df2_ids_proc = []
for a in df2_ids:
    df2_ids_proc.append(a.split(".")[0])
    df2_genes.append(a.split("_")[-1])

df2.insert(0,'genes',df2_genes)
df2.insert(0,'revised_IDs',df2_ids_proc)
del df2['gene_id']

df_merge = df2.merge(df,on='genes',how='inner')
df_merge = round_numeric_columns(df_merge,0)

os.makedirs("results/step2", exist_ok=True) # Create the directory results/step2
df_merge.to_csv(snakemake.output[0], index=False, sep="\t", compression="gzip")
