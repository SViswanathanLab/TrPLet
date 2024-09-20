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

df_ids = df['Unnamed: 0'].tolist()
#ID format: ENST00000000233.9|ENSG00000004059.10|OTTHUMG00000023246.6|OTTHUMT00000059567.2|ARF5-001|ARF5|1103|protein_coding|

df2 = pd.read_csv(snakemake.input[1], sep="\t")
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



os.makedirs("results/step2", exist_ok=True) # Create the directory results/step2
df_merge.to_csv(snakemake.output[0], index=False, sep="\t", compression="gzip")

