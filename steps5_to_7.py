import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
import xgboost as xgb
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
import torch
import torch.nn as nn
import torch.optim as optim
from scipy.stats import pearsonr
from sklearn.linear_model import Lasso
from sklearn.preprocessing import StandardScaler
import random
import math
from operator import itemgetter
import statistics

#Step 1: Merge RNA-seq and mutation calls

dam_mut_df = pd.read_csv("Hotspot_Mutations.csv") #Downloadable from DepMap Portal
dam_mut_df['depmap_id'] = dam_mut_df['Unnamed: 0']
del dam_mut_df['Unnamed: 0']

new_cols = []
cols = dam_mut_df.columns.tolist()
for a in cols:
    if a != 'depmap_id':
        new_cols.append(a + "_y")
    else:
        new_cols.append(a)
dam_mut_df.columns = new_cols

exp_df = pd.read_csv("Expression_Public_23Q2.csv") #Downloadable from DepMap Portal
del exp_df['cell_line_display_name']
del exp_df['lineage_2']
del exp_df['lineage_3']
del exp_df['lineage_4']
del exp_df['lineage_5']
del exp_df['lineage_6']

new_cols = []
cols = exp_df.columns.tolist()
for a in cols:
    if a != 'depmap_id' and a != 'lineage_1':
        new_cols.append(a + "_x")
    else:
        new_cols.append(a)
exp_df.columns = new_cols

#Uncomment the line below if you want mutation calls included
merged_df = exp_df #.merge(dam_mut_df, how='inner', on='depmap_id')
lineages = merged_df['lineage_1'].tolist()
uq_lineages = sorted(list(set(lineages)))
numeric_lineages = []
for a in lineages:
    numeric_lineages.append(uq_lineages.index(a))
del merged_df['lineage_1']
merged_df['lineage'] = numeric_lineages
RNA_IDs = merged_df['depmap_id'].tolist()

#Step 2: Read the Chronos scores and match order between IV and DV dataframes

df2 = pd.read_csv("CRISPR_DepMap_Public_23Q2_Score_Chronos.csv") #Downloadable from DepMap Portal
IDs = df2['Unnamed: 0'].tolist()
shared_ids = list(set(RNA_IDs) & set(IDs))
print("Number of cell lines with complete data: " + str(len(shared_ids)))

merged_df = merged_df[merged_df['depmap_id'].isin(shared_ids)]
df2 = df2[df2['Unnamed: 0'].isin(shared_ids)]
df = merged_df.sort_values('depmap_id')
df2 = df2.sort_values('Unnamed: 0')

#Step 3: Separate lineage column and final preprocessing steps

final_IDs_list = df['depmap_id'].tolist()
del df['depmap_id']
del df2['Unnamed: 0']

lineages = df['lineage'].tolist()
del df['lineage']

uq_lineages = sorted(list(set(lineages)))
lineage_cols = []

for a in uq_lineages:
    col = []
    for b in lineages:
        if a == b:
            col.append(1)
        else:
            col.append(0)
    lineage_val = 'Lineage' + str(a)
    lineage_cols.append(lineage_val)
    df[lineage_val] = col

#Step 4: Z-score the RNA-seq, hotspot mutation calls, and lineage annotations

testing_df = pd.read_csv("Wang_genelog2TPMp1_Zscored_batchCorr_with_TCGA_lineageCOV_final_no_NORMALS.tsv.gz", compression='gzip', header=0, sep='\t', quotechar='"') #Generated in step 4
testing_df = testing_df.dropna(axis=1, how='all')
names_to_test = testing_df['Unnamed: 0'].tolist()
del testing_df['Unnamed: 0']

testing_df_cols = testing_df.columns.tolist()
new_cols = []
for a in testing_df_cols:
    new_cols.append(a + "_x")
testing_df.columns = new_cols

shared_cols = list(set(testing_df.columns.tolist()) & set(df.columns.tolist()))
df = df[shared_cols]
testing_df = testing_df[shared_cols]
scaler = StandardScaler()
df = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

#Step 5: Predict highly predictable DVs

M5000_df = pd.read_csv("SVR_performance/predicted_vs_observed_correl_coeffs_Chronos_5foldCV_N5000_SVR_testdata.csv") #available in this repo
top_df = M5000_df[M5000_df['Correl_Coeffs_5CV']>=0.4]
vals_to_test = top_df['Gene_Name'].tolist() #change vals_to_test to whatever genes you want to test (to force inclusion of certain genes)

N=5000 #Number of features to restrict to

print(len(vals_to_test))

model = SVR()
#model = LinearRegression()
#model = GradientBoostingRegressor()
#model = xgb.XGBRegressor()
#model = DecisionTreeRegressor()
#model = KNeighborsRegressor()
#model = Lasso(alpha=0.01)
#model = RandomForestRegressor(n_estimators=100, random_state=42)

correl_coeffs = []
abc=0

output_vals = []

for q in vals_to_test:

    try:
        dv_list = df2[q].tolist()
    except:
        output_vals.append([])
        print("Gene does not have dependency data")
        continue

    X_train = df
    y_train = dv_list
    X_test = testing_df

    #Restrict to top N features
    cols = X_train.columns.tolist()

    corr_list = []
    aa=0
    while aa<len(cols):
        corr_list.append(abs(np.corrcoef(X_train[cols[aa]].tolist(),y_train)[0, 1]))
        aa=aa+1

    sublist = [x for x in corr_list if x == x]
    sorted_list = sorted(sublist, reverse=True)
    nth_highest = sorted_list[N-1]

    top_features = []

    aa=0
    while aa<len(corr_list):
        if corr_list[aa] >= nth_highest:
            top_features.append(cols[aa])
        aa=aa+1

    X_train = X_train[top_features]
    X_test = X_test[top_features]

    try:
        model.fit(X_train, y_train)
    except:
        output_vals.append([])
        print("Model unfittable")
        continue

    """
    #Output coefficients:
    support_vectors = model.support_vectors_
    support_vector_indices = model.support_

    dual_coef = np.transpose(model.dual_coef_)
    coefficients = np.sum(dual_coef * support_vectors, axis=0)

    coefficients_df = pd.DataFrame({"Features": top_features, "Coefficients": coefficients})
    coefficients_df.to_csv("coefficients.csv")
    """

    avg_value = [statistics.mean(y_train)]
    y_pred = model.predict(X_test)
    output_vals.append(list(y_pred) + list(avg_value))
    print(abc)
    abc=abc+1

output_df = pd.DataFrame(output_vals)
output_df.columns = names_to_test + ['CCLE_mean']
output_df['Gene'] = vals_to_test
output_df.to_csv("tRCC_dep_predictions_SVR_Wang_batchCorrwTCGA_lineageCOV_Zscored_RNAseq_M5000.csv", index=False)
