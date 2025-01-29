import pandas as pd
from sklearn.svm import SVR
import numpy as np
from scipy.stats import pearsonr
from sklearn.preprocessing import StandardScaler
import statistics
import random
import math
from operator import itemgetter
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler
import os



#Step 1: Deleted
exp_df = pd.read_csv(snakemake.input[1])
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

merged_df = exp_df #.merge(dam_mut_df, how='inner', on='depmap_id')
lineages = merged_df['lineage_1'].tolist()
uq_lineages = sorted(list(set(lineages)))
numeric_lineages = []
for a in lineages:
    numeric_lineages.append(uq_lineages.index(a))
del merged_df['lineage_1']
merged_df['lineage'] = numeric_lineages
RNA_IDs = merged_df['depmap_id'].tolist()



#Step 3: Read the Chronos scores and match order between IV and DV dataframess
df2 = pd.read_csv(snakemake.input[2])
IDs = df2['Unnamed: 0'].tolist()
shared_ids = [x for x in list(set(RNA_IDs) & set(IDs)) if x == x]
RNA_specific_IDs = [x for x in list(np.setdiff1d(RNA_IDs, IDs)) if x == x]

print("Number of cell lines with complete data: " + str(len(shared_ids)))
print("Number of cell lines with RNA only: " + str(len(RNA_specific_IDs)))

merged_df = merged_df[merged_df['depmap_id'].isin(shared_ids)]
df2 = df2[df2['Unnamed: 0'].isin(shared_ids)]
df = merged_df.sort_values('depmap_id')
df2 = df2.sort_values('Unnamed: 0')



#Step 4: Separate lineage column and final preprocessing steps
final_IDs_list = df['depmap_id'].tolist()
del df['depmap_id']
del df2['Unnamed: 0']

lineages = df['lineage'].tolist()
del df['lineage']

"""
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
"""

#generate vals to test

#M5000_df = pd.read_csv("/home/au409/TCGA_count_matrices/predicted_vs_observed_correl_coeffs_Chronos_5foldCV_N5000_SVR_testdata.csv")
#top_df = M5000_df[M5000_df['Correl_Coeffs_5CV']>=0.2]
#vals_to_test = top_df['Gene_Name'].tolist()

#Step 5: Z-score the RNA-seq, hotspot mutation calls, and lineage annotations, merge with tRCC RNA
testing_df = pd.read_csv(snakemake.input[0], compression='gzip', header=0, sep='\t', quotechar='"')
tumors_tested_on = testing_df['Unnamed: 0'].tolist()
del testing_df['Unnamed: 0']

testing_df = testing_df.fillna(0)

gene_id_unp = testing_df.columns.tolist()
gene_id = []
for a in gene_id_unp:
    gene_id.append(str(a) + "_x")
testing_df.columns = gene_id

testing_df_cols = testing_df.columns.tolist()
df_cols = df.columns.tolist()
common_cols = [x for x in list(set(testing_df_cols) & set(df_cols)) if x == x]

testing_df = testing_df[common_cols]
df = df[common_cols]

scaler = StandardScaler()
df = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

#Step 6: Predict  DVs restricting to those with features in our dataset

best_model_df = pd.read_csv(snakemake.input[3])
best_model_df = best_model_df[best_model_df['Correl_Coeffs_5CV']>=0.2]
vals_to_test = best_model_df['Gene_Name'].tolist()
model_types_list = best_model_df['Model_Type'].tolist()

print(len(vals_to_test))

#model = SVR() #R-bar = 0.50
#model = SVR(kernel='linear')
#model = LinearRegression()
#model = GradientBoostingRegressor()
#model = xgb.XGBRegressor()
#model = DecisionTreeRegressor() #0.32
#model = KNeighborsRegressor() #0.42
#model = Lasso(alpha=0.01) #0.45
#model = RandomForestRegressor(n_estimators=100, random_state=42)

correl_coeffs = []
abc=0

output_vals = []

for q,temp_model_type in zip(vals_to_test,model_types_list):

    try:
        dv_list = df2[q].tolist()
    except:
        output_vals.append([])
        print("Gene does not have dependency data")
        continue

    model_portion = temp_model_type.split("_")[1]
    N = int(temp_model_type.split("_")[0].split("N")[-1])

    if model_portion == "linear":
        model = SVR(kernel='linear')
    elif model_portion == "SVR":
        model = SVR()
    elif model_portion == "ridge":
        model = Ridge(alpha=0.01)
    elif model_portion == "lasso":
        model = Lasso(alpha=0.01)
    elif model_portion == "knn":
        model = KNeighborsRegressor()
    elif model_portion == "elasticnet":
        model = ElasticNet(alpha=0.01)
    else:
        print("ERROR NO MODEL UPDATED")
        print(q)

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
    # Get the support vectors
    support_vectors = model.support_vectors_
    support_vector_indices = model.support_

    # Get the support vector coefficients
    dual_coef = np.transpose(model.dual_coef_)

    # Compute the approximate coefficients
    coefficients = np.sum(dual_coef * support_vectors, axis=0)

    # Create a DataFrame with features and coefficients
    coefficients_df = pd.DataFrame({"Features": top_features, "Coefficients": coefficients})
    coefficients_df.to_csv("/Users/ananthansadagopan/Downloads/%s_coefficients_N1_linear.csv" % (q))
    """

    avg_value = [statistics.mean(y_train)]
    y_pred = model.predict(X_test)
    output_vals.append(list(y_pred) + list(avg_value))
    print(abc)
    abc=abc+1

    if abc % 100 == 0:
        output_df = pd.DataFrame(output_vals)
        output_df.columns = tumors_tested_on + ['CCLE_mean']
        output_df['Gene'] = vals_to_test[0:abc]
        output_df.to_csv(snakemake.output[0], index=False)

output_df = pd.DataFrame(output_vals)
output_df.columns = tumors_tested_on + ['CCLE_mean']
output_df['Gene'] = vals_to_test
output_df.to_csv(snakemake.output[0], index=False)
