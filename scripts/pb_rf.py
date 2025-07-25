import pandas as pd
import numpy as np

from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler

def run_pb_rf(adata, sample_key, condition_key, n_splits, params, **kwargs):
    from utils import custom_pseudobulk_aggregation
    
    # Use custom pseudobulk aggregation to avoid decoupler compatibility issues
    adata_ = custom_pseudobulk_aggregation(adata, sample_key=sample_key, groups_col=None)
    rename_dict = {name: number for number, name in enumerate(np.unique(adata_.obs[condition_key]))}

    if params['norm'] is True:
        scaler = StandardScaler()
        adata_.X = scaler.fit_transform(adata_.X)
    adata_.obs[condition_key] = adata_.obs[condition_key].astype('category')

    val_accuracies = []
    val_avg = []
    dfs = []

    for i in range(n_splits):
        print(f'Processing split = {i}...')
        df = adata.obs[[f'split{i}', sample_key]].drop_duplicates()
        train = list(df[df[f'split{i}'] == 'train'][sample_key])
        val = list(df[df[f'split{i}'] == 'val'][sample_key])
        # train data
        x = pd.DataFrame(adata_[adata_.obs_names.isin(train)].X).to_numpy()
        y = adata_[adata_.obs_names.isin(train)].obs[condition_key].cat.rename_categories(rename_dict)
        y = y.to_numpy()
        print("Train shapes:")
        print(f"x.shape = {x.shape}")
        print(f"y.shape = {y.shape}")
        # val data
        x_val = pd.DataFrame(adata_[adata_.obs_names.isin(val)].X).to_numpy()
        y_val = adata_[adata_.obs_names.isin(val)].obs[condition_key].cat.rename_categories(rename_dict)
        y_val = y_val.to_numpy()
        print("Val shapes:")
        print(f"x_val.shape = {x_val.shape}")
        print(f"y_val.shape = {y_val.shape}")
        # fit
        X = x
        Y = y
        clf = RandomForestClassifier()
        clf.fit(X, Y)
        print(f'Train accuracy = {np.sum(clf.predict(X) == Y)/len(Y)}.')
        y_pred = clf.predict(x_val)
        val_accuracy = np.sum(y_pred == y_val)/len(y_val)
        
        df = classification_report(y_val, y_pred, output_dict=True)
        df = pd.DataFrame(df).T
        df['split'] = i
        df['method'] = 'pb_rf'
        dfs.append(df)
        val_accuracies.append(df["f1-score"]["accuracy"])
        val_avg.append(df["f1-score"]["weighted avg"])
        
        val_accuracy = df["f1-score"]["accuracy"]
        
        print(f'Val accuracy = {val_accuracy}.')
        print('===========================')

    df = pd.concat(dfs)
    print(f"Mean validation accuracy across 5 CV splits for a random forest model = {np.mean(np.array(val_accuracies))}.")
    print(f"Mean validation weighted avg across 5 CV splits for a random forest model = {np.mean(np.array(val_avg))}.")
    return df
