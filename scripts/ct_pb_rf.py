import scanpy as sc
import pandas as pd
import numpy as np

from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler

def run_ct_pb_rf(adata, sample_key, condition_key, n_splits, params, label_key, **kwargs):
    from utils import custom_pseudobulk_aggregation
    
    adata.obs[condition_key] = adata.obs[condition_key].astype('category')
    all_cts = set(adata.obs[label_key].cat.categories)
    all_cts_sorted = sorted(list(all_cts))

    # Use custom pseudobulk aggregation to avoid decoupler compatibility issues
    adata_ = custom_pseudobulk_aggregation(adata, sample_key=sample_key, groups_col=label_key)
    
    rename_dict = {name: number for number, name in enumerate(np.unique(adata_.obs[condition_key]))}

    if params['norm'] is True:
        scaler = StandardScaler()
        adata_.X = scaler.fit_transform(adata_.X)
    adata_.obs[condition_key] = adata_.obs[condition_key].astype('category')
    adata_.obs[sample_key] = adata_.obs[sample_key].astype('category')
    adata_.obs[label_key] = adata_.obs[label_key].astype('category')

    df = {}
    for donor in adata_.obs[sample_key].cat.categories:
        df[donor] = {}
        for ct in adata_.obs[label_key].cat.categories:
            tmp = adata_[adata_.obs[sample_key] == donor].copy()
            tmp = tmp[tmp.obs[label_key] == ct].copy()
            if len(tmp) > 0:
                df[donor][ct] = tmp.X[0]
            else:
                # If no cells for this donor-celltype combination, use zeros
                df[donor][ct] = np.zeros(adata_.shape[1])
    
    # Ensure all donors have all cell types
    for donor in df.keys():
        for ct in all_cts_sorted:
            if ct not in df[donor]:
                df[donor][ct] = np.zeros(adata_.shape[1])
    
    # Convert to DataFrame with consistent structure
    df_processed = {}
    for donor in df.keys():
        donor_data = []
        for ct in all_cts_sorted:
            donor_data.append(df[donor][ct])
        df_processed[donor] = np.concatenate(donor_data)
    
    df = pd.DataFrame(df_processed)
    df = df.reindex(sorted(df.columns), axis=1)
    df = df.T
    adata_ = sc.AnnData(df)
    
    obs_to_keep = [sample_key, condition_key]
    obs_to_keep.extend([f'split{i}' for i in range(n_splits)])
    obs = adata.obs[obs_to_keep].drop_duplicates().sort_values(sample_key).set_index(sample_key)
    adata_.obs = adata_.obs.join(obs)
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
        
        df = classification_report(y_val, y_pred, output_dict=True)
        df = pd.DataFrame(df).T
        df['split'] = i
        df['method'] = 'ct_pb_rf'
        dfs.append(df)

        val_accuracy = df["f1-score"]["accuracy"]

        val_accuracies.append(val_accuracy)
        val_avg.append(df["f1-score"]["weighted avg"])
        
        
        print(f'Val accuracy = {val_accuracy}.')
        print('===========================')

    df = pd.concat(dfs)
    print(f"Mean validation accuracy across 5 CV splits for a random forest model = {np.mean(np.array(val_accuracies))}.")
    print(f"Mean validation weighted avg across 5 CV splits for a random forest model = {np.mean(np.array(val_avg))}.")
    return df
