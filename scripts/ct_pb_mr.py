import scanpy as sc
import pandas as pd
import numpy as np

from sklearn.metrics import classification_report
from sklearn.preprocessing import OneHotEncoder
from scipy.special import softmax
from sklearn.preprocessing import StandardScaler
from scipy.special import logsumexp

onehot_encoder = OneHotEncoder(sparse_output=False)


def loss(X, Y, W):
    """
    Y: onehot encoded
    """
    Z = - X @ W
    N = X.shape[0]
    loss = 1/N * (np.trace(X @ W @ Y.T) + np.sum(logsumexp(Z, axis=1)))
    return loss

def gradient(X, Y, W, mu):
    """
    Y: onehot encoded 
    """
    Z = - X @ W
    P = softmax(Z, axis=1)
    N = X.shape[0]
    gd = 1/N * (X.T @ (Y - P)) + 2 * mu * W
    return gd

def gradient_descent(X, Y, max_iter=2000, eta=0.1, mu=0.01):
    """
    Very basic gradient descent algorithm with fixed eta and mu
    """
    Y_onehot = onehot_encoder.fit_transform(Y.reshape(-1,1))
    W = np.zeros((X.shape[1], Y_onehot.shape[1]))
    step = 0
    step_lst = [] 
    loss_lst = []
    W_lst = []
 
    while step < max_iter:
        step += 1
        W -= eta * gradient(X, Y_onehot, W, mu)
        step_lst.append(step)
        W_lst.append(W)
        loss_lst.append(loss(X, Y_onehot, W))

    df = pd.DataFrame({
        'step': step_lst, 
        'loss': loss_lst
    })
    return df, W

class Multiclass:
    def fit(self, X, Y):
        self.loss_steps, self.W = gradient_descent(X, Y)

    def loss_plot(self):
        return self.loss_steps.plot(
            x='step', 
            y='loss',
            xlabel='step',
            ylabel='loss'
        )

    def predict(self, H):
        Z = - H @ self.W
        P = softmax(Z, axis=1)
        return np.argmax(P, axis=1)

def run_ct_pb_mr(adata, sample_key, condition_key, n_splits, params, label_key, method, **kwargs):
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
        unique_conditions = np.unique(adata_.obs[condition_key])
        num_of_classes = len(unique_conditions)
        y = adata_[adata_.obs_names.isin(train)].obs[condition_key].map(rename_dict).to_numpy()
        print("Train shapes:")
        print(f"x.shape = {x.shape}")
        print(f"y.shape = {y.shape}")
        # val data
        x_val = pd.DataFrame(adata_[adata_.obs_names.isin(val)].X).to_numpy()
        y_val = adata_[adata_.obs_names.isin(val)].obs[condition_key].map(rename_dict).to_numpy()
        print("Val shapes:")
        print(f"x_val.shape = {x_val.shape}")
        print(f"y_val.shape = {y_val.shape}")
        # fit
        X = x
        Y = y
        model = Multiclass()
        model.fit(X, Y)
        print(f'Train accuracy = {np.sum(model.predict(X) == Y)/len(Y)}.')
        y_pred = model.predict(x_val)
        val_accuracy = np.sum(y_pred == y_val)/len(y_val)

        df = classification_report(y_val, y_pred, output_dict=True)
        df = pd.DataFrame(df).T
        df['split'] = i
        df['method'] = 'ct_pb_mr'
        dfs.append(df)
        
        val_accuracy = df["f1-score"]["accuracy"]

        val_accuracies.append(val_accuracy)
        val_avg.append(df["f1-score"]["weighted avg"])
        
        print(f'Val accuracy = {val_accuracy}.')
        print('===========================')

    df = pd.concat(dfs)

    print(f"Mean validation accuracy across 5 CV splits for a {method} model = {np.mean(np.array(val_accuracies))}.")
    print(f"Mean validation weighted avg across 5 CV splits for a {method} model = {np.mean(np.array(val_avg))}.")
    return df