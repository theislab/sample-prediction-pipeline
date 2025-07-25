import seaborn as sns
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import os
from tqdm import tqdm
import scipy
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

from sklearn.preprocessing import StandardScaler    
from sklearn.metrics import classification_report


class ClassifierDataset(Dataset):
    
    def __init__(self, X_data, y_data):
        self.X_data = X_data
        self.y_data = y_data
        
    def __getitem__(self, index):
        return self.X_data[index], self.y_data[index]
        
    def __len__ (self):
        return len(self.X_data)

class MulticlassClassification(nn.Module):
    def __init__(self, num_feature, num_class):
        super(MulticlassClassification, self).__init__()
        
        self.layer_1 = nn.Linear(num_feature, 64)
        self.layer_out = nn.Linear(64, num_class) 
        
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(p=0.2)
        self.batchnorm1 = nn.BatchNorm1d(64)
        
    def forward(self, x):
        x = self.layer_1(x)
        if x.shape[0] > 1:
            x = self.batchnorm1(x)
        x = self.relu(x)
        x = self.dropout(x)
        x = self.layer_out(x)
        
        return x
    
def multi_acc(y_pred, y_test):
    y_pred_softmax = torch.log_softmax(y_pred, dim = 1)
    _, y_pred_tags = torch.max(y_pred_softmax, dim = 1)    
    
    correct_pred = (y_pred_tags == y_test).float()
    acc = correct_pred.sum() / len(correct_pred)
    
    acc = torch.round(acc * 100)
    
    return acc

def run_gex_nn(adata, sample_key, condition_key, n_splits, params, hash, method, task, **kwargs):

    adata.obs[condition_key] = adata.obs[condition_key].astype('category')
    rename_dict = {name: number for number, name in enumerate(np.unique(adata.obs[condition_key]))}
    
    if params['norm'] is True:
        scaler = StandardScaler()
        adata.X = scaler.fit_transform(adata.X)

    num_features = adata.shape[1]
    num_classes = len(adata.obs[condition_key].cat.categories)
    batch_size = params['batch_size']
    learning_rate = params['lr']
    epochs = params['epochs']

    val_accuracies = []
    val_avg = []
    dfs = []

    for i in range(n_splits):
        print(f"Split {i}...")
        df = adata.obs[[f'split{i}', sample_key]].drop_duplicates()
        train = list(df[df[f'split{i}'] == 'train'][sample_key])
        test = list(df[df[f'split{i}'] == 'val'][sample_key])
        if scipy.sparse.issparse(adata.X):
            x = pd.DataFrame(adata[adata.obs[sample_key].isin(train)].X.A).to_numpy()
        else:
            x = pd.DataFrame(adata[adata.obs[sample_key].isin(train)].X).to_numpy()
        y = adata[adata.obs[sample_key].isin(train)].obs[condition_key].cat.rename_categories(rename_dict)
        y = y.to_numpy()
        print(f'Train shape = ({x.shape}, {y.shape}).')

        if scipy.sparse.issparse(adata.X):
            x_test = pd.DataFrame(adata[adata.obs[sample_key].isin(test)].X.A).to_numpy()
        else:
            x_test = pd.DataFrame(adata[adata.obs[sample_key].isin(test)].X).to_numpy()
        y_test = adata[adata.obs[sample_key].isin(test)].obs[condition_key].cat.rename_categories(rename_dict)
        y_test = y_test.to_numpy()
        print(f'Test shape = ({x_test.shape}, {y_test.shape}).')

        # Use all training data
        X_train = x
        y_train = y
        X_test = x_test
        y_test = y_test

        train_dataset = ClassifierDataset(torch.from_numpy(X_train).float(), torch.from_numpy(y_train).long())
        test_dataset = ClassifierDataset(torch.from_numpy(X_test).float(), torch.from_numpy(y_test).long())

        train_loader = DataLoader(dataset=train_dataset, batch_size=batch_size)
        test_loader = DataLoader(dataset=test_dataset, batch_size=1)

        model = MulticlassClassification(num_feature = num_features, num_class=num_classes)
        
        criterion = nn.CrossEntropyLoss()
        optimizer = optim.Adam(model.parameters(), lr=learning_rate)

        accuracy_stats = {
            'train': []
        }
        loss_stats = {
            'train': []
        }

        print("Begin training.")
        for e in tqdm(range(1, epochs+1)):
            # TRAINING
            train_epoch_loss = 0
            train_epoch_acc = 0
            model.train()
            for X_train_batch, y_train_batch in train_loader:
                X_train_batch, y_train_batch = X_train_batch, y_train_batch
                optimizer.zero_grad()
        
                y_train_pred = model(X_train_batch)
        
                train_loss = criterion(y_train_pred, y_train_batch)
                train_acc = multi_acc(y_train_pred, y_train_batch)
        
                train_loss.backward()
                optimizer.step()
        
                train_epoch_loss += train_loss.item()
                train_epoch_acc += train_acc.item()

            loss_stats['train'].append(train_epoch_loss/len(train_loader))
            accuracy_stats['train'].append(train_epoch_acc/len(train_loader))

        print(f'Epoch {e+0:03}: | Train Loss: {train_epoch_loss/len(train_loader):.5f} | Train Acc: {train_epoch_acc/len(train_loader):.3f}')
        
        # Create dataframes
        train_acc_df = pd.DataFrame.from_dict(accuracy_stats).reset_index().melt(id_vars=['index']).rename(columns={"index":"epochs"})
        train_loss_df = pd.DataFrame.from_dict(loss_stats).reset_index().melt(id_vars=['index']).rename(columns={"index":"epochs"})
        # Plot the dataframes
        _, axes = plt.subplots(nrows=1, ncols=2, figsize=(20,7))

        fig_path = f'data/reports/{task}/{method}/{hash}/figures/'
        os.makedirs(fig_path, exist_ok = True)
        sns.lineplot(data=train_acc_df, x = "epochs", y="value", hue="variable",  ax=axes[0]).set_title('Train Accuracy/Epoch')
        sns.lineplot(data=train_loss_df, x = "epochs", y="value", hue="variable", ax=axes[1]).set_title('Train Loss/Epoch').get_figure().savefig(fig_path + f'plot_loss_split_{i}.png')
        
        # predict
        y_pred_list = []
        with torch.no_grad():
            model.eval()
            for X_batch, _ in test_loader:
                X_batch = X_batch
                y_test_pred = model(X_batch)
                _, y_pred_tags = torch.max(y_test_pred, dim = 1)
                y_pred_list.append(y_pred_tags.cpu().numpy())
        y_pred_list = [a.squeeze().tolist() for a in y_pred_list]
        
        df = classification_report(y_test, y_pred_list, output_dict=True)
        df = pd.DataFrame(df).T
        df['split'] = i
        df['method'] = 'gex_nn'
        dfs.append(df)
        
        val_accuracy = df["f1-score"]["accuracy"]

        val_accuracies.append(val_accuracy)
        val_avg.append(df["f1-score"]["weighted avg"])
        
        print(f'Test accuracy = {val_accuracy}.')
        print('===========================')
        
    df = pd.concat(dfs)
    
    print(f"Mean test accuracy across {n_splits} CV splits for a {method} model = {np.mean(np.array(val_accuracies))}.")
    print(f"Mean test weighted avg across {n_splits} CV splits for a {method} model = {np.mean(np.array(val_avg))}.")
    return df