from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

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

def run_pb_nn(adata, sample_key, condition_key, n_splits, params, hash, method, task, **kwargs):
    from utils import custom_pseudobulk_aggregation
    
    # Use custom pseudobulk aggregation to avoid decoupler compatibility issues
    adata_ = custom_pseudobulk_aggregation(adata, sample_key=sample_key, groups_col=None)

    rename_dict = {name: number for number, name in enumerate(np.unique(adata_.obs[condition_key]))}

    if params['norm'] is True:
        scaler = StandardScaler()
        adata_.X = scaler.fit_transform(adata_.X)
        
    adata_.obs[condition_key] = adata_.obs[condition_key].astype('category')

    num_features = adata_.shape[1]
    num_classes = len(adata_.obs[condition_key].cat.categories)
    batch_size = params['batch_size']
    learning_rate = params['lr']
    epochs = params['epochs']

    test_accuracies = []
    test_avg = []
    dfs = []

    for i in range(n_splits):
        print(f'Processing split = {i}...')
        df = adata.obs[[f'split{i}', sample_key]].drop_duplicates()
        train = list(df[df[f'split{i}'] == 'train'][sample_key])
        test = list(df[df[f'split{i}'] == 'val'][sample_key])
        # train data
        x = pd.DataFrame(adata_[adata_.obs_names.isin(train)].X).to_numpy()
        y = adata_[adata_.obs_names.isin(train)].obs[condition_key].cat.rename_categories(rename_dict)
        y = y.to_numpy()
        
        # Check for NaN/Inf values in data
        print(f"Data validation - x contains NaN: {np.isnan(x).any()}")
        print(f"Data validation - x contains Inf: {np.isinf(x).any()}")
        print(f"Data validation - x range: [{x.min():.3f}, {x.max():.3f}]")
        print(f"Data validation - x mean: {x.mean():.3f}, std: {x.std():.3f}")

        print("Train shapes:")
        print(f"x.shape = {x.shape}")
        print(f"y.shape = {y.shape}")
        # test data
        x_test = pd.DataFrame(adata_[adata_.obs_names.isin(test)].X).to_numpy()
        y_test = adata_[adata_.obs_names.isin(test)].obs[condition_key].cat.rename_categories(rename_dict)
        y_test = y_test.to_numpy()
        print("Test shapes:")
        print(f"x_test.shape = {x_test.shape}")
        print(f"y_test.shape = {y_test.shape}")
        
        # create datasets
        train_dataset = ClassifierDataset(torch.from_numpy(x).float(), torch.from_numpy(y).long())
        test_dataset = ClassifierDataset(torch.from_numpy(x_test).float(), torch.from_numpy(y_test).long())
        
        # create loaders
        train_loader = DataLoader(dataset=train_dataset, batch_size=batch_size)
        test_loader = DataLoader(dataset=test_dataset, batch_size=1)
        
        # init model
        model = MulticlassClassification(num_feature = num_features, num_class=num_classes)
        # define loss
        criterion = nn.CrossEntropyLoss()
        # define optimizer
        optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        
        # loss recorder
        train_losses = []
        train_accuracies = []
        
        # train
        print("Begin training.")
        for e in tqdm(range(1, epochs+1)):
            train_epoch_loss = 0
            train_epoch_acc = 0
            model.train()
            for X_train_batch, y_train_batch in train_loader:
                optimizer.zero_grad()
                y_train_pred = model(X_train_batch)
                train_loss = criterion(y_train_pred, y_train_batch)
                train_acc = multi_acc(y_train_pred, y_train_batch)
                train_loss.backward()
                optimizer.step()
                train_epoch_loss += train_loss.item()
                train_epoch_acc += train_acc.item()
            train_losses.append(train_epoch_loss/len(train_loader))
            train_accuracies.append(train_epoch_acc/len(train_loader))
            print(f'Epoch {e:03}: | Train Loss: {train_epoch_loss/len(train_loader):.5f} | Train Acc: {train_epoch_acc/len(train_loader):.3f}')
        
        # Plot training loss and accuracy
        fig_path = f'data/reports/{task}/{method}/{hash}/figures/'
        os.makedirs(fig_path, exist_ok = True)
        plt.figure(figsize=(10,4))
        plt.subplot(1,2,1)
        plt.plot(train_accuracies, label='Train Acc')
        plt.title('Train Accuracy/Epoch')
        plt.xlabel('Epoch')
        plt.ylabel('Accuracy')
        plt.legend()
        plt.subplot(1,2,2)
        plt.plot(train_losses, label='Train Loss')
        plt.title('Train Loss/Epoch')
        plt.xlabel('Epoch')
        plt.ylabel('Loss')
        plt.legend()
        plt.tight_layout()
        plt.savefig(fig_path + f'plot_loss_split_{i}.png')
        plt.close()
        
        # predict
        y_pred_list = []
        with torch.no_grad():
            model.eval()
            for X_batch, _ in test_loader:
                y_test_pred = model(X_batch)
                _, y_pred_tags = torch.max(y_test_pred, dim = 1)
                y_pred_list.append(y_pred_tags.cpu().numpy())
        y_pred_list = [a.squeeze().tolist() for a in y_pred_list]
        
        df = classification_report(y_test, y_pred_list, output_dict=True)
        df = pd.DataFrame(df).T
        df['split'] = i
        df['method'] = 'pb_nn'
        dfs.append(df)
        
        test_accuracies.append(df["f1-score"]["accuracy"])
        test_avg.append(df["f1-score"]["weighted avg"])
        
        print(f'Accuracy on the test set = {df["f1-score"]["accuracy"]}.')
        print('===========================')
        
    df = pd.concat(dfs)
    print(f"Mean test accuracy across {n_splits} CV splits for a NN model = {np.mean(np.array(test_accuracies))}.")
    print(f"Mean test weighted avg across {n_splits} CV splits for a NN model = {np.mean(np.array(test_avg))}.")
    return df