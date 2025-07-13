import scanpy as sc
import pandas as pd
import numpy as np
import multimil as mtm
from pathlib import Path
import ast
import os.path

import scanpy as sc
import pandas as pd
import anndata as ad
import torch
import scvi
from pprint import pprint
from utils import get_existing_checkpoints

import warnings
warnings.filterwarnings('ignore')


input_file = snakemake.input.tsv
output_file = snakemake.output.txt
config = snakemake.params['config']

df = pd.read_csv(input_file, sep='\t', index_col=None)

task = snakemake.wildcards.task
method = snakemake.wildcards.method
df_best = df[(df['task'] == task) & (df['method'] == method)]

if len(df_best) == 0:
    with open(output_file, 'w') as f:
        f.write(f'No anndatas to save for {method} for {task}!')
    exit(0)

assert(len(df_best) == 1)

input = config['TASKS'][task]['input']

condition_key = config['TASKS'][task]['condition_key']
sample_key = config['TASKS'][task]['sample_key']
n_splits = config['TASKS'][task]['n_splits']

donor = sample_key
condition = condition_key
n_splits = n_splits

row = df_best.iloc[0]
h = row['hash']
train_epoch = row['best_epoch']
query_epoch = row['best_query_epoch']
params = ast.literal_eval(row['best_params'])

print(f'Processing {task} and {method}...')

path = f'data/{method}/{task}/{h}_adata_both.h5ad'
if os.path.isfile(path):
    print('AnnData already exists. Exiting...')
    with open(output_file, 'w') as f:
        f.write(f'Saved anndata for {method} for {task}!')
    exit(0)

if 'multimil' in method:
    
    regression = False
    if 'reg' in method:
        regression = True
    
    original_adata = sc.read_h5ad(input)

    setup_params = {}
    cat_cov = params.get('categorical_covariate_keys', None)
    if cat_cov is None or cat_cov.strip('[]').replace("'", '').replace('"', '').strip() == '':
        setup_params["categorical_covariate_keys"] = [sample_key, condition_key]
    else:
        setup_params["categorical_covariate_keys"] = cat_cov.strip('][').replace("'", '').replace('"', '').split(', ')
        if sample_key not in setup_params["categorical_covariate_keys"]:
            setup_params["categorical_covariate_keys"].append(sample_key)
        if condition_key not in setup_params["categorical_covariate_keys"]:
            setup_params["categorical_covariate_keys"].append(condition_key)

    model_params = {}
    if regression is True:
        model_params["regression_loss_coef"] = params['regression_loss_coef']
    else:
        model_params["class_loss_coef"] = params['class_loss_coef']

    lr = params['lr']
    batch_size = params['batch_size']
    seed = params['seed']

    scvi.settings.seed = seed

    for i in range(n_splits):
        # Start from a fresh copy of the original AnnData for each split
        adata = original_adata.copy()
            
        query = adata[adata.obs[f"split{i}"] == "val"].copy()
        adata = adata[adata.obs[f"split{i}"] == "train"].copy()

        idx = adata.obs[donor].sort_values().index
        adata = adata[idx].copy()

        mtm.model.MILClassifier.setup_anndata(
            adata, 
            **setup_params
        )

        print('Initializing the model...')
        
        path_to_train_checkpoints = f'data/{method}/{task}/{h}/{i}/checkpoints/'
        
        # First, let's check what type of model was saved by examining the checkpoint
        train_checkpoints = get_existing_checkpoints(path_to_train_checkpoints)
        best_ckpt = None
        for ckpt in train_checkpoints:
            if str(int(train_epoch)) in ckpt:
                best_ckpt = ckpt
                break
        
        # Load the checkpoint to determine the model type
        checkpoint = torch.load(path_to_train_checkpoints + f'{best_ckpt}.ckpt')
        state_dict = checkpoint['state_dict']
        
        # # Check if it's a classification or regression model based on state dict keys
        # has_classifiers = any('classifiers' in key for key in state_dict.keys())
        # has_regressors = any('regressors' in key for key in state_dict.keys())
        
        # if has_classifiers:
        #     # Check if the number of classes matches between saved model and current data
        #     saved_num_classes = None
        #     for key in state_dict.keys():
        #         if 'classifiers.0.1.weight' in key:
        #             saved_num_classes = state_dict[key].shape[0]
        #             break
            
        #     current_num_classes = len(adata.obs[condition].cat.categories)
            
        #     if saved_num_classes != current_num_classes:
        #         continue
            
        #     # It's a classification model
        #     mil = mtm.model.MILClassifier(
        #         adata, 
        #         classification=[
        #             condition
        #         ],
        #         sample_key=donor,
        #         **model_params,
        #     )
        # elif has_regressors:
        #     # It's a regression model
        #     mil = mtm.model.MILClassifier(
        #         adata, 
        #         ordinal_regression=[
        #             condition
        #         ],
        #         sample_key=donor,
        #         **model_params,
        #     )
        # else:
        #     # Fallback to the original logic
        if regression is True:
            mil = mtm.model.MILClassifier(
                adata, 
                ordinal_regression=[
                    condition
                ],
                sample_key=donor,
                **model_params,
            )
        else:
            mil = mtm.model.MILClassifier(
                adata, 
                classification=[
                    condition
                ],
                sample_key=donor,
                **model_params,
            )

        path_to_train_checkpoints = f'data/{method}/{task}/{h}/{i}/checkpoints/'
        
        # Use the already loaded state_dict from earlier
        train_state_dict = state_dict
        for key in list(train_state_dict.keys()):
            train_state_dict[key.replace('module.', '')] = train_state_dict.pop(key)

        # Load the compatible parts of the state dict
        mil.module.load_state_dict(train_state_dict, strict=False)
        
        # Reinitialize the model to ensure all components are properly set up
        mil.is_trained_ = True
        
        # mil.get_model_output(adata, batch_size=batch_size)
        
        idx = query.obs[donor].sort_values().index
        query = query[idx].copy()

        new_model = mtm.model.MILClassifier.load_query_data(query, reference_model=mil)

        new_model.is_trained_ = True
        new_model.get_model_output(query, batch_size=batch_size)
        
        adata_both = ad.concat([adata, query])
        
        # Only save cell attention if it exists
        if 'cell_attn' in adata_both.obs.columns:
            original_adata.obs[f'cell_attn_{i}'] = adata_both.obs['cell_attn']
        else:
            print(f"No cell_attn found in adata_both for split {i}")
        
        print(f"Completed split {i}")
    
    if 'cell_attn' in adata_both.obs.columns:
        cell_attn_cols = [f'cell_attn_{i}' for i in range(n_splits)]
        original_adata.obs['cell_attn'] = np.mean([original_adata.obs[col] for col in cell_attn_cols], axis=0)
        
    original_adata.write(f'data/{method}/{task}/{h}_adata_both.h5ad')
else:
    raise ValueError(f'Unknown method: {method}')

with open(output_file, 'w') as f:
    f.write(f'Saved anndata for {method} for {task}!')