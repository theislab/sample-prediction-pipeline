import time
start_time = time.time()

import scanpy as sc
import pandas as pd
import multimil as mtm
from pathlib import Path
import shutil
import json

from matplotlib import pyplot as plt
import scanpy as sc
import pandas as pd
import anndata as ad
import os
import random
import torch
import scvi
from sklearn.metrics import classification_report
from utils import get_existing_checkpoints

import warnings
warnings.filterwarnings('ignore')


def run_multimil(adata, sample_key, condition_key, n_splits, params, hash, task, regression, **kwargs):
    # Store original data to avoid modifying it across splits
    original_adata = adata.copy()

    print('============ MultiMIL training ============')
    torch.set_float32_matmul_precision('medium')

    json_config = {}

    if regression is True:
        method = 'multimil_reg'
    else:
        method = 'multimil'
    donor = sample_key
    condition = condition_key
    n_splits = n_splits

    json_config['method'] = method
    json_config['donor'] = donor
    json_config['condition'] = condition
    json_config['n_splits'] = n_splits
    json_config['params'] = params
    json_config['hash'] = hash
    json_config['task'] = task
    json_config['regression'] = regression
    json_config['adata'] = str(adata)

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

    if regression is True:
        model_params = {
            "regression_loss_coef": params['regression_loss_coef'],
        }
    else:
        model_params = {
            "class_loss_coef": params['class_loss_coef'],
        }
    
    subset_umap = params['subset_umap']
    umap_colors = params['umap_colors'].strip('][').replace('\'', '').replace('\"', '').split(', ')
    train_params = {
        "max_epochs": params['train_max_epochs'],
        "save_checkpoint_every_n_epochs": params['train_save_checkpoint_every_n_epochs'],
    }

    # train params
    lr = params['lr']
    batch_size = params['batch_size']
    seed = params['seed']

    scvi.settings.seed = seed

    dfs = []

    for i in range(n_splits):
        print(f'Split {i}...')

        ########################
        ######## TRAIN #########
        ########################

        print('Organizing multimodal anndatas...')
        # Use fresh copy of original data for each split
            
        query = original_adata[original_adata.obs[f"split{i}"] == "val"].copy()
        adata = original_adata[original_adata.obs[f"split{i}"] == "train"].copy()

        idx = adata.obs[donor].sort_values().index
        adata = adata[idx].copy()

        print('Setting up anndata...')
        mtm.model.MILClassifier.setup_anndata(
            adata, 
            **setup_params
        )

        print('Initializing the model...')
        
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

        ###############################
        ###### TRAIN CLASSIFIER #######
        ###############################

        path_to_train_checkpoints = f'data/{method}/{task}/{hash}/{i}/checkpoints/'
        dirpath=Path(path_to_train_checkpoints)
        if dirpath.exists():
            shutil.rmtree(dirpath)

        os.makedirs(path_to_train_checkpoints, exist_ok=True)
        print('Starting training...')
        mil.train(
            lr=lr,
            batch_size=batch_size, 
            path_to_checkpoints=path_to_train_checkpoints,
            **train_params
        )

        mil.is_trained_ = True

        print('Starting inference...')
        mil.get_model_output(batch_size=batch_size)

        mil.save(f'data/{method}/{task}/{hash}/{i}/model/', overwrite=True)

        if subset_umap is not None:
            print(f'Subsetting to {subset_umap}...')
            actual_subset = min(subset_umap, len(adata.obs_names))
            idx = random.sample(list(adata.obs_names), actual_subset)
            adata = adata[idx].copy()

        if "X_umap" not in adata.obsm.keys():
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)

        sc.pl.umap(
            adata,
            color=umap_colors+["cell_attn"],
            ncols=1,
            show=False,
        )

        plt.savefig(f'data/{method}/{task}/{hash}/{i}/train_umap.png', bbox_inches="tight")
        plt.close()

        mil.plot_losses(save=f'data/{method}/{task}/{hash}/{i}/train_losses.png')

        checkpoints  = get_existing_checkpoints(path_to_train_checkpoints)

        print(f"Found {len(checkpoints)} in {path_to_train_checkpoints}.")

        for ckpt in checkpoints:

            train_state_dict = torch.load(path_to_train_checkpoints + f'{ckpt}.ckpt')['state_dict']
            for key in list(train_state_dict.keys()):
                train_state_dict[key.replace('module.', '')] = train_state_dict.pop(key)

            mil.module.load_state_dict(train_state_dict)

            idx = query.obs[donor].sort_values().index
            query = query[idx].copy()

            new_model = mtm.model.MILClassifier.load_query_data(query, reference_model=mil)

            new_model.is_trained_ = True
            new_model.get_model_output(query, batch_size=batch_size)

            report = classification_report(
                query.obs[condition], query.obs[f"predicted_{condition}"], output_dict=True
            )
            df = pd.DataFrame(report).T
            df.to_csv(path_to_train_checkpoints + f'{ckpt}.csv')

            df['split'] = i
            df['method'] = method
            df['epoch'] = ckpt.split('-')[0].split('=')[-1]
            
            dfs.append(df)

            adata.obs['reference'] = 'reference'
            query.obs['reference'] = 'query'
            adata_both = ad.concat([adata, query])

            if subset_umap is not None:
                print(f'Subsetting to {subset_umap}...')
                actual_subset = min(subset_umap, len(adata_both.obs_names))
                idx = random.sample(list(adata_both.obs_names), actual_subset)
                adata_both = adata_both[idx].copy()

            if "X_umap" not in adata_both.obsm.keys():
                sc.pp.neighbors(adata_both)
                sc.tl.umap(adata_both)

            sc.pl.umap(
                adata_both,
                color=umap_colors+["cell_attn", "reference"],
                ncols=1,
                show=False,
            )
            plt.savefig(path_to_train_checkpoints + f'{ckpt}_umap.png', bbox_inches="tight")
            plt.close()

    with open(f'data/{method}/{task}/{hash}/config.json', 'w') as f:
        json.dump(json_config, f)

    df = pd.concat(dfs)
    return df
