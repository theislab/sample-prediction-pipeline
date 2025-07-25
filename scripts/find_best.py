import pandas as pd
import numpy as np


input_file = snakemake.input.tsv
config = snakemake.params['config']

df = pd.read_csv(input_file, sep='\t', index_col=0)

best_methods = {}

for task in config['TASKS']:
    n_splits = config['TASKS'][task]['n_splits']
    df_task = df[df['task'] == task]
    best_methods_per_task = {}

    for method in np.unique(df_task['method']):
        df_method = df_task[df_task['method'] == method]
        best_accuracy = -1
        best_hash = None
        best_accuracies = None

        for h in np.unique(df_method['hash']):
            df_tmp = df_method[df_method['hash'] == h]
            # for multimil need to check if all epochs and query_epochs are present and group by those
            # check only if all splits are present
            if 'multimil' in method:
                if 'query_epoch' in df_tmp.columns:
                    if df_tmp.drop(['epoch', 'query_epoch'], axis=1).isnull().values.any():
                        continue
                else:
                    if df_tmp.drop(['epoch'], axis=1).isnull().values.any():
                        continue
                for epoch in np.unique(df_tmp['epoch']):
                    df_tmp_epoch = df_tmp[df_tmp['epoch'] == epoch]
                    if len(np.unique(df_tmp_epoch['split'])) != n_splits:
                        continue
                    accuracies = []
                    for i in range(n_splits):
                        accuracies.append(df_tmp_epoch[df_tmp_epoch['split'] == i]['f1-score']['accuracy'])
                    accuracy = np.mean(accuracies)
                    if accuracy > best_accuracy:
                        best_accuracy = accuracy
                        best_hash = h
                        best_params = df_tmp_epoch['method_params'].iloc[0]
                        best_epoch = epoch
                        best_query_epoch = np.nan
                        best_accuracies = accuracies
            else:
                if len(np.unique(df_tmp['split'])) != n_splits:
                    continue
                accuracies = []
                for i in range(n_splits):
                    accuracies.append(df_tmp[df_tmp['split'] == i]['f1-score']['accuracy'])
                accuracy = np.mean(accuracies)
                if accuracy > best_accuracy:
                    best_accuracy = accuracy
                    best_hash = h
                    best_params = df_tmp['method_params'].iloc[0]
                    best_epoch = np.nan
                    best_query_epoch = np.nan
                    best_accuracies = accuracies

        if best_hash is not None: # this can happen if e.g. all multimil jobs failed because of a too high lr
            best_methods_per_task[method] = {
                'hash': best_hash,
                'accuracy': best_accuracy,
                'best_params': best_params,
                'best_epoch': best_epoch,
                'best_query_epoch': best_query_epoch,
                'accuracies': best_accuracies,
                'method': method,
            }
    best_df = pd.DataFrame.from_dict(best_methods_per_task).T
    best_df['task'] = task
    best_methods[task] = best_df
    best_df.to_csv(f'data/reports/{task}/best.csv', sep='\t', index=False)

best_all = pd.concat(best_methods)
best_all.to_csv(snakemake.output.tsv, sep='\t', index=False)