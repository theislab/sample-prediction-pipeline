import pandas as pd
from snakemake.utils import Paramspace
from scripts.utils import create_tasks_df
from pprint import pprint
import numpy as np
import os

configfile: 'config.yaml'
tasks_df = create_tasks_df('config.yaml', save='data/tasks.tsv')
hashes = np.unique(tasks_df['hash'])
matched_methods = [tasks_df.loc[tasks_df['hash'] == h, 'method'].values[0] for h in hashes]
matched_tasks = [tasks_df.loc[tasks_df['hash'] == h, 'task'].values[0] for h in hashes]
tasks = np.unique(tasks_df['task'])
multimil_tasks_df = tasks_df.loc[(tasks_df['method'] == 'multimil') | (tasks_df['method'] == 'multimil_reg'), :]
multimil_tasks_df = multimil_tasks_df[['method', 'task']].drop_duplicates()

rule run_method:
    output:
        tsv='data/reports/{task}/{method}/{hash}/accuracy.tsv'
    params:
        params=lambda wildcards: tasks_df.loc[tasks_df['hash'] == wildcards.hash, :].iloc[0, :].to_dict()
    conda:
        'envs/sample_prediction_pipeline.yaml'
    script:
        'scripts/run_method.py'

rule merge:
    input:
        expand(
            rules.run_method.output.tsv,
            zip,
            task=matched_tasks,
            method=matched_methods,
            hash=hashes,
        )
    output:
        tsv='data/reports/methods.tsv'
    run:
        input_files = [input] if isinstance(input, str) else input
        dfs = [pd.read_table(file) for file in input_files if os.path.exists(file)]
        metrics_df = pd.concat(dfs)
        metrics_df.to_csv(output.tsv, sep='\t', index=False)

rule find_best_per_task:
    input:
        tsv='data/reports/methods.tsv'
    output:
        tsv='data/reports/best.tsv'
    params:
        config=config
    conda:
    	'envs/sample_prediction_pipeline.yaml'
    script:
        'scripts/find_best.py'

rule plot_accuracy:
    input:
        tsv='data/reports/best.tsv'
    output:
        expand('data/reports/{task}_accuracy.png', task=tasks)
    conda:
    	'envs/sample_prediction_pipeline.yaml'
    script:
        'scripts/plot_accuracy.py'

rule save_anndata:
    input:
        tsv='data/reports/best.tsv'
    output:
        txt='data/reports/{task}/best_{method}.txt'
    conda:
    	'envs/sample_prediction_pipeline.yaml'
    params:
        config=config
    script:
        'scripts/save_anndata.py'  

rule all:
    input:
        expand('data/reports/{task}_accuracy.png', task=tasks),
        expand('data/reports/{task}/best_{method}.txt', zip, task=multimil_tasks_df['task'].values, method = multimil_tasks_df['method'].values)
    default_target: True
