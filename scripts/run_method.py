import time 
start_time = time.time()
import scanpy as sc
import ast

from pb_rf import run_pb_rf
from gex_rf import run_gex_rf
from pb_nn import run_pb_nn
from gex_nn import run_gex_nn
from pb_mr import run_pb_mr
from gex_mr import run_gex_mr
from ct_pb_nn import run_ct_pb_nn
from ct_pb_rf import run_ct_pb_rf
from ct_pb_mr import run_ct_pb_mr
from freq_mr import run_freq_mr
from freq_rf import run_freq_rf
from freq_nn import run_freq_nn
from run_multimil import run_multimil

print('--- %s seconds ---' % (time.time() - start_time))

METHOD_MAP = dict(
    pb_rf=dict(function=run_pb_rf),
    gex_rf=dict(function=run_gex_rf),
    pb_nn=dict(function=run_pb_nn),
    gex_nn=dict(function=run_gex_nn),
    pb_mr=dict(function=run_pb_mr),
    gex_mr=dict(function=run_gex_mr),
    ct_pb_nn=dict(function=run_ct_pb_nn),
    ct_pb_rf=dict(function=run_ct_pb_rf),
    ct_pb_mr=dict(function=run_ct_pb_mr),
    freq_mr=dict(function=run_freq_mr),
    freq_rf=dict(function=run_freq_rf),
    freq_nn=dict(function=run_freq_nn),
    multimil=dict(function=run_multimil),
    multimil_reg=dict(function=run_multimil),
)

params = snakemake.params.params

method_params = ast.literal_eval(params['params']) # this is dict
input = params['input']
label_key = params['label_key']
condition_key = params['condition_key']
condition_regression_key = params.get('condition_regression_key', None)
sample_key = params['sample_key']
n_splits = params['n_splits']
h = params['hash']
method = params['method']
task = params['task']
output_file = snakemake.output.tsv

method_function = METHOD_MAP[method]['function']

regression = False
if (method == 'multimil_reg'):
    regression = True

adata = sc.read_h5ad(input)

# Handle different function signatures based on method requirements
if method in ['multimil', 'multimil_reg']:
    df = method_function(
        adata=adata, 
        sample_key=sample_key, 
        condition_key=condition_key, 
        n_splits=n_splits, 
        params=method_params,
        hash=h,
        method=method,
        task=task,
        regression=regression,
    )
elif method == 'gex_rf':
    df = method_function(
        adata=adata, 
        condition_key=condition_key, 
        n_splits=n_splits, 
        params=method_params,
    )
elif method in ['ct_pb_rf', 'ct_pb_mr']:
    df = method_function(
        adata=adata, 
        sample_key=sample_key, 
        condition_key=condition_key, 
        n_splits=n_splits, 
        params=method_params,
        label_key=label_key,
        method=method,
    )
elif method in ['freq_rf', 'freq_mr']:
    df = method_function(
        adata=adata, 
        sample_key=sample_key, 
        condition_key=condition_key, 
        n_splits=n_splits, 
        params=method_params,
        label_key=label_key,
        method=method,
    )
elif method in ['freq_nn']:
    df = method_function(
        adata=adata, 
        sample_key=sample_key, 
        condition_key=condition_key, 
        n_splits=n_splits, 
        params=method_params,
        label_key=label_key,
        method=method,
        task=task,
    )
elif method in ['pb_mr', 'gex_mr']:
    df = method_function(
        adata=adata, 
        sample_key=sample_key, 
        condition_key=condition_key, 
        n_splits=n_splits, 
        params=method_params,
        method=method,
    )
elif method in ['pb_nn', 'gex_nn']:
    df = method_function(
        adata=adata, 
        sample_key=sample_key, 
        condition_key=condition_key, 
        n_splits=n_splits, 
        params=method_params,
        hash=h,
        method=method,
        task=task,
    )
elif method in ['ct_pb_nn']:
    df = method_function(
        adata=adata, 
        sample_key=sample_key, 
        condition_key=condition_key, 
        n_splits=n_splits, 
        params=method_params,
        hash=h,
        method=method,
        task=task,
        label_key=label_key,
    )
else:
    # Default case for pb_rf and other simple methods
    df = method_function(
        adata=adata, 
        sample_key=sample_key, 
        condition_key=condition_key, 
        n_splits=n_splits, 
        params=method_params,
    )

df['hash'] = h
df['method_params'] = params['params']
df['task'] = task
df.to_csv(output_file, sep='\t')
