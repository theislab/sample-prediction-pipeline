import pandas as pd
import yaml
import hashlib
import os
import numpy as np
import anndata as ad

# from https://github.com/HCA-integration/hca_integration_toolbox/blob/main/workflow/utils/misc.py#L129
def create_hash(string: str, digest_size: int = 5):
    string = string.encode('utf-8')
    return hashlib.blake2b(string, digest_size=digest_size).hexdigest()

def create_tasks_df(config, save=None):
    tasks_df = []
    with open(config, "r") as stream:
        params = yaml.safe_load(stream)
    for task in params['TASKS']:
        task_dict = params['TASKS'][task]
        params_list = []
        method_dfs = []
        for method in task_dict['methods']:
            method_params = task_dict['methods'][method]
            # Read parameter file - don't use index_col=0 for multi-column files
            df_params = pd.read_csv(method_params, sep='\t')
            params_list = [str(row) for row in df_params.to_dict(orient='records')]
            method_df = {}
            method_df['params'] = params_list
            method_df['hash'] = [create_hash(row + method + task) for row in params_list]
            method_df['method'] = method
            method_dfs.append(pd.DataFrame(method_df))
        method_dfs = pd.concat(method_dfs)
        method_dfs['task'] = task
        for key in task_dict:
            if key != 'methods':
                method_dfs[key] = task_dict[key] 
        tasks_df.append(method_dfs)
    tasks_df = pd.concat(tasks_df)
    if save is not None:
        tasks_df.to_csv(save, sep='\t')
    return tasks_df

def get_existing_checkpoints(rootdir):

    checkpoints = []

    for _, _, files in os.walk(rootdir):
        for filename in files:
            if filename.endswith('.ckpt'):
                checkpoints.append(filename.strip('.ckpt'))

    return checkpoints

def safe_standardize(X, axis=0, eps=1e-8):
    """
    Safely standardize data to prevent NaN values from division by zero.
    
    Args:
        X: Input array to standardize
        axis: Axis along which to compute mean and std
        eps: Small value to prevent division by zero
    
    Returns:
        Standardized array with NaN values replaced
    """
    if X.size == 0:
        return X
    
    # Compute mean and std
    mean_val = np.mean(X, axis=axis, keepdims=True)
    std_val = np.std(X, axis=axis, keepdims=True)
    
    # Add small epsilon to prevent division by zero
    std_val = np.maximum(std_val, eps)
    
    # Standardize
    X_std = (X - mean_val) / std_val
    
    # Replace any remaining NaN values with 0
    X_std = np.nan_to_num(X_std, nan=0.0, posinf=0.0, neginf=0.0)
    
    return X_std

def validate_and_clean_data(X, y, min_samples=2, min_features=1, verbose=True):
    """
    Validate and clean input data to prevent NaN issues.
    
    Args:
        X: Feature matrix
        y: Target labels
        min_samples: Minimum number of samples required
        min_features: Minimum number of features required
        verbose: Whether to print diagnostic information
    
    Returns:
        tuple: (cleaned_X, cleaned_y, is_valid)
    """
    if verbose:
        print(f"Input data shape: X={X.shape}, y={y.shape}")
        print(f"X has NaN: {np.isnan(X).any()}")
        print(f"y has NaN: {np.isnan(y).any()}")
    
    # Check for sufficient samples
    if len(y) < min_samples:
        if verbose:
            print(f"Insufficient samples: {len(y)} < {min_samples}")
        return X, y, False
    
    # Check for sufficient features
    if X.shape[1] < min_features:
        if verbose:
            print(f"Insufficient features: {X.shape[1]} < {min_features}")
        return X, y, False
    
    # Check for NaN values
    if np.isnan(X).any():
        if verbose:
            print("Found NaN values in X, attempting to clean...")
        # Replace NaN with 0
        X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    
    if np.isnan(y).any():
        if verbose:
            print("Found NaN values in y, removing corresponding samples...")
        # Remove samples with NaN labels
        valid_mask = ~np.isnan(y)
        X = X[valid_mask]
        y = y[valid_mask]
        
        if len(y) < min_samples:
            if verbose:
                print(f"After removing NaN labels, insufficient samples: {len(y)} < {min_samples}")
            return X, y, False
    
    # Check for zero variance features
    feature_vars = np.var(X, axis=0)
    zero_var_features = np.where(feature_vars == 0)[0]
    
    if len(zero_var_features) > 0:
        if verbose:
            print(f"Found {len(zero_var_features)} zero-variance features, removing...")
        # Remove zero variance features
        X = np.delete(X, zero_var_features, axis=1)
        
        if X.shape[1] < min_features:
            if verbose:
                print(f"After removing zero-variance features, insufficient features: {X.shape[1]} < {min_features}")
            return X, y, False
    
    if verbose:
        print(f"Cleaned data shape: X={X.shape}, y={y.shape}")
        print(f"Final X has NaN: {np.isnan(X).any()}")
        print(f"Final y has NaN: {np.isnan(y).any()}")
    
    return X, y, True

def custom_pseudobulk_aggregation(adata, sample_key, groups_col=None):
    """
    Custom pseudobulk aggregation function to replace decoupler's get_pseudobulk.
    
    Args:
        adata: AnnData object
        sample_key: Column name for patient/sample grouping
        groups_col: Column name for additional grouping (e.g., cell types), optional
    
    Returns:
        AnnData object with aggregated pseudobulk data
    """
    print("Performing custom pseudobulk aggregation...")
    
    # Ensure sample_key is string type
    adata.obs[sample_key] = adata.obs[sample_key].astype(str)
    
    if groups_col is not None:
        # Group by both sample and cell type
        adata.obs[groups_col] = adata.obs[groups_col].astype(str)
        groups = adata.obs.groupby([sample_key, groups_col])
        print(f"Found {len(groups)} sample-celltype combinations")
        
        pseudobulk_data = []
        pseudobulk_obs = []
        
        for (patient, celltype), group_indices in groups.groups.items():
            # Get the subset of adata for this patient-celltype combination
            group_adata = adata[group_indices]
            
            # Average expression across cells for this patient-celltype combination
            # Convert to numpy array if it's a DataFrame
            if hasattr(group_adata.X, 'to_numpy'):
                group_X = group_adata.X.to_numpy()
            else:
                group_X = group_adata.X
                
            # Handle NaN values in the data
            if np.isnan(group_X).any():
                print(f"Warning: Found NaN values in group {patient}_{celltype}, replacing with 0")
                group_X = np.nan_to_num(group_X, nan=0.0, posinf=0.0, neginf=0.0)
            
            avg_expression = group_X.mean(axis=0)
            
            # Check for NaN in averaged expression
            if np.isnan(avg_expression).any():
                print(f"Warning: Found NaN in averaged expression for {patient}_{celltype}, replacing with 0")
                avg_expression = np.nan_to_num(avg_expression, nan=0.0, posinf=0.0, neginf=0.0)
            
            # Get the first occurrence of metadata for this group
            group_obs = group_adata.obs.iloc[0].copy()
            group_obs.name = f"{patient}_{celltype}"  # Set observation name
            
            pseudobulk_data.append(avg_expression)
            pseudobulk_obs.append(group_obs)
    else:
        # Group by patient only
        patients = adata.obs[sample_key].unique()
        print(f"Found {len(patients)} patients")
        
        pseudobulk_data = []
        pseudobulk_obs = []
        
        for patient in patients:
            patient_mask = adata.obs[sample_key] == patient
            patient_data = adata[patient_mask]
            
            # Average expression across cells for this patient
            # Convert to numpy array if it's a DataFrame
            if hasattr(patient_data.X, 'to_numpy'):
                patient_X = patient_data.X.to_numpy()
            else:
                patient_X = patient_data.X
                
            avg_expression = patient_X.mean(axis=0)

            # Get the first occurrence of metadata for this patient
            patient_obs = patient_data.obs.iloc[0].copy()
            patient_obs.name = patient  # Set observation name to patient ID
            
            pseudobulk_data.append(avg_expression)
            pseudobulk_obs.append(patient_obs)
    
    # Create new AnnData object
    adata_ = ad.AnnData(
        X=np.vstack(pseudobulk_data),
        obs=pd.DataFrame(pseudobulk_obs),
        var=adata.var.copy()
    )
    
    print(f"Pseudobulk shape: {adata_.shape}")
    print(f"Pseudobulk observations: {adata_.obs_names.tolist()}")
    
    return adata_
                