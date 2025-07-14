# Sample Prediction Pipeline Guide

This guide covers installation, setup, and usage of the Sample Prediction Pipeline for sample-level prediction benchmarking.

## Prerequisites

Before installing the pipeline, ensure you have:

- **Python 3.10 or higher**
- **Mamba** (recommended, faster than conda) or **Conda/Miniconda**
- **Git** for cloning the repository

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/theislab/sample-prediction-pipeline.git
cd sample-prediction-pipeline
```

### 2. Create Mamba Environment

The pipeline uses a mamba environment to manage dependencies:

```bash
mamba env create -f envs/sample_prediction_pipeline.yaml
```

This installs all required packages including:
- Python 3.10
- Scanpy
- Pandas
- Matplotlib
- Snakemake
- MultiMIL
- Other dependencies

### 3. Activate the Environment

```bash
mamba activate sample_prediction_pipeline
```

## Quick Start

The repository includes ready-to-use example data and configuration files.

### 1. Use the Provided Example Data

1. **Run the pipeline:**
   ```bash
   snakemake --cores 1
   ```

This will run all methods on the example dataset and produce results in `data/reports/`.

## Configuration

### Task Configuration

Create or modify `config.yaml` to define your analysis:

```yaml
TASKS:
  your_task_name:
    input: path/to/your_data.h5ad
    label_key: cell_type  # column name for target variable
    condition_key: condition  # column name for condition
    sample_key: sample_id     # column name for sample IDs
    n_splits: 5  # number of cross-validation splits
    methods:
      multimil: params/your_task/multimil.tsv
      pb_rf: params/your_task/pb_rf.tsv
      gex_nn: params/your_task/gex_nn.tsv
```

### Parameter Files

For each method, create a TSV file with parameter combinations (see `params/example_task/` for examples).

## Data Requirements

Your input AnnData (.h5ad) file **must** contain:

1. **Expression Matrix** (`adata.X`): Normalized gene expression (cells × genes)
2. **Cell Metadata** (`adata.obs`): DataFrame with columns for:
   - Sample IDs
   - Cell types/labels
   - Condition information
   - **Split columns**: For each cross-validation split, you must have columns named `split0`, `split1`, ..., `split{n_splits-1}` in `adata.obs`, with values 'train' or 'val' indicating the assignment of each cell/sample to the training or validation set for that split.
3. **Gene Metadata** (`adata.var`): Optional gene information

> **Note:** If you use the provided example data, these columns are already included. If you use your own data, you must generate these columns yourself.

## Supported Methods

### MultiMIL Methods

- **multimil**: MultiMIL for sample-level prediction
- **multimil_reg**: MultiMIL for regression tasks

### Traditional Methods

- **pb_rf**: Random Forest on pseudobulk data
- **pb_nn**: Neural Network on pseudobulk data
- **gex_rf**: Random Forest on gene expression
- **gex_nn**: Neural Network on gene expression

### Frequency-based Methods

- **freq_rf**: Random Forest on frequency data
- **freq_nn**: Neural Network on frequency data
- **freq_mr**: Multiple Regression on frequency data

### Cell Type Methods

- **ct_pb_rf**: Cell type-aware Random Forest
- **ct_pb_nn**: Cell type-aware Neural Network
- **ct_pb_mr**: Cell type-aware Multiple Regression

## Parameter Documentation

Parameter files (TSV format) in `params/recommended_params/` specify the recommended hyperparameter search space for each method. Each file contains a header row and one or more rows of parameter values. The pipeline will run all parameter combinations listed in these files. You can edit these files to control which hyperparameters are tested for each method. For MultiMIL methods, you need to change `'your_sample_column_name', 'your_phenotype_column_name', 'your_celltype_column_name'` to the coresponding column names in your data. 

### MultiMIL Methods

#### multimil.tsv / multimil_reg.tsv
Parameters for MultiMIL classification and regression:

- **categorical_covariate_keys**: List of categorical covariates to include in the model (should include the prediction covariate)
- **class_loss_coef** (multimil.tsv) / **regression_loss_coef** (multimil_reg.tsv): Coefficient for loss function
- **train_max_epochs**: Maximum number of training epochs
- **train_save_checkpoint_every_n_epochs**: Frequency of checkpoint saving during training
- **lr**: Learning rate for model training
- **batch_size**: Batch size for training
- **seed**: Random seed for reproducibility
- **subset_umap**: Number of cells to subset for UMAP visualization
- **umap_colors**: List of variables to color UMAP plots by

We recommend to start the optimization with changing the learning rate as suggested in the `recommened_params`, and if the performance is unsatisfactory also optimize the loss coefficients.

### Random Forest and Multiple Regression Methods

#### pb_rf.tsv / pb_mr.tsv / gex_rf.tsv / gex_mr.tsv / freq_rf.tsv / freq_mr.tsv / ct_pb_rf.tsv / ct_pb_mr.tsv
Parameters for all Random Forest and Multiple Regression methods:

- **norm**: Whether to normalize the data using StandardScaler

### Neural Network Methods

#### pb_nn.tsv / gex_nn.tsv / freq_nn.tsv / ct_pb_nn.tsv
Parameters for all Neural Network methods:

- **norm**: Whether to normalize the data using StandardScaler
- **batch_size**: Batch size for training
- **lr**: Learning rate for model training
- **epochs**: Number of training epochs

## Running the Pipeline

The pipeline uses Snakemake for workflow management. If you're new to Snakemake, check out the [Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) for detailed guidance.

### Basic Usage

```bash
# Activate environment
mamba activate sample_prediction_pipeline

# Run the pipeline
snakemake --cores 4

# Or run specific targets
snakemake data/reports/your_task_accuracy.png --cores 4
```

## Output Interpretation

### Results Files

The pipeline generates several output files:

1. **Individual Results** (`data/reports/{task}/{method}/{hash}/accuracy.tsv`):
   - Performance metrics for each parameter combination
   - Cross-validation results

2. **Merged Results** (`data/reports/methods.tsv`):
   - Combined results from all methods and tasks

3. **Visualization** (`data/reports/{task}_accuracy.png`):
   - Performance comparison plots

## Getting Help

If you encounter issues:

1. Check the [GitHub Issues](https://github.com/theislab/sample-prediction-pipeline/issues) page
2. Create a new issue with detailed error messages