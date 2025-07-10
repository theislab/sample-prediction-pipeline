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
- Scanpy for single-cell analysis
- Pandas for data manipulation
- Matplotlib for plotting
- Other dependencies

### 3. Activate the Environment

```bash
mamba activate sample_prediction_pipeline
```

### 4. Verify Installation

```bash
# Test imports
python -c "import scanpy as sc; import pandas as pd; import matplotlib.pyplot as plt; print('All packages imported successfully!')"

# Test Snakemake
snakemake --version
```

## Quick Start

The repository includes ready-to-use example data and configuration files.

### 1. Use the Provided Example Data

1. **Copy the example config:**
   ```bash
   cp config_example.yaml config.yaml
   ```
2. **Run the pipeline:**
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

For each method, create a TSV file with parameter combinations (see `params/example_task/` for examples):

```tsv
learning_rate	epochs
0.01	10
0.001	20
```

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

## Running the Pipeline

### Basic Usage

```bash
# Activate environment
mamba activate sample_prediction_pipeline

# Run the pipeline
snakemake --cores 4

# Or run specific targets
snakemake data/reports/your_task_accuracy.png --cores 4
```

### Example Analyses

#### Classification Task

For cell type classification:

```yaml
TASKS:
  cell_type_classification:
    input: data/cell_types.h5ad
    label_key: cell_type
    condition_key: condition
    sample_key: sample_id
    n_splits: 5
    methods:
      multimil: params/cell_types/multimil.tsv
      pb_rf: params/cell_types/pb_rf.tsv
```

#### Regression Task

For continuous outcome prediction:

```yaml
TASKS:
  regression_analysis:
    input: data/regression_data.h5ad
    label_key: continuous_outcome
    condition_key: condition
    sample_key: sample_id
    n_splits: 5
    methods:
      multimil_reg: params/regression/multimil_reg.tsv
      pb_rf: params/regression/pb_rf.tsv
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

## Advanced Usage

### Regenerating Example Data

If you want to regenerate the example data (e.g., with different parameters), run:
```bash
python scripts/create_example_data.py
```
This will overwrite the files in `data/example/` and `params/example_task/`.

### Using Your Own Data

If you want to use your own data, ensure it is in H5AD format (AnnData object) with the following structure:

```python
import scanpy as sc
import pandas as pd

# Your data should have:
# - X: expression matrix (cells x genes)
# - obs: cell metadata (sample_id, condition, cell_type, etc.)
# - var: gene metadata (optional)

adata = sc.read_h5ad("your_data.h5ad")
```

## Getting Help

If you encounter issues:

1. Check the [GitHub Issues](https://github.com/theislab/sample-prediction-pipeline/issues) page
2. Create a new issue with detailed error messages