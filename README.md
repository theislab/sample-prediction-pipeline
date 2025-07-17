# Sample Prediction Pipeline

A Snakemake-based pipeline for benchmarking sample-level prediction methods from single-cell data, with a focus on Multiple Instance Learning (MultiMIL) and other machine learning approaches.

## Overview

This pipeline provides a comprehensive framework for comparing different methods for sample-level prediction from single-cell data, including:

- MultiMIL
- Random Forest (RF)
- Multiclass Regression (MR)
- Feed-forward Neural Network (NN)

The pipeline supports both classification and regression tasks on single-cell data, and compares the performance on:
- bulk
- cell-type bulk
- frequency vector
- cell embedding
representations.

## Features

- Modular Design: Easy to add new methods and datasets
- Reproducible: Snakemake ensures reproducible workflows
- Comprehensive Evaluation: Multiple metrics and visualization options
- Conda Environment: Isolated environment management

## Quick Start

The repository includes ready-to-use example data and configuration files.

1. **Clone and setup:**
```bash
git clone https://github.com/theislab/sample-prediction-pipeline.git
cd sample-prediction-pipeline
mamba env create -f envs/sample_prediction_pipeline.yaml
mamba activate sample_prediction_pipeline
```

2. **Run the example:**
```bash
snakemake --cores 1
```

This will run all methods on the example dataset with minimal test parameters and produce results in `data/reports/`.

## Documentation

For detailed installation instructions, configuration options, and advanced usage, see the **[Complete Guide](docs/usage.md)**.

## Supported Methods

### MultiMIL
- **multimil**: MultiMIL for sample-level classification prediction
- **multimil_reg**: MultiMIL for sample-level regression prediction

### Pseudo-bulk Methods
- **pb_rf**: Random Forest on pseudobulk
- **pb_nn**: Neural Network on pseudobulk
- **pb_mr**: Multi-class Regression on pseudobulk
- 
### Cell Type  Methods
- **ct_pb_rf**: Random Forest on cell type-aware pseudobulk
- **ct_pb_nn**: Neural Network on cell type-aware pseudobulk 
- **ct_pb_mr**: Multi-class Regression on cell type-aware pseudobulk 

### Frequency-based Methods
- **freq_rf**: Random Forest on frequency data
- **freq_nn**: Neural Network on frequency data
- **freq_mr**: Multi-class Regression on frequency data

### Cell Embedding Methods
- **gex_rf**: Random Forest on cell embeddings
- **gex_nn**: Neural Network on cell embeddings
- **gex_mr**: Multi-class Regression on cell embeddings

## Output Structure

```
data/
├── reports/
│   ├── {task}/
│   │   ├── {method}/
│   │   │   └── {hash}/
│   │   │       └── accuracy.tsv
│   │   └── best_{method}.txt
│   ├── methods.tsv
│   ├── best.tsv
│   └── {task}_accuracy.png
└── tasks.tsv
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

<!-- ## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{sample_prediction_pipeline,
  title={Sample Prediction Pipeline: A Snakemake-based framework for sample-level prediction benchmarking},
  author={Anastasia Litinetskaya},
  year={2025},
  url={https://github.com/theislab/sample-prediction-pipeline}
}
``` -->

## License

This project is licensed under the BSD-3-Clause License - see the [LICENSE](LICENSE) file for details.

## Support

For questions and support, please open an issue on GitHub. 
