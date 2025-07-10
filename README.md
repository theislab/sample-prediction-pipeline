# Sample Prediction Pipeline

A Snakemake-based pipeline for benchmarking sample-level prediction methods from single-cell data, with a focus on Multiple Instance Learning (MultiMIL) and other machine learning approaches.

## Overview

This pipeline provides a comprehensive framework for comparing different methods for sample-level prediction from single-cell data, including:

- MultiMIL
- Random Forest (RF)
- Neural Networks (NN)
- Multiclass Regression (MR)

The pipeline supports both classification and regression tasks on single-cell data, and compares the performance on:
- bulk
- cell-type bulk
- MIL pooled
- mean pooled
sample representations.

## Features

- **Modular Design**: Easy to add new methods and datasets
- **Reproducible**: Snakemake ensures reproducible workflows
- **Flexible**: Supports various data formats and analysis tasks
- **Comprehensive Evaluation**: Multiple metrics and visualization options
- **Conda Environment**: Isolated environment management

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
cp config_example.yaml config.yaml
snakemake --cores 1
```

This will run all methods on the example dataset and produce results in `data/reports/`.

## Documentation

For detailed installation instructions, configuration options, and advanced usage, see the **[Complete Guide](docs/usage.md)**.

## Supported Methods

### MultiMIL
- **multimil**: MultiMIL for sample-level classification prediction
- **multimil_reg**: MultiMIL for sample-level regression prediction

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

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{sample_prediction_pipeline,
  title={Sample Prediction Pipeline: A Snakemake-based framework for sample-level prediction benchmarking},
  author={Anastasia Litinetskaya},
  year={2025},
  url={https://github.com/theislab/sample-prediction-pipeline}
}
```

## License

This project is licensed under the BSD-3-Clause License - see the [LICENSE](LICENSE) file for details.

## Support

For questions and support, please open an issue on GitHub or contact the maintainers. 