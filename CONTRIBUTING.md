# Contributing to Sample Prediction Pipeline

Thank you for your interest in contributing to the Sample Prediction Pipeline! This document provides guidelines for contributing to the project.

## How to Contribute

### Reporting Issues

Before creating a new issue, please:

1. Check if the issue has already been reported
2. Search the existing issues for similar problems
3. Provide detailed information about your environment and the problem

When reporting an issue, include:

- **Operating system** and version
- **Python version**
- **Conda/Mamba version**
- **Error messages** (full traceback)
- **Steps to reproduce** the problem
- **Expected vs actual behavior**

### Suggesting Enhancements

For feature requests:

1. Describe the feature you'd like to see
2. Explain why this feature would be useful
3. Provide examples of how it might work
4. Consider if you can implement it yourself

### Code Contributions

#### Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/theislab/sample-prediction-pipeline.git
   cd sample-prediction-pipeline
   ```
3. **Create a new branch** for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```
4. **Set up the development environment**:
   ```bash
   conda env create -f envs/sample_prediction_pipeline.yaml
   conda activate sample_prediction_pipeline
   ```

#### Development Guidelines

##### Code Style

- Follow **PEP 8** for Python code
- Use **type hints** where appropriate
- Write **docstrings** for all functions and classes
- Keep functions **small and focused**

##### Testing

- Write tests for new functionality
- Ensure existing tests pass
- Add integration tests for pipeline components

##### Documentation

- Update relevant documentation files
- Add docstrings to new functions
- Update README if needed
- Add examples for new features

#### Making Changes

1. **Make your changes** following the guidelines above
2. **Test your changes**:
   ```bash
   # Run the pipeline with example data
   python scripts/create_example_data.py
   snakemake --cores 2 --configfile config_example.yaml
   ```
3. **Commit your changes** with clear commit messages:
   ```bash
   git add .
   git commit -m "Add new feature: description of what was added"
   ```
4. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

#### Pull Request Process

1. **Create a pull request** from your fork to the main repository
2. **Fill out the PR template** with:
   - Description of changes
   - Related issues
   - Testing performed
   - Screenshots (if applicable)
3. **Wait for review** from maintainers
4. **Address feedback** and make requested changes
5. **Maintainers will merge** when ready

## Development Setup

### Local Development Environment

   ```bash
   # Clone the repository
   git clone https://github.com/theislab/sample-prediction-pipeline.git
   cd sample-prediction-pipeline

# Create development environment
conda env create -f envs/sample_prediction_pipeline.yaml
conda activate sample_prediction_pipeline

# Install development dependencies
pip install pytest black flake8 mypy

# Install in development mode
pip install -e .
```

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_utils.py

# Run with coverage
pytest --cov=scripts
```

### Code Quality Checks

```bash
# Format code
black scripts/

# Check code style
flake8 scripts/

# Type checking
mypy scripts/
```

## Project Structure

```
sample-prediction-pipeline/
├── config.yaml              # Main configuration
├── Snakefile                # Snakemake workflow
├── scripts/                 # Pipeline scripts
│   ├── utils.py            # Utility functions
│   ├── run_method.py       # Main execution script
│   └── ...                 # Method-specific scripts
├── envs/                   # Conda environments
├── params/                 # Parameter files
├── data/                   # Data directory
├── docs/                   # Documentation
└── tests/                  # Test files
```

## Adding New Methods

To add a new method to the pipeline:

1. **Create a new script** in `scripts/`:
   ```python
   # scripts/your_method.py
   def run_your_method(adata, **kwargs):
       """Run your method on the data."""
       # Your implementation here
       return results_df
   ```

2. **Add to METHOD_MAP** in `scripts/run_method.py`:
   ```python
   METHOD_MAP = dict(
       # ... existing methods ...
       your_method=dict(function=run_your_method, mode='rna'),
   )
   ```

3. **Create parameter files** in `params/`:
   ```tsv
   param1	param2	param3
   value1	value2	value3
   value4	value5	value6
   ```

4. **Update documentation** with method description

5. **Add tests** for the new method

## Code Review Guidelines

### For Contributors

- **Respond promptly** to review comments
- **Be open to feedback** and suggestions
- **Explain your reasoning** when defending design decisions
- **Keep commits focused** and atomic

### For Reviewers

- **Be constructive** and respectful
- **Focus on the code**, not the person
- **Provide specific feedback** with examples
- **Consider the bigger picture** and maintainability

## Release Process

### Versioning

We use [Semantic Versioning](https://semver.org/):

- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality (backward compatible)
- **PATCH**: Bug fixes (backward compatible)

### Release Checklist

Before releasing:

- [ ] All tests pass
- [ ] Documentation is updated
- [ ] CHANGELOG is updated
- [ ] Version number is updated
- [ ] Release notes are prepared

## Communication

### Getting Help

- **GitHub Issues**: For bugs and feature requests
- **GitHub Discussions**: For questions and general discussion
- **Email**: For private or sensitive matters

### Community Guidelines

- **Be respectful** and inclusive
- **Help others** when you can
- **Share knowledge** and experiences
- **Follow the code of conduct**

## Recognition

Contributors will be recognized in:

- **README.md** contributors section
- **Release notes** for significant contributions
- **GitHub contributors** page

## License

By contributing to this project, you agree that your contributions will be licensed under the same license as the project (BSD-3-Clause).

## Questions?

If you have questions about contributing, please:

1. Check the documentation
2. Search existing issues
3. Create a new issue with the "question" label
4. Start a discussion in GitHub Discussions

Thank you for contributing to the Sample Prediction Pipeline! 