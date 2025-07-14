# Changelog

All notable changes to the Sample Prediction Pipeline project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive parameter files in `params/recommended_params/` for all supported methods
- Parameter optimization framework with multiple combinations for each method
- Documentation for parameter files and hyperparameter tuning in usage guide
- Support for NN methods

### Changed
- Updated `norm` parameter from `0` to `False` in all parameter files for consistency
- Duplicated parameter combinations with `norm=True` to test both normalization settings
- Filtered MultiMIL parameter files to keep only `gated_attn` scoring method
- Updated MultiMIL parameter files with generic column names for easier customization
- Enhanced parameter documentation with usage instructions and optimization recommendations

### Fixed
- N/A

### Security
- N/A

---

## [0.1.0] - 2025-07-09

### Added
- Initial public release
- Snakemake-based pipeline for sample-level prediction benchmarking
- Support for MultiMIL, Random Forest, and Multiclass Regression
- Comprehensive documentation and examples
- Parameter optimization framework
- Cross-validation support
- Performance visualization tools

### Changed
- N/A

### Deprecated
- N/A

### Removed
- N/A

### Fixed
- N/A

### Security
- N/A

---

## Version History

- **0.1.0**: Initial development version with basic functionality
- **Unreleased**: First public release with comprehensive documentation and examples

## Contributing

To add entries to this changelog:

1. Add your changes under the appropriate section in `[Unreleased]`
2. Use the following categories:
   - **Added**: New features
   - **Changed**: Changes in existing functionality
   - **Deprecated**: Soon-to-be removed features
   - **Removed**: Removed features
   - **Fixed**: Bug fixes
   - **Security**: Security-related changes

3. When releasing, move `[Unreleased]` content to a new version section
4. Update the version number and date

## Format

Each version section should include:

- **Version number** and release date
- **Categories** of changes
- **Descriptive entries** for each change
- **Links** to relevant issues or pull requests (if applicable) 