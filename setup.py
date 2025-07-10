#!/usr/bin/env python3
"""
Setup script for Sample Prediction Pipeline.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# Read requirements from conda environment file
def get_requirements():
    """Extract requirements from conda environment file."""
    env_file = this_directory / "envs" / "sample_prediction_pipeline.yaml"
    requirements = []
    
    if env_file.exists():
        with open(env_file, 'r') as f:
            lines = f.readlines()
        
        in_dependencies = False
        for line in lines:
            line = line.strip()
            if line == 'dependencies:':
                in_dependencies = True
                continue
            elif line.startswith('-') and in_dependencies:
                if line.startswith('- pip:'):
                    break
                elif line.startswith('- python='):
                    continue
                else:
                    # Extract package name from conda format
                    package = line[2:]  # Remove '- '
                    if '=' in package:
                        package = package.split('=')[0]
                    requirements.append(package)
    
    return requirements

setup(
    name="sample-prediction-pipeline",
    version="0.1.0",
    author="Anastasia Litinetskaya",
    description="A Snakemake-based pipeline for benchmarking sample-level prediction methods from single-cell data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/theislab/sample-prediction-pipeline",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    python_requires=">=3.10",
    install_requires=get_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.800",
        ],
    },
    entry_points={
        "console_scripts": [
            "sample-prediction-pipeline=scripts.run_method:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.yaml", "*.yml", "*.tsv", "*.txt"],
    },
    keywords="bioinformatics, single-cell, machine-learning, snakemake, multiple-instance-learning",
    project_urls={
        "Bug Reports": "https://github.com/theislab/sample-prediction-pipeline/issues",
        "Source": "https://github.com/theislab/sample-prediction-pipeline",
        "Documentation": "https://github.com/theislab/sample-prediction-pipeline#readme",
    },
) 