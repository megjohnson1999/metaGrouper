#!/usr/bin/env python3
"""
Setup script for MetaGrouper
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# Read requirements
def parse_requirements(filename):
    """Parse requirements file, excluding comments and dev dependencies."""
    requirements = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and not line.startswith("pytest"):
                requirements.append(line)
    return requirements

# Core requirements
install_requires = parse_requirements("requirements.txt")

# Development requirements
dev_requires = [
    "pytest>=7.0.0",
    "pytest-cov>=4.0.0",
    "black>=22.0.0",
    "flake8>=5.0.0",
    "mypy>=0.991",
    "pre-commit>=2.20.0"
]

setup(
    name="metagrouper",
    version="1.0.0",
    description="K-mer-based analysis for optimal metagenomic assembly grouping",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Megan Johnson",
    author_email="meganjohnson1w@gmail.com",
    url="https://github.com/megjohnson1999/metaGrouper",
    project_urls={
        "Bug Reports": "https://github.com/megjohnson1999/metaGrouper/issues",
        "Source": "https://github.com/megjohnson1999/metaGrouper",
        "Documentation": "https://github.com/megjohnson1999/metaGrouper/blob/main/README.md",
    },
    packages=find_packages(),
    py_modules=[
        "metagrouper",
        "metadata_analyzer", 
        "assembly_recommender",
        "cli",
        "config"
    ],
    python_requires=">=3.8",
    install_requires=install_requires,
    extras_require={
        "dev": dev_requires,
        "config": ["PyYAML>=6.0,<7.0"],
        "system": ["psutil>=5.8.0,<6.0.0"],
        "all": dev_requires + ["PyYAML>=6.0,<7.0", "psutil>=5.8.0,<6.0.0"]
    },
    entry_points={
        "console_scripts": [
            "metagrouper=metagrouper:main",
            "metagrouper-cli=cli:enhanced_main",
            "metagrouper-config=config:main",
        ],
    },
    scripts=[
        "metagrouper.py",
        "example_usage.py",
        "example_phase2.py", 
        "example_phase3.py"
    ],
    include_package_data=True,
    package_data={
        "": ["*.md", "*.txt", "*.yml", "*.yaml"]
    },
    keywords=[
        "bioinformatics",
        "metagenomics", 
        "assembly",
        "k-mer",
        "microbiome",
        "genomics",
        "sequencing",
        "PERMANOVA",
        "clustering"
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Natural Language :: English",
    ],
    license="MIT",
    zip_safe=False,
)