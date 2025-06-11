# Changelog

All notable changes to MetaGrouper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-XX-XX

### Added

#### Phase 1: K-mer Profiling and Sample Similarity Analysis
- Complete k-mer profiling system with canonical k-mer representation
- Support for single and paired-end FASTQ files (including gzipped)
- Multiple distance metrics: Bray-Curtis, Jaccard, Cosine, Euclidean
- PCA and MDS dimensionality reduction for visualization
- Interactive distance matrix heatmaps
- Memory-efficient k-mer extraction for large datasets
- Configurable k-mer sizes (11-31, odd numbers only)
- Comprehensive error handling and logging

#### Phase 2: Metadata Variable Analysis
- PERMANOVA (Permutational Multivariate Analysis of Variance) implementation
- Statistical significance testing with permutation-based p-values
- Variable importance ranking by R-squared values
- Support for both categorical and numerical metadata variables
- Automatic data type detection and preprocessing
- K-means and hierarchical clustering analysis
- Silhouette score optimization for cluster number selection
- PCA visualization colored by metadata variables
- Comprehensive statistical reporting in Markdown format

#### Phase 3: Assembly Strategy Recommendations
- Intelligent assembly strategy engine with confidence scoring
- Three strategy types: Individual, Grouped, and Global assembly
- Hybrid recommendation system combining similarity and metadata analysis
- Multi-assembler support: MEGAHIT, SPAdes, Flye
- Customizable similarity thresholds and group size constraints
- Performance prediction with contamination risk assessment
- Automated assembly command generation with optimized parameters
- Interactive strategy visualization with 4-panel overview
- Executable shell scripts for each assembly tool

#### Command Line Interface
- Comprehensive CLI with organized argument groups
- Performance optimization presets (quick, default, high-quality)
- Enhanced error handling and user-friendly error messages
- Progress reporting and runtime estimation
- System resource optimization recommendations
- Batch processing capabilities
- Configuration file support (JSON/YAML)

#### Documentation and Testing
- Complete README with installation and usage instructions
- Comprehensive tutorial with real-world examples
- Detailed API documentation
- Unit test suite covering all major components
- Integration tests for end-to-end workflows
- Performance benchmarking tools
- Example datasets and use cases

#### Output Formats
- Human-readable Markdown reports
- Machine-readable JSON outputs
- High-quality publication-ready plots (PNG, SVG)
- CSV data tables for further analysis
- Executable shell scripts for assembly tools
- Compressed output options for large datasets

### Technical Features

#### Performance Optimizations
- Memory-efficient k-mer processing with streaming
- Configurable read limits for large dataset testing
- Parallel distance matrix computation
- Optimized permutation testing algorithms
- Intelligent caching of intermediate results
- Resource usage monitoring and optimization

#### Robustness
- Comprehensive input validation
- Graceful handling of missing or malformed data
- Automatic detection of file formats and encodings
- Recovery from partial failures
- Detailed logging and debugging information
- Cross-platform compatibility (Windows, macOS, Linux)

#### Extensibility
- Modular architecture for easy extension
- Plugin system for custom distance metrics
- Configurable assembly tool parameters
- Custom visualization themes
- API for programmatic usage
- Integration with existing bioinformatics pipelines

### Dependencies
- Python 3.8+ support
- NumPy >= 1.21.0 for numerical computations
- Pandas >= 1.3.0 for data manipulation
- Scikit-learn >= 1.0.0 for machine learning algorithms
- SciPy >= 1.7.0 for statistical functions
- Matplotlib >= 3.5.0 for plotting
- Seaborn >= 0.11.0 for statistical visualization

### Known Issues
- Large k-mer sizes (>25) may require significant memory
- PERMANOVA with >10,000 permutations can be time-intensive
- Some assembly tools may require specific input formats

### Breaking Changes
- None (initial release)

## [Unreleased]

### Planned Features
- Integration with sourmash for scalable k-mer analysis
- Support for long-read sequencing data (PacBio, Oxford Nanopore)
- Machine learning-based assembly quality prediction
- Interactive web dashboard for result exploration
- Cloud computing integration (AWS, Google Cloud)
- Support for additional assembly tools (Unicycler, OPERA-MS)
- Taxonomic profiling integration
- Contamination detection algorithms
- Benchmarking against known reference datasets

### Under Consideration
- GPU acceleration for k-mer computation
- Distributed computing support for very large datasets
- Real-time analysis streaming
- Integration with sequence databases (NCBI, ENA)
- Phylogenetic analysis integration
- Multi-omics data integration capabilities