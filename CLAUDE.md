# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MetaGrouper is a bioinformatics tool for analyzing metagenomic samples through k-mer composition analysis, metadata statistical testing, and assembly strategy recommendations. The project is organized into 4 phases:

- **Phase 1**: Core k-mer profiling and similarity analysis
- **Phase 2**: Metadata association testing using PERMANOVA
- **Phase 3**: Assembly strategy recommendations for multiple assemblers
- **Phase 4**: Interactive HTML visualizations

## Key Commands

### Environment Setup
```bash
# Conda (recommended)
conda env create -f env.yaml
conda activate metagrouper

# Development installation
pip install -e ".[dev]"
pre-commit install
```

### Running Tests
```bash
# Run all tests
python -m pytest

# Run with coverage
python -m pytest --cov=metagrouper

# Run specific test
python -m pytest tests/test_kmer_profiling.py

# Run integration tests
python test_metagrouper.py
```

### Code Quality
```bash
# Format code with Black
black .

# Lint with Flake8
flake8 .

# Type checking
mypy metagrouper
```

### Building and Running
```bash
# Basic analysis
python metagrouper.py /path/to/fastq/files/ -o results/

# Full analysis with metadata and assembly recommendations
python metagrouper.py /path/to/fastq/files/ \
  --metadata samples_metadata.csv \
  --output results/ \
  --assembly-tools megahit spades \
  --similarity-threshold 0.45 \
  --permutations 999

# Preprocessing raw data
python preprocess.py raw_data/ -o clean_data/
```

## Recent Improvements (2025)

### Threshold Updates
MetaGrouper has been updated with research-based similarity thresholds:
- **similarity_threshold_high**: 0.15 → 0.25 (more permissive)
- **similarity_threshold_medium**: 0.30 → 0.45 (research-validated)
- **default_similarity_threshold**: 0.30 → 0.45 (better biological grouping)
- **max_group_size**: 10 → 20 (supports larger co-assemblies)

### Statistical Validation
- Added statistical power validation for PERMANOVA analysis
- Warns when groups have <10 samples for reliable results
- Improves scientific rigor of metadata associations

### Tool Documentation
- Clarified SPAdes limitations vs MEGAHIT for co-assembly
- Updated all configuration files and presets
- Improved assembly strategy recommendations

### Benefits
- More biologically meaningful sample groupings
- Better detection of similar samples for co-assembly
- Reduced false negative groupings
- Improved statistical robustness

## Architecture

### Core Modules
- `metagrouper.py` - Main entry point and k-mer profiling orchestration
- `metadata_analyzer.py` - Statistical analysis of metadata associations
- `assembly_recommender.py` - Assembly strategy recommendations
- `cli.py` - Enhanced command-line interface with progress tracking
- `config.py` - Configuration management and validation
- `preprocess.py` - FASTQ preprocessing pipeline

### Package Structure
The modular package (`metagrouper_package/metagrouper/`) contains:
- `profiler.py` - K-mer profiling implementation
- `analyzer.py` - Similarity analysis algorithms
- `visualizer.py` - Visualization generation
- `interactive_visualizer.py` - Interactive HTML visualizations (Phase 4)
- `utils.py` - Utility functions
- `config.py` - Configuration classes

### Key Design Patterns
1. **Worker Pool Pattern**: Used for parallel processing of FASTQ files
2. **Builder Pattern**: Configuration objects for complex analysis parameters
3. **Strategy Pattern**: Different assembly recommendation strategies
4. **Modular Architecture**: Each phase can be run independently

### Performance Optimizations
- Memory-efficient k-mer counting using canonical k-mers
- Parallel processing with multiprocessing.Pool
- Sparse matrix operations for large datasets
- Chunked file reading for memory efficiency

## Development Guidelines

### Code Standards
- Python 3.8+ required
- Use type hints for function signatures
- Follow Black formatting (line length 88)
- Docstrings for all public functions/classes
- Keep functions under 50 lines when possible

### Testing Requirements
- Write tests for new functionality
- Maintain test coverage above 80%
- Use pytest fixtures for common test data
- Test both success and error cases

### Common Development Tasks
When modifying core algorithms:
1. Check `test_metagrouper.py` for existing test patterns
2. Run performance tests with `test_parallel_performance.py`
3. Validate results against demo scripts (`demo*.py`)

When adding new features:
1. Update relevant phase documentation
2. Add command-line options to `cli.py` if needed
3. Update configuration schema in `config.py`
4. Add example usage to appropriate demo script