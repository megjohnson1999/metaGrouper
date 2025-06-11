# Contributing to MetaGrouper

We welcome contributions to MetaGrouper! This document provides guidelines for contributing to the project.

## Code of Conduct

This project adheres to a code of conduct. By participating, you are expected to uphold this code. Please report unacceptable behavior to the project maintainers.

## How to Contribute

### Reporting Bugs

Before creating bug reports, please check the existing issues to avoid duplicates. When creating a bug report, include:

- **Clear title and description**
- **Steps to reproduce** the issue
- **Expected vs actual behavior**
- **Environment details** (OS, Python version, MetaGrouper version)
- **Sample data** or minimal reproducible example
- **Log files** if applicable

### Suggesting Features

Feature requests are welcome! Please include:

- **Clear description** of the feature
- **Use case** and motivation
- **Proposed implementation** if you have ideas
- **Alternatives considered**

### Pull Requests

1. **Fork the repository** and create your branch from `main`
2. **Install development dependencies**:
   ```bash
   pip install -e ".[dev]"
   ```
3. **Make your changes** following our coding standards
4. **Add tests** for new functionality
5. **Run the test suite**:
   ```bash
   python -m pytest tests/
   ```
6. **Update documentation** if needed
7. **Commit with clear messages**
8. **Submit a pull request**

## Development Setup

### Prerequisites

- Python 3.8 or higher
- Git
- Virtual environment (recommended)

### Installation

```bash
# Clone your fork
git clone https://github.com/yourusername/metaGrouper.git
cd metaGrouper

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### Running Tests

```bash
# Run all tests
python -m pytest

# Run with coverage
python -m pytest --cov=metagrouper

# Run specific test file
python -m pytest tests/test_kmer_profiling.py

# Run integration tests
python test_metagrouper.py
```

### Code Style

We use several tools to maintain code quality:

- **Black** for code formatting
- **Flake8** for linting
- **MyPy** for type checking
- **Pre-commit hooks** for automated checks

```bash
# Format code
black .

# Check linting
flake8 .

# Type checking
mypy metagrouper.py metadata_analyzer.py assembly_recommender.py
```

## Coding Standards

### Python Style

- Follow PEP 8 with Black formatting
- Use type hints where appropriate
- Write docstrings for all public functions and classes
- Keep functions focused and under 50 lines when possible
- Use descriptive variable names

### Documentation

- Update README.md for user-facing changes
- Update TUTORIAL.md for new features
- Include docstrings with examples:

```python
def example_function(param1: str, param2: int = 10) -> bool:
    """
    Brief description of what the function does.
    
    Args:
        param1: Description of parameter
        param2: Description with default value
        
    Returns:
        Description of return value
        
    Example:
        >>> result = example_function("test", 5)
        >>> print(result)
        True
    """
    pass
```

### Testing

- Write unit tests for new functions
- Include integration tests for new workflows
- Test edge cases and error conditions
- Use meaningful test names:

```python
def test_kmer_profiler_handles_empty_fastq_file():
    """Test that KmerProfiler gracefully handles empty FASTQ files."""
    pass
```

## Project Structure

```
metagrouper/
├── metagrouper.py          # Main module (Phase 1)
├── metadata_analyzer.py   # Phase 2 functionality
├── assembly_recommender.py # Phase 3 functionality
├── cli.py                 # Enhanced CLI
├── config.py              # Configuration management
├── test_metagrouper.py    # Test suite
├── demo.py                # Demonstration script
├── example_*.py           # Example scripts
├── requirements.txt       # Dependencies
├── setup.py              # Package setup
├── README.md             # Main documentation
├── TUTORIAL.md           # User tutorial
├── CONTRIBUTING.md       # This file
├── CHANGELOG.md          # Version history
└── LICENSE               # MIT license
```

## Release Process

1. Update version in `setup.py`
2. Update `CHANGELOG.md`
3. Run full test suite
4. Create release branch
5. Submit pull request
6. Tag release after merge
7. Update PyPI package

## Areas for Contribution

### High Priority

- **Performance optimization** for large datasets
- **Memory usage** improvements
- **Additional distance metrics**
- **More assembly tool integrations**
- **Better error messages**

### Medium Priority

- **Interactive visualizations** (Plotly/Bokeh)
- **Web interface** development
- **Cloud integration** (AWS, GCP)
- **Database backends** for large studies
- **Plugin system** for custom analyzers

### Documentation

- **Video tutorials**
- **Best practices guide**
- **Performance benchmarks**
- **Use case studies**
- **API documentation**

## Questions?

- **GitHub Issues**: For bugs and feature requests
- **GitHub Discussions**: For questions and general discussion
- **Email**: Contact maintainers directly

## Recognition

Contributors will be recognized in:
- `CHANGELOG.md` for their contributions
- GitHub contributors page
- Release notes for significant contributions

Thank you for contributing to MetaGrouper!