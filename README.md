# MetaGrouper

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**MetaGrouper** is a comprehensive tool for analyzing metagenomic samples and recommending optimal assembly strategies based on k-mer composition similarity and metadata analysis.

## Overview

MetaGrouper helps researchers make data-driven decisions about how to group metagenomic samples for assembly by:

1. **Analyzing k-mer composition** to identify similar samples
2. **Testing metadata variables** to find biologically meaningful groupings
3. **Recommending assembly strategies** with confidence scores and ready-to-run commands
4. **Supporting multiple assemblers** (MEGAHIT, SPAdes, Flye)

## Key Features

- 🧬 **K-mer profiling** with canonical k-mer representation
- 📊 **Statistical analysis** using PERMANOVA for metadata variables  
- 🎯 **Assembly strategy recommendations** (individual, grouped, or global)
- 🛠️ **Multi-assembler support** with optimized parameters
- 📈 **Comprehensive visualizations** and reports
- ⚡ **Memory-efficient** processing for large datasets
- 🔧 **Configurable parameters** for different use cases

## Installation

### Requirements

- Python 3.8 or higher
- Required packages listed in `requirements.txt`

### Quick Install

```bash
# Clone the repository
git clone https://github.com/megjohnson1999/metaGrouper.git
cd metaGrouper

# Install dependencies
pip install -r requirements.txt

# Test installation
python metagrouper.py --help
```

### Development Install

```bash
# Install in development mode
pip install -e .
```

## Quick Start

### Basic Usage (Phase 1 only)

```bash
# Analyze k-mer similarities only
python metagrouper.py /path/to/fastq/files/ -o results/
```

### With Metadata Analysis (Phases 1+2)

```bash
# Include metadata for comprehensive analysis
python metagrouper.py /path/to/fastq/files/ \
  --metadata samples_metadata.csv \
  --output results/ \
  --sample-id-column sample_id
```

### Full Analysis with Assembly Recommendations (All Phases)

```bash
# Complete analysis with assembly strategy
python metagrouper.py /path/to/fastq/files/ \
  --metadata samples_metadata.csv \
  --output results/ \
  --assembly-tools megahit spades \
  --similarity-threshold 0.25 \
  --permutations 999
```

## Input Requirements

### FASTQ Files
- Single or paired-end reads
- Gzipped files supported (`.fastq.gz`, `.fq.gz`)
- Files can be in subdirectories

### Metadata File (Optional)
CSV or TSV file with sample information:

```csv
sample_id,patient_id,timepoint,treatment,location
sample_001,P001,baseline,control,site_A
sample_002,P001,week4,treatment,site_A
sample_003,P002,baseline,control,site_B
```

**Required column:** Sample identifier (default: `sample_id`)

## Output Files

### Phase 1: K-mer Analysis
- `distance_heatmap.png` - Sample similarity heatmap
- `pca_plot.png` - PCA of k-mer profiles
- `mds_plot.png` - MDS of sample distances
- `distance_matrix.csv` - Pairwise distance matrix
- `kmer_profiles.pkl` - Serialized k-mer profiles

### Phase 2: Metadata Analysis
- `analysis_report.md` - Comprehensive statistical report
- `variable_importance.png` - PERMANOVA results visualization
- `permanova_results.csv` - Statistical test results
- `pca_by_*.png` - PCA plots colored by metadata variables
- `clustering_*.png` - Clustering analysis results

### Phase 3: Assembly Recommendations
- `assembly_recommendations/` - Directory with detailed recommendations
  - `assembly_strategy.md` - Human-readable strategy summary
  - `assembly_recommendations.json` - Machine-readable results
  - `run_megahit_assemblies.sh` - MEGAHIT assembly commands
  - `run_spades_assemblies.sh` - SPAdes assembly commands
  - `run_flye_assemblies.sh` - Flye assembly commands (if requested)
- `assembly_strategy_overview.png` - Visual strategy summary

## Command Line Options

### Input/Output
- `input_dir` - Directory containing FASTQ files (required)
- `-o, --output` - Output directory (default: `metagrouper_output`)
- `-v, --verbose` - Enable verbose logging

### K-mer Analysis (Phase 1)
- `-k, --kmer-size` - K-mer size (default: 21)
- `--max-reads` - Maximum reads per sample (for testing)
- `--distance-metric` - Distance metric: `braycurtis`, `jaccard`, `cosine`, `euclidean` (default: `braycurtis`)

### Metadata Analysis (Phase 2)
- `-m, --metadata` - Metadata file (CSV/TSV)
- `--sample-id-column` - Sample ID column name (default: `sample_id`)
- `--variables` - Specific variables to analyze (default: all)
- `--permutations` - Number of permutations for PERMANOVA (default: 999)
- `--cluster-range` - Range for cluster numbers to test (default: 2 8)

### Assembly Recommendations (Phase 3)
- `--assembly-tools` - Tools to generate commands for: `megahit`, `spades`, `flye`, `all` (default: `megahit spades`)
- `--similarity-threshold` - Distance threshold for grouping (default: 0.30)
- `--min-group-size` - Minimum samples per group (default: 2)
- `--max-group-size` - Maximum samples per group (default: 10)

## Examples

### Example 1: Quick Analysis
```bash
# Fast analysis with reduced parameters
python metagrouper.py example_data/ \
  --kmer-size 15 \
  --max-reads 1000 \
  --output quick_test/
```

### Example 2: Comprehensive Study
```bash
# Full analysis for publication
python metagrouper.py fastq_files/ \
  --metadata patient_metadata.csv \
  --output comprehensive_analysis/ \
  --kmer-size 21 \
  --permutations 9999 \
  --assembly-tools all \
  --similarity-threshold 0.20
```

### Example 3: Focus on Specific Variables
```bash
# Analyze only treatment and timepoint effects
python metagrouper.py samples/ \
  --metadata clinical_data.csv \
  --variables treatment timepoint \
  --permutations 999 \
  --output treatment_analysis/
```

## Interpreting Results

### Assembly Strategy Types

1. **Individual Assembly** - Each sample assembled separately
   - Recommended when samples are very diverse
   - No clear grouping patterns found
   - Low risk of contamination

2. **Grouped Assembly** - Samples grouped for co-assembly
   - Based on k-mer similarity and/or metadata
   - Balances benefits and contamination risk
   - Most common recommendation

3. **Global Assembly** - All samples co-assembled together
   - When all samples are very similar
   - Maximizes coverage and contiguity
   - Higher contamination risk

### Confidence Scores
- **0.8-1.0:** High confidence - strong evidence for strategy
- **0.6-0.8:** Medium confidence - reasonable evidence
- **0.4-0.6:** Low confidence - weak evidence, consider alternatives
- **<0.4:** Very low confidence - manual review recommended

### Statistical Significance
- **p < 0.05:** Statistically significant association
- **R² > 0.20:** High explanatory power
- **Silhouette > 0.5:** Good clustering quality

## Performance Optimization

### For Large Datasets
```bash
# Reduce computational load
python metagrouper.py large_dataset/ \
  --kmer-size 17 \
  --max-reads 5000 \
  --permutations 499 \
  --similarity-threshold 0.25
```

### For Testing/Development
```bash
# Fast testing parameters
python metagrouper.py test_data/ \
  --kmer-size 13 \
  --max-reads 500 \
  --permutations 99
```

## Troubleshooting

### Common Issues

**Memory errors with large datasets:**
- Reduce `--kmer-size` (try 17 or 15)
- Use `--max-reads` to limit input
- Process samples in smaller batches

**Slow performance:**
- Reduce `--permutations` for testing
- Use smaller k-mer sizes
- Check available RAM and CPU cores

**No significant metadata associations:**
- Check metadata formatting
- Ensure sample IDs match between files
- Consider different variables or groupings
- Increase sample size if possible

**Assembly groups seem suboptimal:**
- Adjust `--similarity-threshold`
- Review metadata variables
- Check distance matrix for patterns
- Consider manual grouping based on biology

### Error Messages

**"No FASTQ files found"**
- Check file extensions (`.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`)
- Verify directory path
- Ensure files are not empty

**"Metadata column not found"**
- Check `--sample-id-column` parameter
- Verify CSV/TSV format
- Ensure no extra spaces in column names

**"Phase X not available"**
- Check Python dependencies
- Install missing packages: `pip install -r requirements.txt`

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup
```bash
git clone https://github.com/megjohnson1999/metaGrouper.git
cd metaGrouper
pip install -e ".[dev]"
python test_metagrouper.py
```

## Citation

If you use MetaGrouper in your research, please cite:

```
MetaGrouper: K-mer-based analysis for optimal metagenomic assembly grouping
[Authors, Year, Journal]
DOI: [DOI]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- 📚 **Documentation:** Full documentation at [GitHub Repository](https://github.com/megjohnson1999/metaGrouper)
- 🐛 **Bug Reports:** [GitHub Issues](https://github.com/megjohnson1999/metaGrouper/issues)
- 💬 **Discussions:** [GitHub Discussions](https://github.com/megjohnson1999/metaGrouper/discussions)
- 📧 **Contact:** meganjohnson1w@gmail.com

## Acknowledgments

- Inspired by community needs in metagenomic analysis
- Built with scikit-learn, pandas, and matplotlib
- Thanks to all contributors and beta testers