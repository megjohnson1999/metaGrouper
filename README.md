# MetaGrouper

<p align="center">
  <img src="metagrouper_logo.png" alt="MetaGrouper Logo" width="300"/>
</p>

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**MetaGrouper** analyzes metagenomic samples and recommends optimal assembly strategies based on k-mer composition similarity and metadata analysis.

## Key Features

- 🧬 **K-mer profiling** with fast sourmash MinHash sketching
- 📊 **Smart metadata filtering** focuses on biologically relevant variables
- 🎯 **Assembly recommendations** with ready-to-run commands
- 📈 **Interactive HTML reports** with dynamic visualizations
- ⚡ **Memory-efficient** processing for large datasets
- 🔧 **Multi-assembler support** (MEGAHIT, SPAdes, Flye)

## Recent Improvements (2025)

### 🚀 **High Sensitivity Sourmash Analysis** (NEW!)
- **10x higher sensitivity**: scaled=100 vs 1000 (retains 10x more k-mers)
- **Multi-scale analysis**: k=21,31,51 captures different similarity patterns
- **Presence/absence mode**: Robust to PCR bias (no abundance tracking by default)
- **Resolves "97-99% dissimilar" issues** from previous analyses
- **Scientifically robust** for metagenomic data with technical artifacts

### ✨ **Sample Name Normalization**
- **Automatic suffix removal** (`_hr`, `_trimmed`, `_filtered`) for proper metadata matching
- **Robust sample ID alignment** between FASTQ files and metadata
- **Resolves common data integration issues**

### 🔍 **Smart Metadata Variable Filtering**
- **Auto-filter biological variables** while excluding technical noise
- **Preserves patient IDs** (PID, GEMM) while removing lab identifiers
- **Configurable filtering** with `--auto-filter-variables` and `--exclude-variables`
- **Detailed filtering reports** explain decisions

### ⚡ **Performance Optimizations**
- **10-100x faster** k-mer analysis with sourmash MinHash
- **Constant memory usage** for large datasets
- **Parallel processing** with auto-detected CPU cores

## Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/megjohnson1999/metaGrouper.git
cd metaGrouper

# Install with conda (recommended)
conda env create -f env.yaml
conda activate metagrouper

# Test installation
python metagrouper.py --help
```

### Basic Usage

```bash
# Simple analysis
python metagrouper.py /path/to/fastq/files -o results/

# High sensitivity analysis (recommended)
python metagrouper.py /path/to/fastq/files \
    --scaled 100 \
    --additional-k-sizes 31 51 \
    -o results/

# With metadata and auto-filtering
python metagrouper.py /path/to/fastq/files \
    --metadata samples_metadata.csv \
    --auto-filter-variables \
    --output results/

# Full analysis with interactive report
python metagrouper.py /path/to/fastq/files \
    --metadata samples_metadata.csv \
    --auto-filter-variables \
    --assembly-tools megahit spades \
    --comprehensive-report \
    --output results/
```

## Input Requirements

### FASTQ Files
- Single or paired-end reads
- Gzipped files supported (`.fastq.gz`, `.fq.gz`)
- Files can be in subdirectories

### Metadata File (Optional)
CSV file with sample information:

```csv
Sample_ID,patient_id,case_control,Sex,month,Delivery_Mode
NovaSeq_N983_I13380_39894_Sample_01,P001,Case,Female,12,Vaginal
NovaSeq_N983_I13381_39895_Sample_02,P001,Case,Female,18,Vaginal
NovaSeq_N983_I13382_39896_Sample_03,P002,Control,Male,12,C-section
```

**Key Points:**
- Sample IDs must match FASTQ filenames (processing suffixes like `_hr` are automatically handled)
- Use `--sample-id-column Sample_ID` to specify the correct column name
- MetaGrouper will automatically filter out technical variables (Plate, Well, etc.)

## Command Line Options

### Core Options
- `input_dir` - Directory containing FASTQ files (required)
- `-o, --output` - Output directory (default: `metagrouper_output`)
- `-m, --metadata` - Metadata CSV file
- `-v, --verbose` - Verbose logging

### Smart Metadata Filtering
- `--auto-filter-variables` - **Automatically focus on biological variables**
- `--variables` - Specify variables manually (e.g., `--variables case_control Sex age`)
- `--exclude-variables` - Exclude specific variables (e.g., `--exclude-variables Plate Well`)

### K-mer Analysis
- `-k, --kmer-size` - K-mer size (default: 21)
- `--scaled` - Sourmash scaled parameter (default: 100 for high sensitivity)
- `--track-abundance` - Track k-mer abundances (disabled by default, more robust to PCR bias)
- `--additional-k-sizes` - Additional k-mer sizes for multi-scale analysis (e.g., 31 51)
- `--save-signatures` - Save sourmash signatures

### Assembly Recommendations
- `--assembly-tools` - Tools: `megahit`, `spades`, `flye` (default: `megahit spades`)
- `--similarity-threshold` - Grouping threshold (default: 0.45)
- `--min-group-size` - Minimum samples per group (default: 2)
- `--max-group-size` - Maximum samples per group (default: 20)

### Reports
- `--comprehensive-report` - Generate interactive HTML report
- `--html-title` - Title for HTML report
- `--permutations` - PERMANOVA permutations (default: 999)

## Output Files

### Core Results
- `distance_matrix.csv` - Sample similarity matrix
- `pca_plot.png` - PCA visualization
- `distance_heatmap.png` - Similarity heatmap
- `kmer_profiles.pkl` - K-mer profiles

### Metadata Analysis (if provided)
- `permanova_results.csv` - Statistical test results
- `variable_filtering_report.md` - Filtering decisions explained
- `variable_importance.png` - PERMANOVA results
- `analysis_report.md` - Summary report

### Assembly Recommendations
- `assembly_recommendations/` - Directory with detailed recommendations
  - `assembly_strategy.md` - Strategy summary
  - `run_megahit_assemblies.sh` - MEGAHIT commands
  - `run_spades_assemblies.sh` - SPAdes commands
- `assembly_strategy_overview.png` - Visual summary

### Interactive Report
- `interactive_report.html` - **Comprehensive interactive analysis**
  - Dynamic visualizations with zoom/pan/hover
  - Assembly strategy explanations
  - Professional publication-ready layout

## Examples

### Example 1: Auto-Filter Biological Variables
```bash
# Focus on biological variables, exclude technical noise
python metagrouper.py samples/ \
    --metadata patient_data.csv \
    --auto-filter-variables \
    --output clean_analysis/
```

### Example 2: Manual Variable Selection
```bash
# Analyze specific variables only
python metagrouper.py samples/ \
    --metadata patient_data.csv \
    --variables case_control Sex age delivery_mode \
    --output focused_analysis/
```

### Example 3: Large Dataset Analysis
```bash
# Efficient analysis with interactive report
python metagrouper.py large_dataset/ \
    --metadata samples_metadata.csv \
    --auto-filter-variables \
    --assembly-tools megahit spades \
    --comprehensive-report \
    --processes 8 \
    --output large_analysis/
```

### Example 4: Exclude Unwanted Variables
```bash
# Auto-filter but exclude specific variables
python metagrouper.py samples/ \
    --metadata data.csv \
    --auto-filter-variables \
    --exclude-variables batch_id processing_date \
    --output filtered_analysis/
```

## Interpreting Results

### Variable Filtering Report
The `variable_filtering_report.md` explains which variables were included/excluded:

- **Included**: Biological variables (disease, demographics, genetics)
- **Excluded**: Technical variables (Plate, Well, Barcode) and low-quality data

### Assembly Strategies
1. **Individual Assembly** - Each sample assembled separately (diverse samples)
2. **Grouped Assembly** - Samples grouped by similarity/metadata (balanced approach)
3. **Global Assembly** - All samples together (very similar samples)

### Confidence Scores
- **>0.8**: High confidence - strong recommendation
- **0.6-0.8**: Medium confidence - reasonable approach
- **0.4-0.6**: Low confidence - consider alternatives
- **<0.4**: Very low confidence - manual review needed

### Statistical Results
- **p < 0.05**: Significant metadata association
- **R² > 0.20**: Strong explanatory power
- **Multiple variables**: Compare R² values to prioritize

## Troubleshooting

### Common Issues

**"0 unique values, 0 non-null values" in metadata:**
- **Fixed!** Sample name normalization now handles this automatically
- Processing suffixes (`_hr`, `_trimmed`) are automatically stripped

**Metadata column not found:**
- Use `--sample-id-column` to specify correct column (e.g., `--sample-id-column Sample_ID`)
- Common alternatives: `sample_id`, `Sample_ID`, `sample`, `accession`

**Too many/few variables analyzed:**
- Use `--auto-filter-variables` to focus on biological variables
- Use `--variables` to specify exactly what you want
- Use `--exclude-variables` to remove unwanted variables

**Assembly recommendations seem poor:**
- Adjust `--similarity-threshold` (try 0.35-0.55)
- Check if metadata variables are meaningful
- Review the filtering report for excluded variables

### Getting Help

- 📚 **Documentation**: Check `RECOMMENDED_VARIABLES.md` for variable selection guidance
- 🐛 **Issues**: [GitHub Issues](https://github.com/megjohnson1999/metaGrouper/issues)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/megjohnson1999/metaGrouper/discussions)

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

- Built with scikit-learn, pandas, matplotlib, and sourmash
- Thanks to the sourmash team for fast MinHash implementation
- Thanks to all contributors and beta testers