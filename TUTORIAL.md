# MetaGrouper Tutorial

This tutorial provides step-by-step guidance for using MetaGrouper to analyze metagenomic samples and determine optimal assembly strategies.

## Table of Contents

1. [Getting Started](#getting-started)
2. [Understanding Your Data](#understanding-your-data)
3. [Phase 1: K-mer Analysis](#phase-1-k-mer-analysis)
4. [Phase 2: Metadata Analysis](#phase-2-metadata-analysis)
5. [Phase 3: Assembly Recommendations](#phase-3-assembly-recommendations)
6. [Real-World Examples](#real-world-examples)
7. [Best Practices](#best-practices)

## Getting Started

### Prerequisites

Ensure you have:
- Python 3.8+ installed
- MetaGrouper dependencies: `pip install -r requirements.txt`
- FASTQ files from your metagenomic samples
- (Optional) Metadata file with sample information

### Test Installation

```bash
# Run the example to verify installation
python example_usage.py
```

This should generate example data and produce analysis results in `example_output/`.

## Understanding Your Data

### FASTQ File Organization

MetaGrouper expects FASTQ files in a directory structure like:

```
fastq_files/
├── sample_001.fastq.gz
├── sample_002.fastq.gz
├── sample_003.fastq
└── subdirectory/
    ├── sample_004.fq.gz
    └── sample_005.fq
```

**Supported formats:**
- `.fastq`, `.fq` (uncompressed)
- `.fastq.gz`, `.fq.gz` (gzipped)
- Files can be in subdirectories

### Metadata File Format

Create a CSV or TSV file with sample information:

```csv
sample_id,patient_id,timepoint,treatment_group,sampling_site,age_group
sample_001,P001,baseline,control,gut,adult
sample_002,P001,week4,control,gut,adult
sample_003,P002,baseline,treatment,gut,adult
sample_004,P002,week4,treatment,gut,adult
sample_005,P003,baseline,control,skin,elderly
```

**Key requirements:**
- Must include a sample ID column matching FASTQ filenames
- Can include categorical variables (treatment, site, etc.)
- Can include numerical variables (age, BMI, etc.)
- Missing values are handled automatically

## Phase 1: K-mer Analysis

### Basic K-mer Analysis

Start with a simple analysis to understand sample similarities:

```bash
python metagrouper.py fastq_files/ --output phase1_results/
```

### Interpreting Phase 1 Results

#### Distance Heatmap (`distance_heatmap.png`)
- **Dark colors:** Similar samples (low distance)
- **Light colors:** Dissimilar samples (high distance)
- **Clusters:** Groups of similar samples appear as dark blocks

#### PCA Plot (`pca_plot.png`)
- **Close points:** Similar k-mer composition
- **Distant points:** Different k-mer composition
- **Variance explained:** Shown in axis labels (higher is better)

#### Key Files
- `distance_matrix.csv`: Quantitative similarity data
- `kmer_profiles.pkl`: Raw k-mer data for further analysis

### Optimizing K-mer Parameters

#### Choosing K-mer Size

```bash
# Smaller k-mers (faster, less specific)
python metagrouper.py fastq_files/ --kmer-size 15 --output k15_results/

# Default k-mers (balanced)
python metagrouper.py fastq_files/ --kmer-size 21 --output k21_results/

# Larger k-mers (slower, more specific)
python metagrouper.py fastq_files/ --kmer-size 27 --output k27_results/
```

**Guidelines:**
- **k=15-17:** Fast analysis, good for initial exploration
- **k=21:** Default, works well for most datasets
- **k=25-31:** More specific, better for closely related samples

#### Distance Metrics

```bash
# Bray-Curtis (default, good for abundance data)
python metagrouper.py fastq_files/ --distance-metric braycurtis

# Jaccard (presence/absence, good for diverse samples)
python metagrouper.py fastq_files/ --distance-metric jaccard

# Cosine (angle-based, good for normalized data)
python metagrouper.py fastq_files/ --distance-metric cosine
```

## Phase 2: Metadata Analysis

### Adding Metadata

```bash
python metagrouper.py fastq_files/ \
  --metadata sample_metadata.csv \
  --sample-id-column sample_id \
  --output phase2_results/
```

### Understanding PERMANOVA Results

#### Variable Importance Plot (`variable_importance.png`)
- **Red bars:** Significant variables (p < 0.05)
- **Orange bars:** Marginally significant (p < 0.1)
- **Gray bars:** Non-significant variables
- **Bar length:** Proportion of variation explained (R²)

#### Statistical Interpretation

From `permanova_results.csv`:

```csv
variable,r_squared,p_value,f_statistic,n_samples,n_groups
treatment_group,0.45,0.001,8.2,20,2
sampling_site,0.32,0.003,6.1,20,3
timepoint,0.18,0.065,2.8,20,4
patient_id,0.15,0.124,1.9,20,8
```

**Interpretation:**
- `treatment_group`: Explains 45% of variation, highly significant
- `sampling_site`: Explains 32% of variation, significant
- `timepoint`: Explains 18% of variation, marginally significant
- `patient_id`: Explains 15% of variation, not significant

### PCA by Metadata

The `pca_by_*.png` plots show how samples cluster by each variable:

- **Clear separation:** Variable strongly affects sample composition
- **Overlapping clusters:** Variable has weak effect
- **Gradual transitions:** Continuous variables with gradual effects

### Advanced Metadata Analysis

#### Focusing on Specific Variables

```bash
python metagrouper.py fastq_files/ \
  --metadata sample_metadata.csv \
  --variables treatment_group sampling_site \
  --output focused_analysis/
```

#### Increasing Statistical Power

```bash
python metagrouper.py fastq_files/ \
  --metadata sample_metadata.csv \
  --permutations 9999 \
  --output high_power_analysis/
```

## Phase 3: Assembly Recommendations

### Getting Assembly Recommendations

```bash
python metagrouper.py fastq_files/ \
  --metadata sample_metadata.csv \
  --assembly-tools megahit spades \
  --similarity-threshold 0.25 \
  --output assembly_analysis/
```

### Understanding Assembly Strategies

#### Individual Assembly
```
Strategy: Individual
Confidence: 0.8
Rationale: Samples are too diverse for effective co-assembly
```

**When recommended:**
- High inter-sample distances
- No significant metadata associations
- Very diverse sample collection

**Commands generated:**
```bash
megahit -r sample_001.fastq -o sample_001_assembly --min-contig-len 500
megahit -r sample_002.fastq -o sample_002_assembly --min-contig-len 500
```

#### Grouped Assembly
```
Strategy: Grouped
Confidence: 0.7
Primary Criterion: treatment_group
Groups: 2 co-assembly groups covering 18/20 samples
```

**When recommended:**
- Clear metadata-based groupings
- Moderate inter-sample distances
- Balance between benefits and contamination risk

**Example groups:**
- Group 1: Control samples (n=10)
- Group 2: Treatment samples (n=8)
- Individual: Outlier samples (n=2)

#### Global Assembly
```
Strategy: Global
Confidence: 0.9
Rationale: All samples show strong similarity
```

**When recommended:**
- Very low inter-sample distances
- All samples from similar conditions
- Maximizing coverage and contiguity is priority

### Customizing Assembly Parameters

#### Similarity Threshold

```bash
# More stringent grouping (fewer, more similar groups)
python metagrouper.py fastq_files/ \
  --similarity-threshold 0.15 \
  --output stringent_groups/

# More permissive grouping (more, less similar groups)
python metagrouper.py fastq_files/ \
  --similarity-threshold 0.40 \
  --output permissive_groups/
```

#### Group Size Constraints

```bash
python metagrouper.py fastq_files/ \
  --min-group-size 3 \
  --max-group-size 8 \
  --output constrained_groups/
```

#### Assembly Tool Selection

```bash
# MEGAHIT only (fast, good for most cases)
python metagrouper.py fastq_files/ --assembly-tools megahit

# SPAdes only (slower, higher quality)
python metagrouper.py fastq_files/ --assembly-tools spades

# All tools (comprehensive comparison)
python metagrouper.py fastq_files/ --assembly-tools all
```

## Real-World Examples

### Example 1: Human Gut Microbiome Study

**Scenario:** 50 samples from patients before/after antibiotic treatment

```bash
python metagrouper.py gut_samples/ \
  --metadata patient_metadata.csv \
  --variables patient_id timepoint antibiotic_type \
  --kmer-size 21 \
  --assembly-tools megahit spades \
  --similarity-threshold 0.25 \
  --output gut_microbiome_analysis/
```

**Expected results:**
- Strong effect of `timepoint` (pre/post treatment)
- Moderate effect of `antibiotic_type`
- Weak effect of `patient_id` (individual variation)
- Recommendation: Group by timepoint and antibiotic type

### Example 2: Environmental Metagenomes

**Scenario:** Water samples from different lakes and seasons

```bash
python metagrouper.py water_samples/ \
  --metadata environmental_data.csv \
  --variables lake_name season temperature pH \
  --kmer-size 19 \
  --assembly-tools megahit \
  --similarity-threshold 0.30 \
  --output environmental_analysis/
```

**Expected results:**
- Strong effect of `lake_name` (geographic separation)
- Moderate effect of `season` (temporal variation)
- Weak effects of chemical parameters
- Recommendation: Group by lake, possibly subdivide by season

### Example 3: Agricultural Soil Study

**Scenario:** Soil samples from different crop treatments and fields

```bash
python metagrouper.py soil_samples/ \
  --metadata field_experiment.csv \
  --variables crop_type fertilizer_treatment field_location \
  --kmer-size 21 \
  --assembly-tools spades \
  --permutations 9999 \
  --output agriculture_analysis/
```

## Best Practices

### Data Preparation

1. **Quality Control:** Pre-filter low-quality reads before analysis
2. **Consistent Naming:** Use consistent sample naming between FASTQ files and metadata
3. **Metadata Completeness:** Include all relevant experimental variables
4. **Sample Size:** Aim for at least 6-10 samples per group for statistical power

### Parameter Selection

1. **Start Simple:** Begin with default parameters for initial exploration
2. **Biological Relevance:** Choose similarity thresholds based on your system
3. **Computational Resources:** Balance accuracy with available time/memory
4. **Validation:** Use multiple k-mer sizes and distance metrics for robust results

### Interpreting Results

1. **Statistical Significance:** Don't ignore p-values, but consider biological relevance
2. **Effect Sizes:** R² values indicate practical importance
3. **Confidence Scores:** Lower confidence suggests need for manual review
4. **Biological Knowledge:** Combine statistical results with domain expertise

### Assembly Strategy Selection

1. **Conservative Approach:** When in doubt, prefer individual assembly
2. **Contamination Risk:** Consider the cost of cross-sample contamination
3. **Computational Constraints:** Balance quality with available resources
4. **Downstream Analysis:** Consider how assembly strategy affects your research goals

### Troubleshooting

#### Low Confidence Scores
- Check sample quality and diversity
- Review metadata variables for relevance
- Consider increasing sample size
- Adjust similarity thresholds

#### No Clear Groupings
- Examine distance matrix for patterns
- Try different k-mer sizes
- Consider different distance metrics
- Review experimental design

#### Unexpected Groupings
- Validate with biological knowledge
- Check for batch effects in metadata
- Examine individual sample quality
- Consider technical confounding factors

## Advanced Usage

### Batch Processing

For large studies, consider processing in batches:

```bash
# Process by experimental batch
for batch in batch1 batch2 batch3; do
  python metagrouper.py ${batch}_samples/ \
    --metadata ${batch}_metadata.csv \
    --output ${batch}_results/
done
```

### Integration with Assembly Pipelines

Use MetaGrouper outputs in automated pipelines:

```bash
# Run MetaGrouper
python metagrouper.py samples/ --output analysis/

# Execute recommended assemblies
bash analysis/assembly_recommendations/run_megahit_assemblies.sh

# Continue with downstream analysis
```

### Custom Analysis

For specialized needs, use the Python modules directly:

```python
from metagrouper import KmerProfiler, SimilarityAnalyzer
from metadata_analyzer import MetadataAnalyzer
from assembly_recommender import AssemblyRecommender

# Custom analysis workflow
profiler = KmerProfiler(k=21)
# ... custom analysis code
```

This tutorial should help you get started with MetaGrouper and make the most of its capabilities for your metagenomic research!