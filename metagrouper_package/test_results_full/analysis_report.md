# MetaGrouper Phase 2: Metadata Analysis Report

## Variable Analysis Results

### Top Variables by Effect Size (R²)

- collection_date: R² = 0.995, p = 0.280 (23 groups, 24 samples)
- age: R² = 0.939, p = 0.440 (18 groups, 24 samples)
- **patient_id**: R² = 0.882, p = 0.000 (8 groups, 24 samples)
- location: R² = 0.601, p = 0.520 (3 groups, 24 samples)
- treatment: R² = 0.601, p = 0.540 (3 groups, 24 samples)
- timepoint: R² = 0.579, p = 0.760 (3 groups, 24 samples)
- **gender**: R² = 0.486, p = 0.040 (2 groups, 24 samples)
- bmi: R² = nan, p = nan (24 groups, 24 samples)
- viral_load: R² = nan, p = nan (24 groups, 24 samples)

### Statistical Significance

**2 variables** show significant associations (p < 0.05):
- patient_id: R² = 0.882, p = 0.000
- gender: R² = 0.486, p = 0.040

## Clustering Analysis

### Kmeans Clustering
- **Optimal clusters:** 5
- **Silhouette score:** 0.287

### Hierarchical Clustering
- **Optimal clusters:** 8
- **Silhouette score:** 0.215


## Interpretation Guidelines

### Effect Size (R²)
- **R² > 0.20:** Strong association (>20% variance explained)
- **R² = 0.10-0.20:** Moderate association
- **R² < 0.10:** Weak association

### Statistical Significance
- **p < 0.05:** Statistically significant
- **p < 0.01:** Highly significant
- **p < 0.001:** Very highly significant

### Clustering Quality
- **Silhouette > 0.70:** Strong clustering
- **Silhouette = 0.50-0.70:** Reasonable clustering
- **Silhouette < 0.50:** Weak clustering
