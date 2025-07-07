# MetaGrouper Phase 2: Metadata Analysis Report

## Variable Importance (PERMANOVA Results)

Variables ranked by proportion of variation explained (R-squared):

- **collection_date**: R² = 0.958, p = 0.204
  - Type: categorical
  - Valid samples: 24
  - Groups: 23

- **patient_id**: R² = 0.343, p = 0.001**
  - Type: categorical
  - Valid samples: 24
  - Groups: 8

- **viral_load**: R² = 0.088, p = 0.309
  - Type: numerical
  - Valid samples: 24
  - Groups: 3

- **treatment**: R² = 0.085, p = 0.726
  - Type: categorical
  - Valid samples: 24
  - Groups: 3

- **location**: R² = 0.085, p = 0.702
  - Type: categorical
  - Valid samples: 24
  - Groups: 3

- **bmi**: R² = 0.084, p = 0.750
  - Type: numerical
  - Valid samples: 24
  - Groups: 3

- **age**: R² = 0.084, p = 0.785
  - Type: numerical
  - Valid samples: 24
  - Groups: 3

- **timepoint**: R² = 0.082, p = 0.910
  - Type: categorical
  - Valid samples: 24
  - Groups: 3

- **gender**: R² = 0.050, p = 0.032*
  - Type: categorical
  - Valid samples: 24
  - Groups: 2

## Clustering Analysis

### Kmeans Clustering

- **Optimal clusters**: 2
- **Silhouette score**: 0.062

Silhouette scores for different k values:
- k=2: 0.062
- k=3: 0.045
- k=4: 0.041
- k=5: 0.042
- k=6: 0.031
- k=7: 0.024
- k=8: 0.031

### Hierarchical Clustering

- **Optimal clusters**: 2
- **Silhouette score**: 0.065

Silhouette scores for different k values:
- k=2: 0.065
- k=3: 0.020
- k=4: 0.004
- k=5: -0.020
- k=6: -0.014
- k=7: -0.020
- k=8: -0.023

## Recommendations

### Primary Grouping Variable

**collection_date** explains the most variation (95.8%) in sample similarities.

This association is not statistically significant (p ≥ 0.05).

**Recommendation**: Consider individual sample assembly or global co-assembly.

## Interpretation Guide

- **R-squared**: Proportion of variation in sample composition explained by the variable
- **p-value**: Statistical significance (< 0.05 is typically significant)
- **Silhouette score**: Quality of clustering (higher is better, > 0.5 is good)
- Significance codes: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1
