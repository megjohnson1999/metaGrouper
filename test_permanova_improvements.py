#!/usr/bin/env python3
"""
Test script for the improved PERMANOVA analysis with all new features:
1. Categorical variable validation
2. Adaptive continuous variable binning
3. Adjusted R² calculation
4. Enhanced variable filtering
5. Statistical robustness checks
"""

import sys
import numpy as np
import pandas as pd
sys.path.append('phases')

from metadata_analyzer import MetadataAnalyzer, filter_metadata_variables

def create_test_data():
    """Create test data that demonstrates the improvements"""
    np.random.seed(42)
    
    # Create 200 samples
    n_samples = 200
    sample_names = [f"sample_{i:03d}" for i in range(n_samples)]
    
    # Create mock distance matrix
    distance_matrix = np.random.random((n_samples, n_samples))
    distance_matrix = (distance_matrix + distance_matrix.T) / 2  # Make symmetric
    np.fill_diagonal(distance_matrix, 0)  # Diagonal should be 0
    
    # Create test metadata with various variable types
    metadata = pd.DataFrame(index=sample_names)
    
    # 1. Good categorical variable (should pass validation)
    metadata['body_site'] = np.random.choice(['gut', 'oral', 'skin'], n_samples, p=[0.4, 0.3, 0.3])
    
    # 2. Individual-level identifier (should be flagged)
    metadata['patient_id'] = [f"patient_{i:03d}" for i in range(n_samples)]
    
    # 3. Too many small groups (should be flagged)
    metadata['rare_medication'] = np.random.choice([f"drug_{i}" for i in range(50)], n_samples)
    
    # 4. Continuous variable (should be binned adaptively)
    metadata['age'] = np.random.normal(45, 15, n_samples)
    
    # 5. Skewed continuous variable (like FCP)
    metadata['fcp_log'] = np.random.lognormal(4, 1, n_samples)
    
    # 6. Binary variable (should pass validation)
    metadata['disease_status'] = np.random.choice(['healthy', 'disease'], n_samples, p=[0.6, 0.4])
    
    # 7. Too few samples per group
    metadata['small_groups'] = np.random.choice(['A', 'B', 'C'], n_samples, p=[0.95, 0.03, 0.02])
    
    return distance_matrix, sample_names, metadata

def test_variable_filtering():
    """Test the improved variable filtering"""
    print("=== Testing Variable Filtering ===")
    
    _, _, metadata = create_test_data()
    
    # Test auto-filtering
    filtered_vars, exclusion_reasons = filter_metadata_variables(
        metadata, 
        auto_filter=True,
        max_unique_ratio=0.2
    )
    
    print(f"Original variables: {len(metadata.columns)}")
    print(f"Filtered variables: {len(filtered_vars)}")
    print(f"Excluded variables: {len(exclusion_reasons)}")
    
    print("\nIncluded variables:")
    for var in filtered_vars:
        print(f"  - {var}")
    
    print("\nExcluded variables:")
    for var, reason in exclusion_reasons.items():
        print(f"  - {var}: {reason}")
    
    print()

def test_categorical_validation():
    """Test categorical variable validation"""
    print("=== Testing Categorical Variable Validation ===")
    
    distance_matrix, sample_names, metadata = create_test_data()
    analyzer = MetadataAnalyzer(distance_matrix, sample_names)
    
    # Test different categorical variables
    test_vars = ['body_site', 'patient_id', 'rare_medication', 'disease_status', 'small_groups']
    
    for var in test_vars:
        if var in metadata.columns:
            validation = analyzer.validate_categorical_variable(metadata[var], var)
            print(f"{var}: {'VALID' if validation['valid'] else 'INVALID'}")
            print(f"  Reason: {validation['reason']}")
            print(f"  Stats: {validation['stats']}")
            print()

def test_adaptive_binning():
    """Test adaptive binning for continuous variables"""
    print("=== Testing Adaptive Binning ===")
    
    distance_matrix, sample_names, metadata = create_test_data()
    analyzer = MetadataAnalyzer(distance_matrix, sample_names)
    
    # Test binning on different distributions
    test_vars = ['age', 'fcp_log']
    
    for var in test_vars:
        if var in metadata.columns:
            values = metadata[var].values
            bins = analyzer._adaptive_binning(values, min_bin_size=30)
            unique_bins = np.unique(bins)
            
            print(f"{var}:")
            print(f"  - Distribution: mean={np.mean(values):.2f}, std={np.std(values):.2f}")
            print(f"  - Number of bins: {len(unique_bins)}")
            
            for i, bin_val in enumerate(unique_bins):
                count = np.sum(bins == bin_val)
                bin_mean = np.mean(values[bins == bin_val])
                bin_range = f"{np.min(values[bins == bin_val]):.1f}-{np.max(values[bins == bin_val]):.1f}"
                print(f"  - Bin {i}: {count} samples, mean={bin_mean:.2f}, range={bin_range}")
            print()

def test_full_analysis():
    """Test the full analysis pipeline with all improvements"""
    print("=== Testing Full Analysis Pipeline ===")
    
    distance_matrix, sample_names, metadata = create_test_data()
    analyzer = MetadataAnalyzer(distance_matrix, sample_names)
    # Set metadata directly instead of loading from file
    analyzer.metadata = metadata
    
    # Run analysis with auto-filtering
    results = analyzer.analyze_variables(
        auto_filter=True,
        n_permutations=99  # Reduced for testing
    )
    
    print(f"Analysis completed for {len(results)} variables")
    print("\nResults summary:")
    
    for _, row in results.head().iterrows():
        print(f"- {row['variable']}: R²={row['r_squared']:.3f}, ", end="")
        if 'r_squared_adj' in row:
            print(f"R²_adj={row['r_squared_adj']:.3f}, ", end="")
        if 'effect_size' in row:
            print(f"effect={row['effect_size']}, ", end="")
        print(f"p={row['p_value']:.3f}")
        
        if 'power_warning' in row and row['power_warning'] is not None:
            print(f"  ⚠️  {row['power_warning']}")
        if 'validation_warning' in row and pd.notna(row['validation_warning']):
            print(f"  ⚠️  {row['validation_warning']}")
    
    print("\n=== All tests completed successfully! ===")

if __name__ == "__main__":
    test_variable_filtering()
    test_categorical_validation()
    test_adaptive_binning()
    test_full_analysis()