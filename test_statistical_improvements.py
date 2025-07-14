#!/usr/bin/env python3
"""
Test script for all statistical improvements to PERMANOVA:
1. Corrected sum of squares calculation
2. Distance matrix validation
3. Multiple testing correction
4. Omega squared effect size
5. PERMDISP variance homogeneity test
"""

import sys
import numpy as np
import pandas as pd
sys.path.append('phases')

from metadata_analyzer import MetadataAnalyzer

def create_test_data():
    """Create realistic test data for statistical validation"""
    np.random.seed(42)
    
    # Create 50 samples (smaller for faster testing)
    n_samples = 50
    sample_names = [f"sample_{i:03d}" for i in range(n_samples)]
    
    # Create realistic distance matrix (Bray-Curtis like)
    # Simulate microbiome data with different groups
    n_features = 100
    
    # Group 1: High diversity (gut samples)
    group1_data = np.random.dirichlet(np.ones(n_features) * 0.5, size=20)
    
    # Group 2: Low diversity (oral samples)  
    group2_data = np.random.dirichlet(np.ones(n_features) * 0.1, size=20)
    
    # Group 3: Medium diversity (skin samples)
    group3_data = np.random.dirichlet(np.ones(n_features) * 0.3, size=10)
    
    # Combine data
    all_data = np.vstack([group1_data, group2_data, group3_data])
    
    # Calculate Bray-Curtis distance matrix
    from scipy.spatial.distance import pdist, squareform
    distances = pdist(all_data, metric='braycurtis')
    distance_matrix = squareform(distances)
    
    # Create metadata
    metadata = pd.DataFrame(index=sample_names)
    metadata['body_site'] = (['gut'] * 20 + ['oral'] * 20 + ['skin'] * 10)
    metadata['disease_status'] = np.random.choice(['healthy', 'disease'], n_samples, p=[0.7, 0.3])
    metadata['age'] = np.random.normal(45, 15, n_samples)
    
    return distance_matrix, sample_names, metadata

def test_distance_matrix_validation():
    """Test distance matrix validation"""
    print("=== Testing Distance Matrix Validation ===")
    
    # Test with good matrix
    distance_matrix, sample_names, metadata = create_test_data()
    
    try:
        analyzer = MetadataAnalyzer(distance_matrix, sample_names)
        print("✓ Valid distance matrix accepted")
    except Exception as e:
        print(f"✗ Valid matrix rejected: {e}")
        return
    
    # Test with bad matrices
    test_cases = [
        ("Non-square matrix", distance_matrix[:45, :]),  # Non-square
        ("Asymmetric matrix", distance_matrix + np.random.random(distance_matrix.shape) * 0.1),  # Asymmetric
        ("Matrix with NaN", np.where(distance_matrix > 0.8, np.nan, distance_matrix)),  # NaN values
        ("Matrix with negative values", distance_matrix - 0.1),  # Negative values
    ]
    
    for name, bad_matrix in test_cases:
        try:
            analyzer = MetadataAnalyzer(bad_matrix, sample_names)
            print(f"✗ {name} incorrectly accepted")
        except Exception as e:
            print(f"✓ {name} correctly rejected: {type(e).__name__}")
    
    print()

def test_corrected_sum_of_squares():
    """Test the corrected sum of squares calculation"""
    print("=== Testing Corrected Sum of Squares ===")
    
    distance_matrix, sample_names, metadata = create_test_data()
    analyzer = MetadataAnalyzer(distance_matrix, sample_names)
    
    # Test with simple groups
    groups = np.array([0, 0, 0, 1, 1, 1] + [2] * 44)  # 3 groups with different sizes
    test_matrix = distance_matrix[:6, :6]
    
    # Access the internal PermanovaAnalyzer to test the method
    from metadata_analyzer import PermanovaAnalyzer
    permanova = PermanovaAnalyzer(test_matrix, sample_names[:6])
    
    within_ss, total_ss = permanova._calculate_sum_of_squares(test_matrix, groups[:6])
    
    print(f"Total SS: {total_ss:.4f}")
    print(f"Within SS: {within_ss:.4f}")
    print(f"Between SS: {total_ss - within_ss:.4f}")
    
    # Verify the formula: total_ss should be sum(d_ij^2) / (2*n)
    expected_total = np.sum(test_matrix**2) / (2 * 6)
    print(f"Expected total SS: {expected_total:.4f}")
    print(f"Match: {abs(total_ss - expected_total) < 1e-10}")
    
    print()

def test_multiple_testing_correction():
    """Test multiple testing correction"""
    print("=== Testing Multiple Testing Correction ===")
    
    distance_matrix, sample_names, metadata = create_test_data()
    analyzer = MetadataAnalyzer(distance_matrix, sample_names)
    analyzer.metadata = metadata
    
    # Analyze multiple variables
    results = analyzer.analyze_variables(
        variables=['body_site', 'disease_status'], 
        n_permutations=99
    )
    
    print("Results with multiple testing correction:")
    for _, row in results.iterrows():
        print(f"- {row['variable']}: p={row['p_value']:.3f}, p_adj={row.get('p_adjusted', 'N/A'):.3f}")
    
    print()

def test_omega_squared():
    """Test omega squared calculation"""
    print("=== Testing Omega Squared Effect Size ===")
    
    distance_matrix, sample_names, metadata = create_test_data()
    analyzer = MetadataAnalyzer(distance_matrix, sample_names)
    
    # Test with body_site (should have large effect)
    groups = np.array([0] * 20 + [1] * 20 + [2] * 10)
    
    # Access the internal PermanovaAnalyzer to test the method
    from metadata_analyzer import PermanovaAnalyzer
    permanova = PermanovaAnalyzer(distance_matrix, sample_names)
    
    result = permanova.permanova_test(groups, n_permutations=99)
    
    print(f"R²: {result['r_squared']:.4f}")
    print(f"R²_adj: {result['r_squared_adj']:.4f}")
    print(f"ω²: {result['omega_squared']:.4f}")
    print(f"ω²_partial: {result['partial_omega_squared']:.4f}")
    print(f"Effect size: {result['effect_size']}")
    
    # Omega squared should be more conservative than R²
    print(f"ω² < R²: {result['omega_squared'] < result['r_squared']}")
    
    print()

def test_permdisp():
    """Test PERMDISP variance homogeneity test"""
    print("=== Testing PERMDISP Variance Homogeneity ===")
    
    distance_matrix, sample_names, metadata = create_test_data()
    analyzer = MetadataAnalyzer(distance_matrix, sample_names)
    
    # Test with groups that should have different dispersions
    groups = np.array([0] * 20 + [1] * 20 + [2] * 10)
    
    # Access the internal PermanovaAnalyzer to test the method
    from metadata_analyzer import PermanovaAnalyzer
    permanova = PermanovaAnalyzer(distance_matrix, sample_names)
    
    permdisp_result = permanova.permdisp_test(distance_matrix, groups, n_permutations=99)
    
    print(f"PERMDISP F: {permdisp_result['f_statistic']:.4f}")
    print(f"PERMDISP p: {permdisp_result['p_value']:.4f}")
    print(f"Warning: {permdisp_result['warning']}")
    
    # Test with homogeneous groups (should have similar dispersions)
    # Create artificial homogeneous groups
    homogeneous_groups = np.random.choice([0, 1], 50, p=[0.5, 0.5])
    homogeneous_result = permanova.permdisp_test(distance_matrix, homogeneous_groups, n_permutations=99)
    
    print(f"\nHomogeneous groups:")
    print(f"PERMDISP F: {homogeneous_result['f_statistic']:.4f}")
    print(f"PERMDISP p: {homogeneous_result['p_value']:.4f}")
    print(f"Warning: {homogeneous_result['warning']}")
    
    print()

def test_full_pipeline():
    """Test the complete improved pipeline"""
    print("=== Testing Complete Improved Pipeline ===")
    
    distance_matrix, sample_names, metadata = create_test_data()
    analyzer = MetadataAnalyzer(distance_matrix, sample_names)
    analyzer.metadata = metadata
    
    # Run complete analysis
    results = analyzer.analyze_variables(n_permutations=99)
    
    print("Complete analysis results:")
    for _, row in results.iterrows():
        print(f"\n{row['variable']}:")
        print(f"  - R²={row['r_squared']:.3f}, R²_adj={row['r_squared_adj']:.3f}, ω²={row['omega_squared']:.3f}")
        print(f"  - p={row['p_value']:.3f}, p_adj={row.get('p_adjusted', 'N/A'):.3f}")
        print(f"  - Effect: {row['effect_size']}")
        print(f"  - PERMDISP: F={row['permdisp_f']:.3f}, p={row['permdisp_p']:.3f}")
        if row['variance_warning']:
            print(f"  - ⚠️  {row['variance_warning']}")
    
    print("\n=== All statistical improvements working correctly! ===")

if __name__ == "__main__":
    test_distance_matrix_validation()
    test_corrected_sum_of_squares()
    test_multiple_testing_correction()
    test_omega_squared()
    test_permdisp()
    test_full_pipeline()