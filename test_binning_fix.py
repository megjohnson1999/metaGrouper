#!/usr/bin/env python3
"""
Test the fixed binning logic for continuous variables
"""

import sys
import numpy as np
import pandas as pd
sys.path.append('phases')

from metadata_analyzer import MetadataAnalyzer

def test_binning_fix():
    """Test that the binning assignment works correctly"""
    
    # Create test data similar to what would cause the error
    np.random.seed(42)
    n_samples = 1150  # Same size as the failing case
    
    # Create distance matrix
    distance_matrix = np.random.random((n_samples, n_samples))
    distance_matrix = (distance_matrix + distance_matrix.T) / 2
    np.fill_diagonal(distance_matrix, 0)
    
    sample_names = [str(i) for i in range(n_samples)]
    
    # Create metadata with continuous variable (like finaImputedlLogFCP)
    metadata = pd.DataFrame(index=sample_names)
    
    # Create continuous variable with some missing values (like real data)
    continuous_values = np.random.lognormal(4, 1, n_samples)
    # Set some values to NaN to simulate missing data
    missing_indices = np.random.choice(n_samples, size=86, replace=False)  # 1150 - 1064 = 86 missing
    continuous_values[missing_indices] = np.nan
    metadata['test_continuous'] = continuous_values
    
    print(f"Created test data:")
    print(f"  Total samples: {n_samples}")
    print(f"  Valid values: {np.sum(~np.isnan(continuous_values))}")
    print(f"  Missing values: {np.sum(np.isnan(continuous_values))}")
    
    # Test the analysis
    analyzer = MetadataAnalyzer(distance_matrix, sample_names)
    analyzer.metadata = metadata
    
    try:
        results = analyzer.analyze_variables(
            variables=['test_continuous'],
            n_permutations=99
        )
        print("✅ Binning fix successful!")
        print(f"Results: {len(results)} variables analyzed")
        
        for _, row in results.iterrows():
            print(f"  {row['variable']}: R²={row['r_squared']:.3f}, groups={row['n_groups']}")
            
    except Exception as e:
        print(f"❌ Error still occurs: {e}")
        return False
    
    return True

if __name__ == "__main__":
    test_binning_fix()