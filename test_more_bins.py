#!/usr/bin/env python3
"""
Test more aggressive binning parameters for FCP
"""

import sys
import numpy as np
import pandas as pd
sys.path.append('phases')

from metadata_analyzer import MetadataAnalyzer

def test_more_bins():
    """Test more aggressive binning strategies"""
    
    # Load the actual metadata
    metadata = pd.read_csv('/home/megan_johnson/wustl/metaGrouper/21032025_Metadata_freeze4_imputed_fin.txt', 
                          sep='\t', index_col=0)
    
    fcp_values = metadata['finaImputedlLogFCP'].dropna().values
    print(f"FCP data: {len(fcp_values)} samples")
    
    # Create dummy analyzer
    dummy_distance = np.random.random((10, 10))
    analyzer = MetadataAnalyzer(dummy_distance, ['s1']*10)
    
    # Test more aggressive parameters
    aggressive_configs = [
        {"min_bin_size": 15, "max_bins": 10, "name": "More bins (15 min)"},
        {"min_bin_size": 10, "max_bins": 12, "name": "Many bins (10 min)"},
        {"min_bin_size": 25, "max_bins": 6, "name": "Moderate (25 min, 6 max)"},
        {"min_bin_size": 20, "max_bins": 8, "name": "New default (20 min, 8 max)"},
    ]
    
    for config in aggressive_configs:
        print(f"\n=== {config['name']} ===")
        
        bins = analyzer._adaptive_binning(
            fcp_values, 
            min_bin_size=config['min_bin_size'], 
            max_bins=config['max_bins']
        )
        
        unique_bins = np.unique(bins)
        print(f"Result: {len(unique_bins)} bins")
        
        for i, bin_val in enumerate(unique_bins):
            count = np.sum(bins == bin_val)
            bin_mean = np.mean(fcp_values[bins == bin_val])
            bin_min = np.min(fcp_values[bins == bin_val])
            bin_max = np.max(fcp_values[bins == bin_val])
            print(f"  Bin {i}: {count:3d} samples, mean={bin_mean:5.2f}, range=[{bin_min:5.2f}, {bin_max:5.2f}]")
    
    # Test clinical-relevant ranges
    print(f"\n=== Clinical Quartiles (for comparison) ===")
    quartiles = np.percentile(fcp_values, [25, 50, 75])
    print(f"Quartiles: Q1={quartiles[0]:.2f}, Q2={quartiles[1]:.2f}, Q3={quartiles[2]:.2f}")
    
    # Manual quartile binning
    quartile_bins = np.digitize(fcp_values, quartiles)
    unique_q_bins = np.unique(quartile_bins)
    print(f"Quartile binning: {len(unique_q_bins)} bins")
    
    for i, bin_val in enumerate(unique_q_bins):
        count = np.sum(quartile_bins == bin_val)
        bin_mean = np.mean(fcp_values[quartile_bins == bin_val])
        bin_min = np.min(fcp_values[quartile_bins == bin_val])
        bin_max = np.max(fcp_values[quartile_bins == bin_val])
        print(f"  Q-Bin {i}: {count:3d} samples, mean={bin_mean:5.2f}, range=[{bin_min:5.2f}, {bin_max:5.2f}]")

if __name__ == "__main__":
    test_more_bins()