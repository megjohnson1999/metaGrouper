#!/usr/bin/env python3
"""
Test script for the new adaptive binning functionality
"""

import sys
import numpy as np
import pandas as pd
sys.path.append('phases')

from metadata_analyzer import MetadataAnalyzer

def test_adaptive_binning():
    """Test the adaptive binning function"""
    
    # Create a dummy MetadataAnalyzer instance
    dummy_distance_matrix = np.random.random((10, 10))
    dummy_sample_names = [f"sample_{i}" for i in range(10)]
    analyzer = MetadataAnalyzer(dummy_distance_matrix, dummy_sample_names)
    
    # Test case 1: Normal distribution with 1000 samples
    print("Test 1: Normal distribution (1000 samples)")
    np.random.seed(42)
    values = np.random.normal(100, 20, 1000)
    
    bins = analyzer._adaptive_binning(values, min_bin_size=30, max_bins=5)
    unique_bins = np.unique(bins)
    
    print(f"  - Number of bins: {len(unique_bins)}")
    for i, bin_val in enumerate(unique_bins):
        count = np.sum(bins == bin_val)
        bin_mean = np.mean(values[bins == bin_val])
        print(f"  - Bin {i}: {count} samples, mean = {bin_mean:.2f}")
    
    print()
    
    # Test case 2: Skewed distribution (like log-FCP)
    print("Test 2: Log-normal distribution (simulating logFCP)")
    np.random.seed(42)
    values = np.random.lognormal(4, 1, 1000)  # Skewed distribution
    
    bins = analyzer._adaptive_binning(values, min_bin_size=30, max_bins=5)
    unique_bins = np.unique(bins)
    
    print(f"  - Number of bins: {len(unique_bins)}")
    for i, bin_val in enumerate(unique_bins):
        count = np.sum(bins == bin_val)
        bin_mean = np.mean(values[bins == bin_val])
        bin_range = f"{np.min(values[bins == bin_val]):.2f}-{np.max(values[bins == bin_val]):.2f}"
        print(f"  - Bin {i}: {count} samples, mean = {bin_mean:.2f}, range = {bin_range}")
    
    print()
    
    # Test case 3: Small dataset
    print("Test 3: Small dataset (100 samples)")
    np.random.seed(42)
    values = np.random.normal(100, 20, 100)
    
    bins = analyzer._adaptive_binning(values, min_bin_size=30, max_bins=5)
    unique_bins = np.unique(bins)
    
    print(f"  - Number of bins: {len(unique_bins)}")
    for i, bin_val in enumerate(unique_bins):
        count = np.sum(bins == bin_val)
        bin_mean = np.mean(values[bins == bin_val])
        print(f"  - Bin {i}: {count} samples, mean = {bin_mean:.2f}")
    
    print()
    
    # Test case 4: Very small dataset (should return single bin)
    print("Test 4: Very small dataset (20 samples)")
    np.random.seed(42)
    values = np.random.normal(100, 20, 20)
    
    bins = analyzer._adaptive_binning(values, min_bin_size=30, max_bins=5)
    unique_bins = np.unique(bins)
    
    print(f"  - Number of bins: {len(unique_bins)}")
    for i, bin_val in enumerate(unique_bins):
        count = np.sum(bins == bin_val)
        bin_mean = np.mean(values[bins == bin_val])
        print(f"  - Bin {i}: {count} samples, mean = {bin_mean:.2f}")
    
    print("All tests completed successfully!")

if __name__ == "__main__":
    test_adaptive_binning()