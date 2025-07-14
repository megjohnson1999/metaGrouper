#!/usr/bin/env python3
"""
Test forcing more bins to see if meaningful groups emerge
"""

import sys
import numpy as np
import pandas as pd
sys.path.append('phases')

def test_forced_bins():
    """Test forcing 5, 6, 7 bins to see the results"""
    
    # Load actual FCP data
    metadata = pd.read_csv('/home/megan_johnson/wustl/metaGrouper/21032025_Metadata_freeze4_imputed_fin.txt', 
                          sep='\t', index_col=0)
    
    fcp_values = metadata['finaImputedlLogFCP'].dropna().values
    print(f"FCP data: {len(fcp_values)} samples")
    print(f"Range: {fcp_values.min():.2f} to {fcp_values.max():.2f}")
    
    # Test different forced binning strategies
    for n_bins in [3, 4, 5, 6, 7, 8]:
        print(f"\n=== Forced {n_bins} bins (equal size) ===")
        
        # Use quantiles to create n bins
        percentiles = np.linspace(0, 100, n_bins + 1)
        bin_edges = np.percentile(fcp_values, percentiles)
        
        # Create bins
        bins = np.digitize(fcp_values, bin_edges[1:-1])
        
        unique_bins = np.unique(bins)
        print(f"Result: {len(unique_bins)} bins")
        
        for i, bin_val in enumerate(unique_bins):
            mask = bins == bin_val
            count = np.sum(mask)
            bin_mean = np.mean(fcp_values[mask])
            bin_min = np.min(fcp_values[mask])
            bin_max = np.max(fcp_values[mask])
            bin_std = np.std(fcp_values[mask])
            print(f"  Bin {i}: {count:3d} samples, mean={bin_mean:5.2f}±{bin_std:4.2f}, range=[{bin_min:5.2f}, {bin_max:5.2f}]")
        
        # Calculate between-group variance for comparison
        overall_mean = np.mean(fcp_values)
        between_var = 0
        for bin_val in unique_bins:
            mask = bins == bin_val
            bin_mean = np.mean(fcp_values[mask])
            bin_size = np.sum(mask)
            between_var += bin_size * (bin_mean - overall_mean) ** 2
        
        between_var /= len(fcp_values)
        total_var = np.var(fcp_values)
        variance_explained = between_var / total_var
        print(f"  Variance explained: {variance_explained:.4f} ({variance_explained*100:.1f}%)")

if __name__ == "__main__":
    test_forced_bins()