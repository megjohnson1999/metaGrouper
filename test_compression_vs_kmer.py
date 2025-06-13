#!/usr/bin/env python3
"""
Test script to compare compression-based vs k-mer based approaches.
"""

import sys
import time
import tempfile
import shutil
from pathlib import Path

# Add package to path
sys.path.insert(0, 'metagrouper_package')

from metagrouper.compression_profiler import CompressionProfiler
from metagrouper.compression_analyzer import CompressionAnalyzer
from metagrouper.profiler import KmerProfiler
from metagrouper.analyzer import SimilarityAnalyzer
from metagrouper.utils import find_fastq_files


def create_test_data(temp_dir, n_samples=5):
    """Create synthetic FASTQ data for testing."""
    
    # Different sequence patterns to create distinct samples
    patterns = [
        "ATCGATCGATCGATCGATCG",  # Pattern A - AT rich
        "GCTAGCTAGCTAGCTAGCTA",  # Pattern B - GC rich  
        "TTCCGGAATTCCGGAATTCC",  # Pattern C - repeated motifs
        "AAATTTGGGCCCAAATTTGGG", # Pattern D - longer repeats
        "ACTGACTGACTGACTGACTG",  # Pattern E - simple repeat
    ]
    
    fastq_files = []
    
    for i in range(min(n_samples, len(patterns))):
        sample_name = f"sample_{i+1:03d}"
        fastq_file = temp_dir / f"{sample_name}.fastq"
        
        pattern = patterns[i]
        
        with open(fastq_file, 'w') as f:
            # Create 100 reads per sample, each with variations of the pattern
            for read_num in range(100):
                # Add some variation to the pattern
                seq = pattern * 3  # Make longer sequences
                if read_num % 10 == 0:  # Add some noise every 10th read
                    seq = seq.replace('A', 'T', 2)
                
                f.write(f"@read_{read_num}\n")
                f.write(f"{seq}\n")
                f.write("+\n")
                f.write("~" * len(seq) + "\n")
        
        fastq_files.append((str(fastq_file), sample_name))
    
    return fastq_files


def test_compression_approach(fastq_files):
    """Test compression-based approach."""
    print("🗜️  Testing compression-based approach...")
    
    start_time = time.time()
    
    # Initialize compression profiler
    profiler = CompressionProfiler(
        compressor='zstd',
        max_reads=None,  # Use all reads for fair comparison
        compression_level=3
    )
    
    # Process samples
    for filepath, sample_name in fastq_files:
        profiler.profile_sample(filepath, sample_name)
    
    # Compute distance matrix
    distance_matrix = profiler.compute_distance_matrix()
    
    # Analyze
    analyzer = CompressionAnalyzer(distance_matrix, profiler.sample_names)
    stats = analyzer.compute_summary_statistics()
    
    elapsed_time = time.time() - start_time
    
    print(f"   ⏱️  Time: {elapsed_time:.2f}s")
    print(f"   📊 Distance stats: mean={stats['mean_distance']:.3f}, std={stats['std_distance']:.3f}")
    print(f"   💾 Compression sizes: {list(profiler.compression_sizes.values())}")
    
    return elapsed_time, distance_matrix, stats


def test_kmer_approach(fastq_files):
    """Test k-mer based approach."""
    print("🧬 Testing k-mer based approach...")
    
    start_time = time.time()
    
    # Initialize k-mer profiler
    profiler = KmerProfiler(
        k=15,  # Smaller k for speed
        max_reads=None,  # Use all reads for fair comparison
        min_kmer_freq=1
    )
    
    # Process samples
    for filepath, sample_name in fastq_files:
        profiler.profile_sample(filepath, sample_name)
    
    # Analyze
    analyzer = SimilarityAnalyzer(profiler.profiles)
    distance_matrix = analyzer.compute_distance_matrix()
    
    elapsed_time = time.time() - start_time
    
    # Compute stats
    import numpy as np
    upper_triangle = distance_matrix[np.triu_indices_from(distance_matrix, k=1)]
    mean_dist = np.mean(upper_triangle)
    std_dist = np.std(upper_triangle)
    
    print(f"   ⏱️  Time: {elapsed_time:.2f}s")
    print(f"   📊 Distance stats: mean={mean_dist:.3f}, std={std_dist:.3f}")
    print(f"   🧬 K-mer profile sizes: {[len(profile) for profile in profiler.profiles.values()]}")
    
    return elapsed_time, distance_matrix


def compare_distance_matrices(compression_matrix, kmer_matrix, sample_names):
    """Compare the two distance matrices."""
    print("\n🔍 Comparing distance matrices...")
    
    import numpy as np
    from scipy.stats import pearsonr
    
    # Extract upper triangles for comparison
    comp_upper = compression_matrix[np.triu_indices_from(compression_matrix, k=1)]
    kmer_upper = kmer_matrix[np.triu_indices_from(kmer_matrix, k=1)]
    
    # Compute correlation
    correlation, p_value = pearsonr(comp_upper, kmer_upper)
    
    print(f"   📈 Correlation between methods: {correlation:.3f} (p={p_value:.3e})")
    
    # Show pairwise comparisons
    print(f"   📋 Sample-by-sample comparison:")
    print(f"       {'Pair':<20} {'Compression':<12} {'K-mer':<12} {'Diff':<10}")
    print(f"       {'-'*20} {'-'*12} {'-'*12} {'-'*10}")
    
    idx = 0
    for i in range(len(sample_names)):
        for j in range(i+1, len(sample_names)):
            pair_name = f"{sample_names[i][:8]}-{sample_names[j][:8]}"
            comp_dist = compression_matrix[i, j]
            kmer_dist = kmer_matrix[i, j]
            diff = abs(comp_dist - kmer_dist)
            
            print(f"       {pair_name:<20} {comp_dist:<12.3f} {kmer_dist:<12.3f} {diff:<10.3f}")
            idx += 1
            if idx >= 5:  # Limit output
                break
        if idx >= 5:
            break
    
    return correlation


def main():
    """Run comparison test."""
    print("🧪 MetaGrouper: Compression vs K-mer Approach Comparison")
    print("=" * 60)
    
    # Create test data
    temp_dir = Path(tempfile.mkdtemp())
    print(f"📁 Creating test data in {temp_dir}")
    
    try:
        fastq_files = create_test_data(temp_dir, n_samples=5)
        print(f"   Created {len(fastq_files)} test samples")
        
        # Test both approaches
        print(f"\n🏁 Running performance comparison...")
        
        comp_time, comp_matrix, comp_stats = test_compression_approach(fastq_files)
        kmer_time, kmer_matrix = test_kmer_approach(fastq_files)
        
        # Compare results
        sample_names = [name for _, name in fastq_files]
        correlation = compare_distance_matrices(comp_matrix, kmer_matrix, sample_names)
        
        # Summary
        print(f"\n📊 Performance Summary:")
        print(f"   🗜️  Compression approach: {comp_time:.2f}s")
        print(f"   🧬 K-mer approach: {kmer_time:.2f}s")
        print(f"   🚀 Speedup: {kmer_time/comp_time:.1f}x")
        print(f"   📈 Correlation: {correlation:.3f}")
        
        if comp_time < kmer_time and correlation > 0.7:
            print(f"   ✅ Compression approach is faster AND well-correlated!")
        elif comp_time < kmer_time:
            print(f"   ⚡ Compression approach is faster but correlation is low")
        else:
            print(f"   ⚠️  K-mer approach is still faster for this test size")
        
        print(f"\n💡 For your 340 samples:")
        estimated_comp_time = comp_time * (340 / len(fastq_files))
        estimated_kmer_time = kmer_time * (340 / len(fastq_files))
        print(f"   🗜️  Estimated compression time: {estimated_comp_time/60:.1f} minutes")
        print(f"   🧬 Estimated k-mer time: {estimated_kmer_time/60:.1f} minutes")
        
    finally:
        # Cleanup
        shutil.rmtree(temp_dir)
        print(f"\n🧹 Cleaned up test data")


if __name__ == "__main__":
    main()