#!/usr/bin/env python3
"""
Test script for parallel processing performance improvements.

This script creates multiple test FASTQ files and benchmarks sequential vs parallel processing.
"""

import tempfile
import time
import logging
from pathlib import Path
import shutil
import sys
sys.path.insert(0, '.')

from metagrouper import KmerProfiler, find_fastq_files


def create_test_dataset(num_samples: int, reads_per_sample: int = 50, temp_dir: Path = None):
    """Create a test dataset with multiple samples."""
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp())
    
    sequences = [
        "ATCGATCGATCGATCGATCGATCGATCGATCG",
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA", 
        "CGTACGTACGTACGTACGTACGTACGTACGTA",
        "TATGTATGTATGTATGTATGTATGTATGTATG",
        "AAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTT",
        "TTCCGGAATTCCGGAATTCCGGAATTCCGGAA"
    ]
    
    for i in range(num_samples):
        # Create paired-end files for some samples, single-end for others
        if i % 3 == 0:  # Single-end
            filename = f"sample_{i:03d}.fastq"
            filepath = temp_dir / filename
            with open(filepath, "w") as f:
                for j in range(reads_per_sample):
                    seq = sequences[j % len(sequences)]
                    f.write(f"@sample_{i}_read_{j}\n")
                    f.write(f"{seq}\n")
                    f.write("+\n")
                    f.write("~" * len(seq) + "\n")
        else:  # Paired-end
            for direction in ["R1", "R2"]:
                filename = f"sample_{i:03d}_{direction}.fastq"
                filepath = temp_dir / filename
                with open(filepath, "w") as f:
                    for j in range(reads_per_sample):
                        seq = sequences[j % len(sequences)]
                        f.write(f"@sample_{i}_{direction}_read_{j}\n")
                        f.write(f"{seq}\n")
                        f.write("+\n")
                        f.write("~" * len(seq) + "\n")
    
    return temp_dir


def benchmark_processing(test_dir: Path, k_size: int = 17, max_reads: int = None):
    """Benchmark sequential vs parallel processing."""
    print(f"Benchmarking k-mer profiling in {test_dir}")
    
    # Find FASTQ files
    fastq_files = find_fastq_files(str(test_dir))
    print(f"Found {len(fastq_files)} samples")
    
    # Test sequential processing
    print("\n=== Sequential Processing ===")
    profiler_seq = KmerProfiler(k=k_size, max_reads=max_reads)
    
    start_time = time.time()
    seq_failed = []
    
    for filepath, sample_name in fastq_files:
        try:
            profiler_seq.profile_sample(filepath, sample_name)
        except Exception as e:
            seq_failed.append(sample_name)
            print(f"Sequential failed: {sample_name}: {e}")
    
    seq_time = time.time() - start_time
    seq_success = len(fastq_files) - len(seq_failed)
    
    print(f"Sequential: {seq_success}/{len(fastq_files)} samples in {seq_time:.2f}s")
    print(f"Rate: {seq_success/seq_time:.1f} samples/second")
    
    # Test parallel processing
    print("\n=== Parallel Processing ===")
    profiler_par = KmerProfiler(k=k_size, max_reads=max_reads)
    
    start_time = time.time()
    
    try:
        profiles, par_failed = profiler_par.process_samples_parallel(
            fastq_files, 
            n_processes=None,  # Use all available cores
            show_progress=True
        )
        par_time = time.time() - start_time
        par_success = len(fastq_files) - len(par_failed)
        
        print(f"Parallel: {par_success}/{len(fastq_files)} samples in {par_time:.2f}s")
        print(f"Rate: {par_success/par_time:.1f} samples/second")
        
        # Calculate speedup
        if seq_time > 0:
            speedup = seq_time / par_time
            print(f"\nSpeedup: {speedup:.1f}x faster")
        
        # Verify results are similar
        if seq_success == par_success:
            print("✓ Same number of successful samples")
        else:
            print(f"⚠ Different success rates: seq={seq_success}, par={par_success}")
            
    except Exception as e:
        print(f"Parallel processing failed: {e}")
        return seq_time, None, None
    
    return seq_time, par_time, seq_success


def main():
    """Run performance benchmark."""
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    
    print("=== MetaGrouper Parallel Processing Benchmark ===\n")
    
    # Test with different dataset sizes
    test_cases = [
        (5, 20, "Small dataset"),
        (10, 30, "Medium dataset"), 
        (20, 25, "Large dataset")
    ]
    
    for num_samples, reads_per_sample, description in test_cases:
        print(f"\n{description}: {num_samples} samples, {reads_per_sample} reads each")
        print("=" * 60)
        
        # Create test dataset
        test_dir = create_test_dataset(num_samples, reads_per_sample)
        
        try:
            # Run benchmark
            seq_time, par_time, success_count = benchmark_processing(
                test_dir, 
                k_size=15,  # Smaller k for faster testing
                max_reads=reads_per_sample
            )
            
            if par_time:
                efficiency = (seq_time / par_time) / 4  # Assuming 4 cores
                print(f"Parallel efficiency: {efficiency:.1f} (1.0 = perfect scaling)")
        
        finally:
            # Cleanup
            shutil.rmtree(test_dir)
    
    print("\n=== Benchmark Complete ===")
    print("Tips for optimal performance:")
    print("• Use parallel processing for >2 samples")
    print("• Adjust --processes based on CPU cores and memory")
    print("• Use --sequential for debugging or small datasets")


if __name__ == "__main__":
    main()