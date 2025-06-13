#!/usr/bin/env python3
"""
Quick test to check if our parallelization is working.
"""

import sys
import time
import multiprocessing as mp
from pathlib import Path

# Add the package to path
sys.path.insert(0, 'metagrouper_package')

from metagrouper.profiler import KmerProfiler
from metagrouper.utils import find_fastq_files

def test_parallel_vs_sequential(data_dir, max_samples=10):
    """Test if parallel processing is actually faster."""
    
    print("🔍 Testing MetaGrouper parallelization...")
    
    # Find FASTQ files
    fastq_files = find_fastq_files(data_dir)
    
    if len(fastq_files) == 0:
        print(f"❌ No FASTQ files found in {data_dir}")
        return
    
    # Limit to first few samples for testing
    test_files = fastq_files[:min(max_samples, len(fastq_files))]
    print(f"📁 Found {len(fastq_files)} total files, testing with {len(test_files)} samples")
    
    # Test sequential processing
    print("\n🔄 Testing sequential processing...")
    profiler_seq = KmerProfiler(k=15, max_reads=1000)  # Small for speed
    
    start_time = time.time()
    for filepath, sample_name in test_files:
        profiler_seq.profile_sample(filepath, sample_name)
    sequential_time = time.time() - start_time
    
    print(f"⏱️  Sequential: {sequential_time:.1f}s ({sequential_time/len(test_files):.1f}s per sample)")
    
    # Test parallel processing
    print("\n⚡ Testing parallel processing...")
    profiler_par = KmerProfiler(k=15, max_reads=1000)
    
    start_time = time.time()
    try:
        profiles, failed = profiler_par.process_samples_parallel(
            test_files, 
            n_processes=min(4, mp.cpu_count(), len(test_files)),
            show_progress=False
        )
        parallel_time = time.time() - start_time
        
        print(f"⏱️  Parallel: {parallel_time:.1f}s ({parallel_time/len(test_files):.1f}s per sample)")
        print(f"🚀 Speedup: {sequential_time/parallel_time:.1f}x")
        
        if parallel_time < sequential_time:
            print("✅ Parallel processing is working!")
        else:
            print("⚠️  Parallel processing is slower (normal for small datasets)")
            
    except Exception as e:
        print(f"❌ Parallel processing failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python debug_parallel.py <fastq_directory>")
        sys.exit(1)
    
    test_parallel_vs_sequential(sys.argv[1])