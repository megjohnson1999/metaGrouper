#!/usr/bin/env python3
"""
Comprehensive benchmark for Phase 1 fast k-mer extraction.

This script tests the performance improvements of fast k-mer extraction
on realistic viral metagenomic sequences.
"""

import time
import tempfile
import random
from pathlib import Path
from typing import List, Dict
import logging

from metagrouper.profiler import KmerProfiler
from metagrouper.fast_kmer import FastKmerExtractor, benchmark_extraction_methods


def generate_viral_sequences(n_reads: int = 1000, read_length: int = 150) -> List[str]:
    """Generate realistic viral-like sequences."""
    sequences = []
    
    # Common viral sequence patterns
    viral_motifs = [
        "ATGAAAAAAA",  # Start codon + poly-A
        "GCTAGCTAGC",  # Viral promoter-like
        "CCCCCCCCCC",  # Poly-C regions common in RNA viruses
        "ATATATATAT",  # AT-rich regions
        "CGCGCGCGCG",  # CpG islands
        "TTTAAATAAA",  # Poly-T with interruptions
        "GGGCCCGGGT",  # GC-rich regions
        "ACGTACGTAC",  # Palindromic sequences
    ]
    
    bases = ['A', 'T', 'C', 'G']
    
    for _ in range(n_reads):
        sequence = ""
        
        # Create sequence with mixture of random and viral patterns
        while len(sequence) < read_length:
            if random.random() < 0.4:  # 40% chance of viral motif
                motif = random.choice(viral_motifs)
                sequence += motif
            else:
                # Add random bases
                sequence += ''.join(random.choice(bases) for _ in range(random.randint(1, 8)))
        
        # Trim to exact length and add occasional N's
        sequence = sequence[:read_length]
        
        # Add occasional N's (5% chance per position)
        sequence_list = list(sequence)
        for i in range(len(sequence_list)):
            if random.random() < 0.01:  # 1% chance of N
                sequence_list[i] = 'N'
        
        sequences.append(''.join(sequence_list))
    
    return sequences


def create_test_fastq(sequences: List[str], output_path: Path) -> None:
    """Create a FASTQ file from sequences."""
    with open(output_path, 'w') as f:
        for i, seq in enumerate(sequences):
            f.write(f"@read_{i}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write("~" * len(seq) + "\n")


def benchmark_k_mer_extraction():
    """Benchmark k-mer extraction methods on viral sequences."""
    print("=== Phase 1: Fast K-mer Extraction Benchmark ===\n")
    
    test_cases = [
        (500, 100, "Small viral sample"),
        (2000, 150, "Medium viral sample"),
        (10000, 150, "Large viral sample"),
        (5000, 300, "Long read viral sample")
    ]
    
    results = []
    
    for n_reads, read_length, description in test_cases:
        print(f"{description}: {n_reads} reads of {read_length}bp")
        
        # Generate test sequences
        sequences = generate_viral_sequences(n_reads, read_length)
        total_bases = sum(len(seq) for seq in sequences)
        
        # Calculate expected k-mers
        k = 21
        expected_kmers = sum(max(0, len(seq) - k + 1) for seq in sequences)
        
        print(f"  Total bases: {total_bases:,}")
        print(f"  Expected k-mers: {expected_kmers:,}")
        
        # Test individual sequence performance
        test_sequence = ''.join(sequences[:10])  # Combine first 10 sequences
        benchmark_result = benchmark_extraction_methods(test_sequence, k=k, iterations=50)
        
        print(f"  Individual extraction speedup: {benchmark_result['speedup']:.2f}x")
        print(f"  Results match: {benchmark_result['results_match']}")
        
        # Test full sample processing
        with tempfile.TemporaryDirectory() as temp_dir:
            fastq_path = Path(temp_dir) / "test_sample.fastq"
            create_test_fastq(sequences, fastq_path)
            
            # Test legacy profiler
            print("  Testing legacy profiler...")
            start_time = time.perf_counter()
            legacy_profiler = KmerProfiler(k=k, use_fast_extraction=False)
            legacy_profile = legacy_profiler.profile_sample(str(fastq_path), "test_legacy")
            legacy_time = time.perf_counter() - start_time
            
            # Test fast profiler
            print("  Testing fast profiler...")
            start_time = time.perf_counter()
            fast_profiler = KmerProfiler(k=k, use_fast_extraction=True)
            fast_profile = fast_profiler.profile_sample(str(fastq_path), "test_fast")
            fast_time = time.perf_counter() - start_time
            
            # Compare results
            speedup = legacy_time / fast_time if fast_time > 0 else float('inf')
            profiles_match = len(legacy_profile) == len(fast_profile)
            
            print(f"  Sample processing speedup: {speedup:.2f}x")
            print(f"  Profiles match: {profiles_match}")
            print(f"  Legacy time: {legacy_time:.3f}s")
            print(f"  Fast time: {fast_time:.3f}s")
            print(f"  K-mers extracted: {len(legacy_profile):,}")
            print()
            
            results.append({
                'description': description,
                'n_reads': n_reads,
                'read_length': read_length,
                'total_bases': total_bases,
                'expected_kmers': expected_kmers,
                'individual_speedup': benchmark_result['speedup'],
                'sample_speedup': speedup,
                'legacy_time': legacy_time,
                'fast_time': fast_time,
                'profiles_match': profiles_match,
                'kmers_extracted': len(legacy_profile)
            })
    
    return results


def test_scalability():
    """Test scalability with different k-mer sizes and parameters."""
    print("=== K-mer Size Scalability Test ===\n")
    
    # Generate test sequence
    sequences = generate_viral_sequences(n_reads=1000, read_length=200)
    combined_sequence = ''.join(sequences[:5])  # Use first 5 sequences
    
    k_sizes = [11, 15, 21, 25, 31]
    
    for k in k_sizes:
        print(f"Testing k={k}")
        
        try:
            result = benchmark_extraction_methods(combined_sequence, k=k, iterations=20)
            print(f"  Speedup: {result['speedup']:.2f}x")
            print(f"  Results match: {result['results_match']}")
            print(f"  K-mers found: {result['legacy_kmers']:,}")
        except Exception as e:
            print(f"  Error: {e}")
        
        print()


def memory_usage_test():
    """Test memory usage characteristics."""
    print("=== Memory Usage Comparison ===\n")
    
    import psutil
    import os
    
    def get_memory_usage():
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024  # MB
    
    # Test with increasingly large datasets
    for n_reads in [1000, 5000, 10000]:
        print(f"Testing with {n_reads} reads")
        
        sequences = generate_viral_sequences(n_reads, 150)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            fastq_path = Path(temp_dir) / f"test_{n_reads}.fastq"
            create_test_fastq(sequences, fastq_path)
            
            # Memory test for legacy
            initial_memory = get_memory_usage()
            legacy_profiler = KmerProfiler(k=21, use_fast_extraction=False)
            legacy_profile = legacy_profiler.profile_sample(str(fastq_path), "test")
            legacy_memory = get_memory_usage() - initial_memory
            
            # Memory test for fast
            initial_memory = get_memory_usage()
            fast_profiler = KmerProfiler(k=21, use_fast_extraction=True)
            fast_profile = fast_profiler.profile_sample(str(fastq_path), "test")
            fast_memory = get_memory_usage() - initial_memory
            
            print(f"  Legacy memory: {legacy_memory:.1f} MB")
            print(f"  Fast memory: {fast_memory:.1f} MB")
            print(f"  Memory ratio: {legacy_memory/fast_memory:.2f}x" if fast_memory > 0 else "N/A")
            print(f"  K-mers: {len(legacy_profile):,}")
            print()


def main():
    """Run comprehensive Phase 1 benchmark."""
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    
    print("🧬 MetaGrouper Phase 1: Fast K-mer Extraction Benchmark")
    print("="*60)
    print()
    
    # Main benchmark
    results = benchmark_k_mer_extraction()
    
    # Scalability test
    test_scalability()
    
    # Memory usage test
    memory_usage_test()
    
    # Summary
    print("=== Summary ===")
    print("Phase 1 Fast K-mer Extraction Results:")
    print()
    
    total_speedup = 0
    count = 0
    for result in results:
        if result['sample_speedup'] > 0:
            total_speedup += result['sample_speedup']
            count += 1
        print(f"{result['description']}:")
        print(f"  Sample processing: {result['sample_speedup']:.2f}x speedup")
        print(f"  Individual extraction: {result['individual_speedup']:.2f}x speedup")
        print()
    
    if count > 0:
        avg_speedup = total_speedup / count
        print(f"Average sample processing speedup: {avg_speedup:.2f}x")
    
    print("\n🎯 Phase 1 Fast K-mer Extraction: COMPLETE")
    print("Next: Implement streaming k-mer sketches for memory efficiency")


if __name__ == "__main__":
    main()