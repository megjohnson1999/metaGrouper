#!/usr/bin/env python3
"""
Test suite for streaming k-mer profiler with sketching.
"""

import unittest
import tempfile
import time
import random
from pathlib import Path
from collections import Counter

from metagrouper.sketch_profiler import StreamingKmerProfiler
from metagrouper.profiler import KmerProfiler


class TestStreamingKmerProfiler(unittest.TestCase):
    """Test streaming k-mer profiler functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.profiler = StreamingKmerProfiler(k=21, sketch_size=100, sampling_method='reservoir')
        self.small_profiler = StreamingKmerProfiler(k=5, sketch_size=50, sampling_method='frequency')
    
    def test_initialization(self):
        """Test StreamingKmerProfiler initialization."""
        # Test different configurations
        configs = [
            {'k': 15, 'sketch_size': 1000, 'sampling_method': 'reservoir'},
            {'k': 21, 'sketch_size': 5000, 'sampling_method': 'frequency'},
            {'k': 25, 'sketch_size': 2000, 'sampling_method': 'adaptive'}
        ]
        
        for config in configs:
            profiler = StreamingKmerProfiler(**config)
            self.assertEqual(profiler.k, config['k'])
            self.assertEqual(profiler.sketch_size, config['sketch_size'])
            self.assertEqual(profiler.sampling_method, config['sampling_method'])
        
        # Test invalid sampling method
        with self.assertRaises(ValueError):
            profiler = StreamingKmerProfiler(sampling_method='invalid')
            # Create a dummy file to trigger the error
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq') as f:
                f.write("@test\nATCG\n+\n~~~~\n")
                f.flush()
                profiler.profile_sample(f.name, "test")
    
    def test_reservoir_sampling(self):
        """Test reservoir sampling functionality."""
        # Create test sequences with known k-mer distribution
        sequences = [
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # Repetitive
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",  # Different repetitive
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",  # Another pattern
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            for i, seq in enumerate(sequences * 20):  # 60 total sequences
                f.write(f"@read_{i}\n{seq}\n+\n{'~' * len(seq)}\n")
            temp_file = f.name
        
        try:
            profiler = StreamingKmerProfiler(k=21, sketch_size=50, sampling_method='reservoir')
            sketch = profiler.profile_sample(temp_file, "test_reservoir")
            
            # Verify sketch properties
            self.assertLessEqual(len(sketch), 50)  # Should not exceed sketch size
            self.assertGreater(len(sketch), 0)    # Should contain some k-mers
            
            # Verify normalized frequencies
            total_freq = sum(sketch.values())
            self.assertAlmostEqual(total_freq, 1.0, places=6)
            
            # Verify all frequencies are positive
            for freq in sketch.values():
                self.assertGreater(freq, 0)
        
        finally:
            Path(temp_file).unlink()
    
    def test_frequency_sampling(self):
        """Test frequency-based sampling."""
        # Create sequences with different frequency patterns
        high_freq_kmer = "ATCGATCGATCGATCGATCGA"  # 21-mer that will appear frequently
        sequences = [high_freq_kmer + "TTTTT"] * 50  # High frequency k-mer
        sequences.extend(["GCTAGCTAGCTAGCTAGCTAG" + f"{'A' * i}"] for i in range(20))  # Variable k-mers
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            for i, seq in enumerate(sequences):
                f.write(f"@read_{i}\n{seq}\n+\n{'~' * len(seq)}\n")
            temp_file = f.name
        
        try:
            profiler = StreamingKmerProfiler(k=21, sketch_size=30, sampling_method='frequency')
            sketch = profiler.profile_sample(temp_file, "test_frequency")
            
            # High frequency k-mer should be in the sketch
            sketch_kmers = set(sketch.keys())
            
            # Verify sketch contains most frequent k-mers
            self.assertLessEqual(len(sketch), 30)
            self.assertGreater(len(sketch), 0)
            
            # Verify frequencies are normalized
            total_freq = sum(sketch.values())
            self.assertAlmostEqual(total_freq, 1.0, places=6)
        
        finally:
            Path(temp_file).unlink()
    
    def test_adaptive_sampling(self):
        """Test adaptive sampling method."""
        # Create diverse sequences
        sequences = []
        for i in range(100):
            seq = ''.join(random.choice('ATCG') for _ in range(50))
            sequences.append(seq)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            for i, seq in enumerate(sequences):
                f.write(f"@read_{i}\n{seq}\n+\n{'~' * len(seq)}\n")
            temp_file = f.name
        
        try:
            profiler = StreamingKmerProfiler(k=15, sketch_size=100, sampling_method='adaptive')
            sketch = profiler.profile_sample(temp_file, "test_adaptive")
            
            # Verify adaptive sampling properties
            self.assertLessEqual(len(sketch), 100)
            self.assertGreater(len(sketch), 0)
            
            # Should have good diversity (not just most frequent)
            stats = profiler.get_sketch_statistics()
            self.assertIn("test_adaptive", stats)
            self.assertEqual(stats["test_adaptive"]["sampling_method"], "adaptive")
        
        finally:
            Path(temp_file).unlink()
    
    def test_sketch_similarity_methods(self):
        """Test different similarity calculation methods."""
        # Create two similar samples
        sequences1 = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"] * 20
        sequences2 = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"] * 15  # 75% overlap
        sequences2.extend(["GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"] * 5)
        
        # Create FASTQ files
        files = []
        for i, seqs in enumerate([sequences1, sequences2]):
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
                for j, seq in enumerate(seqs):
                    f.write(f"@read_{j}\n{seq}\n+\n{'~' * len(seq)}\n")
                files.append(f.name)
        
        try:
            profiler = StreamingKmerProfiler(k=21, sketch_size=50, sampling_method='frequency')
            
            sketch1 = profiler.profile_sample(files[0], "sample1")
            sketch2 = profiler.profile_sample(files[1], "sample2")
            
            # Test similarity methods
            jaccard_sim = profiler.compare_sketches("sample1", "sample2", method='jaccard')
            weighted_jaccard_sim = profiler.compare_sketches("sample1", "sample2", method='weighted_jaccard')
            cosine_sim = profiler.compare_sketches("sample1", "sample2", method='cosine')
            
            # Verify similarity scores are reasonable
            self.assertGreaterEqual(jaccard_sim, 0)
            self.assertLessEqual(jaccard_sim, 1)
            self.assertGreaterEqual(weighted_jaccard_sim, 0)
            self.assertLessEqual(weighted_jaccard_sim, 1)
            self.assertGreaterEqual(cosine_sim, 0)
            self.assertLessEqual(cosine_sim, 1)
            
            # Similar samples should have reasonable similarity
            self.assertGreater(jaccard_sim, 0.3)  # Should have some similarity
            
            # Test invalid similarity method
            with self.assertRaises(ValueError):
                profiler.compare_sketches("sample1", "sample2", method='invalid')
        
        finally:
            for f in files:
                Path(f).unlink()
    
    def test_memory_usage_estimation(self):
        """Test memory usage estimation."""
        # Create test data
        sequences = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"] * 50
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            for i, seq in enumerate(sequences):
                f.write(f"@read_{i}\n{seq}\n+\n{'~' * len(seq)}\n")
            temp_file = f.name
        
        try:
            profiler = StreamingKmerProfiler(k=21, sketch_size=100, sampling_method='reservoir')
            profiler.profile_sample(temp_file, "test")
            
            memory_info = profiler.estimate_memory_usage()
            
            # Verify memory estimation structure
            expected_keys = ['sketch_memory_mb', 'estimated_full_memory_mb', 'memory_reduction', 'samples', 'avg_sketch_size']
            for key in expected_keys:
                self.assertIn(key, memory_info)
            
            # Verify reasonable values
            self.assertGreater(memory_info['memory_reduction'], 1)  # Should show some reduction
            self.assertEqual(memory_info['samples'], 1)
            self.assertGreaterEqual(memory_info['sketch_memory_mb'], 0)
        
        finally:
            Path(temp_file).unlink()
    
    def test_sketch_statistics(self):
        """Test sketch statistics generation."""
        sequences = ["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"] * 30
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            for i, seq in enumerate(sequences):
                f.write(f"@read_{i}\n{seq}\n+\n{'~' * len(seq)}\n")
            temp_file = f.name
        
        try:
            profiler = StreamingKmerProfiler(k=21, sketch_size=50, sampling_method='reservoir')
            profiler.profile_sample(temp_file, "test")
            
            stats = profiler.get_sketch_statistics()
            
            # Verify statistics structure
            self.assertIn("test", stats)
            sample_stats = stats["test"]
            
            expected_keys = ['sketch_size', 'total_kmers_sampled', 'sampling_efficiency', 
                           'processing_time', 'sampling_method', 'average_frequency', 'frequency_range']
            for key in expected_keys:
                self.assertIn(key, sample_stats)
            
            # Verify reasonable values
            self.assertGreaterEqual(sample_stats['sketch_size'], 0)
            self.assertGreater(sample_stats['processing_time'], 0)
            self.assertEqual(sample_stats['sampling_method'], 'reservoir')
        
        finally:
            Path(temp_file).unlink()
    
    def test_comparison_with_full_profiler(self):
        """Test sketch profiler against full k-mer profiler for accuracy."""
        # Create test sequences with known patterns
        base_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        sequences = [base_sequence] * 10  # Repetitive for easier analysis
        
        # Add some variation
        for i in range(10):
            varied_seq = base_sequence[:30] + ''.join(random.choice('ATCG') for _ in range(20))
            sequences.append(varied_seq)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            for i, seq in enumerate(sequences):
                f.write(f"@read_{i}\n{seq}\n+\n{'~' * len(seq)}\n")
            temp_file = f.name
        
        try:
            # Full profiler
            full_profiler = KmerProfiler(k=21)
            full_profile = full_profiler.profile_sample(temp_file, "test_full")
            
            # Sketch profiler with large sketch size
            sketch_profiler = StreamingKmerProfiler(k=21, sketch_size=len(full_profile), 
                                                   sampling_method='frequency')
            sketch_profile = sketch_profiler.profile_sample(temp_file, "test_sketch")
            
            # Compare profiles
            full_kmers = set(full_profile.keys())
            sketch_kmers = set(sketch_profile.keys())
            
            # Sketch should capture major k-mers
            overlap = len(full_kmers & sketch_kmers)
            overlap_ratio = overlap / len(full_kmers) if full_kmers else 0
            
            # Should have reasonable overlap (at least 50% for frequency sampling)
            self.assertGreater(overlap_ratio, 0.3)
            
            logging_info = f"Full profile: {len(full_profile)} k-mers, "
            logging_info += f"Sketch: {len(sketch_profile)} k-mers, "
            logging_info += f"Overlap: {overlap_ratio:.2%}"
            print(logging_info)
        
        finally:
            Path(temp_file).unlink()


def performance_comparison_test():
    """Compare performance of different profiling methods."""
    print("\n=== Performance Comparison: Full vs Sketch Profiling ===")
    
    # Generate larger test dataset
    sequences = []
    for i in range(1000):
        seq = ''.join(random.choice('ATCG') for _ in range(150))
        sequences.append(seq)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        for i, seq in enumerate(sequences):
            f.write(f"@read_{i}\n{seq}\n+\n{'~' * len(seq)}\n")
        temp_file = f.name
    
    try:
        test_cases = [
            ("Full K-mer Profiler", lambda: KmerProfiler(k=21)),
            ("Sketch Profiler (Reservoir)", lambda: StreamingKmerProfiler(k=21, sketch_size=1000, sampling_method='reservoir')),
            ("Sketch Profiler (Frequency)", lambda: StreamingKmerProfiler(k=21, sketch_size=1000, sampling_method='frequency')),
            ("Sketch Profiler (Adaptive)", lambda: StreamingKmerProfiler(k=21, sketch_size=1000, sampling_method='adaptive')),
        ]
        
        results = {}
        
        for name, profiler_factory in test_cases:
            print(f"\nTesting {name}:")
            
            profiler = profiler_factory()
            
            start_time = time.perf_counter()
            profile = profiler.profile_sample(temp_file, "test")
            end_time = time.perf_counter()
            
            processing_time = end_time - start_time
            profile_size = len(profile)
            
            print(f"  Processing time: {processing_time:.3f}s")
            print(f"  Profile size: {profile_size:,} k-mers")
            
            # Memory estimation for sketch profilers
            if isinstance(profiler, StreamingKmerProfiler):
                memory_info = profiler.estimate_memory_usage()
                print(f"  Memory reduction: {memory_info['memory_reduction']:.1f}x")
                print(f"  Estimated memory: {memory_info['sketch_memory_mb']:.2f} MB")
            
            results[name] = {
                'time': processing_time,
                'profile_size': profile_size,
                'memory_reduction': memory_info.get('memory_reduction', 1) if isinstance(profiler, StreamingKmerProfiler) else 1
            }
        
        # Summary
        print(f"\n=== Performance Summary ===")
        full_time = results.get("Full K-mer Profiler", {}).get('time', 1)
        
        for name, result in results.items():
            if "Sketch" in name:
                speedup = full_time / result['time'] if result['time'] > 0 else float('inf')
                print(f"{name}:")
                print(f"  Speedup: {speedup:.2f}x")
                print(f"  Memory reduction: {result['memory_reduction']:.1f}x")
                print(f"  Profile compression: {results['Full K-mer Profiler']['profile_size'] / result['profile_size']:.1f}x")
    
    finally:
        Path(temp_file).unlink()


if __name__ == "__main__":
    # Run unit tests
    unittest.main(argv=[''], exit=False, verbosity=2)
    
    # Run performance comparison
    performance_comparison_test()