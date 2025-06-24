#!/usr/bin/env python3
"""
Streaming k-mer profiler with sketching for memory-efficient analysis.

This module provides massive memory reduction (10-100x) by maintaining
fixed-size sketches instead of full k-mer profiles.
"""

import logging
import random
import time
import heapq
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Optional, Union, Set
from pathlib import Path

from .profiler import KmerProfiler


class StreamingKmerProfiler(KmerProfiler):
    """Memory-efficient k-mer profiler using sketching techniques."""
    
    def __init__(self, k: int = 21, sketch_size: int = 10000, 
                 sampling_method: str = 'reservoir', **kwargs):
        """
        Initialize streaming k-mer profiler.
        
        Args:
            k: K-mer size
            sketch_size: Maximum number of k-mers to keep per sample
            sampling_method: Sampling strategy ('reservoir', 'frequency', 'adaptive')
            **kwargs: Additional arguments passed to KmerProfiler
        """
        super().__init__(k=k, **kwargs)
        self.sketch_size = sketch_size
        self.sampling_method = sampling_method
        # Use self.profiles from parent class to avoid duplication
        self.sketches = self.profiles  # Alias to avoid duplication
        self.sketch_metadata = {}  # Store metadata about sketches
        
        logging.info(f"StreamingKmerProfiler initialized: k={k}, sketch_size={sketch_size}, "
                    f"method={sampling_method}")
    
    def profile_sample(self, filepath: Union[str, List[str]], sample_name: str, 
                      memory_efficient: bool = True) -> Dict[str, float]:
        """
        Create k-mer sketch for a sample instead of full profile.
        
        Args:
            filepath: Path(s) to FASTQ file(s)
            sample_name: Name of the sample
            memory_efficient: Always True for sketching (kept for compatibility)
            
        Returns:
            Dictionary representing k-mer sketch with frequencies
        """
        if isinstance(filepath, list):
            logging.info(f"Processing paired-end sample: {sample_name} ({len(filepath)} files)")
        else:
            logging.info(f"Processing single-end sample: {sample_name}")
        
        start_time = time.time()
        
        # Create sketch based on sampling method
        if self.sampling_method == 'reservoir':
            sketch = self._reservoir_sampling(filepath)
        elif self.sampling_method == 'frequency':
            sketch = self._frequency_sampling(filepath)
        elif self.sampling_method == 'adaptive':
            sketch = self._adaptive_sampling(filepath)
        else:
            raise ValueError(f"Unknown sampling method: {self.sampling_method}")
        
        # Normalize sketch to frequencies
        total_count = sum(sketch.values())
        if total_count > 0:
            normalized_sketch = {kmer: count / total_count for kmer, count in sketch.items()}
        else:
            normalized_sketch = {}
        
        # Store sketch and metadata
        self.profiles[sample_name] = normalized_sketch  # Store in profiles (sketches is an alias)
        self.sketch_metadata[sample_name] = {
            'total_kmers_sampled': total_count,
            'unique_kmers': len(sketch),
            'processing_time': time.time() - start_time,
            'sampling_method': self.sampling_method,
            'sketch_size_used': len(sketch)
        }
        
        if sample_name not in self.sample_names:
            self.sample_names.append(sample_name)
        
        logging.info(f"Created sketch with {len(normalized_sketch)} k-mers "
                    f"(from {total_count} total) in {time.time() - start_time:.2f}s")
        
        return normalized_sketch
    
    def _reservoir_sampling(self, filepath: Union[str, List[str]]) -> Dict[str, int]:
        """
        Unbiased reservoir sampling of k-mers.
        
        This maintains a representative sample of k-mers while processing
        sequences in a streaming fashion.
        """
        reservoir = []
        kmer_count = 0
        
        # Process all files
        file_paths = [filepath] if isinstance(filepath, str) else filepath
        
        for file_path in file_paths:
            for sequence in self._stream_fastq_records(file_path):
                seq_kmers = self._extract_kmers(sequence)
                
                # Add each k-mer occurrence to reservoir sampling
                for kmer, freq in seq_kmers.items():
                    for _ in range(freq):
                        kmer_count += 1
                        
                        if len(reservoir) < self.sketch_size:
                            reservoir.append(kmer)
                        else:
                            # Reservoir sampling: replace with probability sketch_size/kmer_count
                            j = random.randint(0, kmer_count - 1)
                            if j < self.sketch_size:
                                reservoir[j] = kmer
        
        # Convert to frequency dictionary
        return Counter(reservoir)
    
    def _frequency_sampling(self, filepath: Union[str, List[str]]) -> Dict[str, int]:
        """
        Sample the most frequent k-mers using a single-pass algorithm with min-heap.
        
        This optimized version uses a min-heap to maintain the top-k k-mers
        in a single pass, reducing memory usage and I/O overhead.
        """
        # Use a min-heap to keep top k-mers by frequency
        # Heap items are (count, kmer) tuples
        top_kmers_heap = []
        kmer_counts = {}  # Track counts for k-mers in the heap
        total_sequences = 0
        unique_kmers_seen = 0
        
        file_paths = [filepath] if isinstance(filepath, str) else filepath
        
        for file_path in file_paths:
            for sequence in self._stream_fastq_records(file_path):
                seq_kmers = self._extract_kmers(sequence)
                total_sequences += 1
                
                for kmer, count in seq_kmers.items():
                    if kmer in kmer_counts:
                        # K-mer is already in our top set, update its count
                        old_count = kmer_counts[kmer]
                        kmer_counts[kmer] = old_count + count
                        # Note: heap will have stale entries, but we handle this during extraction
                    else:
                        unique_kmers_seen += 1
                        # New k-mer
                        if len(top_kmers_heap) < self.sketch_size:
                            # Heap not full yet, add the k-mer
                            heapq.heappush(top_kmers_heap, (count, kmer))
                            kmer_counts[kmer] = count
                        else:
                            # Heap is full, check if this k-mer should replace the minimum
                            min_count, min_kmer = top_kmers_heap[0]
                            if count > min_count:
                                # Remove the minimum and add the new k-mer
                                heapq.heapreplace(top_kmers_heap, (count, kmer))
                                del kmer_counts[min_kmer]
                                kmer_counts[kmer] = count
        
        logging.debug(f"Frequency sampling: processed {total_sequences} sequences, "
                     f"saw {unique_kmers_seen} unique k-mers, kept top {len(kmer_counts)}")
        
        # Extract final k-mer counts
        # Note: The heap may have stale entries with old counts, so we use kmer_counts
        return dict(kmer_counts)
    
    def _adaptive_sampling(self, filepath: Union[str, List[str]]) -> Dict[str, int]:
        """
        Adaptive sampling combining frequency and diversity.
        
        Uses 70% most frequent k-mers and 30% diverse k-mers to balance
        signal preservation with diversity.
        """
        frequency_fraction = 0.7
        diversity_fraction = 0.3
        
        frequency_size = int(self.sketch_size * frequency_fraction)
        diversity_size = self.sketch_size - frequency_size
        
        # Get frequency-based sample
        freq_sketch = self._frequency_sampling_limited(filepath, frequency_size)
        
        # Get diversity sample (excluding already selected k-mers)
        diversity_sketch = self._diversity_sampling(filepath, diversity_size, 
                                                   exclude=set(freq_sketch.keys()))
        
        # Combine sketches
        combined_sketch = freq_sketch.copy()
        combined_sketch.update(diversity_sketch)
        
        return combined_sketch
    
    def _frequency_sampling_limited(self, filepath: Union[str, List[str]], 
                                  limit: int) -> Dict[str, int]:
        """Frequency sampling with a specific limit."""
        # Temporarily adjust sketch size
        original_size = self.sketch_size
        self.sketch_size = limit
        
        result = self._frequency_sampling(filepath)
        
        # Restore original size
        self.sketch_size = original_size
        
        return result
    
    def _diversity_sampling(self, filepath: Union[str, List[str]], 
                           limit: int, exclude: Set[str]) -> Dict[str, int]:
        """
        Sample k-mers for diversity (avoiding already selected ones).
        
        Uses a modified reservoir sampling that excludes already selected k-mers.
        """
        reservoir = []
        kmer_count = 0
        
        file_paths = [filepath] if isinstance(filepath, str) else filepath
        
        for file_path in file_paths:
            for sequence in self._stream_fastq_records(file_path):
                seq_kmers = self._extract_kmers(sequence)
                
                for kmer, freq in seq_kmers.items():
                    if kmer in exclude:
                        continue  # Skip already selected k-mers
                    
                    for _ in range(freq):
                        kmer_count += 1
                        
                        if len(reservoir) < limit:
                            reservoir.append(kmer)
                        else:
                            j = random.randint(0, kmer_count - 1)
                            if j < limit:
                                reservoir[j] = kmer
        
        return Counter(reservoir)
    
    def get_sketch_statistics(self) -> Dict[str, Dict]:
        """Get statistics about all sketches."""
        stats = {}
        
        for sample_name, metadata in self.sketch_metadata.items():
            sketch = self.sketches.get(sample_name, {})
            
            stats[sample_name] = {
                'sketch_size': len(sketch),
                'total_kmers_sampled': metadata['total_kmers_sampled'],
                'sampling_efficiency': len(sketch) / self.sketch_size if self.sketch_size > 0 else 0,
                'processing_time': metadata['processing_time'],
                'sampling_method': metadata['sampling_method'],
                'average_frequency': sum(sketch.values()) / len(sketch) if sketch else 0,
                'frequency_range': (min(sketch.values()), max(sketch.values())) if sketch else (0, 0)
            }
        
        return stats
    
    def estimate_memory_usage(self) -> Dict[str, float]:
        """Estimate memory usage of sketches vs full profiles."""
        if not self.sketches:
            return {'sketch_memory_mb': 0, 'estimated_full_memory_mb': 0, 'memory_reduction': 0}
        
        # Estimate sketch memory usage
        sketch_memory = 0
        for sketch in self.sketches.values():
            # Estimate: 21 bytes per k-mer (string) + 8 bytes (float) = 29 bytes
            sketch_memory += len(sketch) * 29
        
        # Estimate full profile memory (assume 100k k-mers per sample average)
        estimated_kmers_per_sample = 100000
        full_memory = len(self.sketches) * estimated_kmers_per_sample * 29
        
        sketch_memory_mb = sketch_memory / (1024 * 1024)
        full_memory_mb = full_memory / (1024 * 1024)
        
        memory_reduction = full_memory_mb / sketch_memory_mb if sketch_memory_mb > 0 else 0
        
        return {
            'sketch_memory_mb': sketch_memory_mb,
            'estimated_full_memory_mb': full_memory_mb,
            'memory_reduction': memory_reduction,
            'samples': len(self.sketches),
            'avg_sketch_size': sum(len(s) for s in self.sketches.values()) / len(self.sketches)
        }
    
    def compare_sketches(self, sample1: str, sample2: str, 
                        method: str = 'jaccard') -> float:
        """
        Compare two sketches using specified similarity metric.
        
        Args:
            sample1: Name of first sample
            sample2: Name of second sample  
            method: Similarity method ('jaccard', 'weighted_jaccard', 'cosine')
            
        Returns:
            Similarity score between 0 and 1
        """
        if sample1 not in self.sketches or sample2 not in self.sketches:
            raise ValueError(f"Sample not found in sketches")
        
        sketch1 = self.sketches[sample1]
        sketch2 = self.sketches[sample2]
        
        if method == 'jaccard':
            return self._jaccard_similarity(sketch1, sketch2)
        elif method == 'weighted_jaccard':
            return self._weighted_jaccard_similarity(sketch1, sketch2)
        elif method == 'cosine':
            return self._cosine_similarity(sketch1, sketch2)
        else:
            raise ValueError(f"Unknown similarity method: {method}")
    
    def _jaccard_similarity(self, sketch1: Dict[str, float], 
                           sketch2: Dict[str, float]) -> float:
        """Jaccard similarity between sketches."""
        set1, set2 = set(sketch1.keys()), set(sketch2.keys())
        intersection = len(set1 & set2)
        union = len(set1 | set2)
        return intersection / union if union > 0 else 0
    
    def _weighted_jaccard_similarity(self, sketch1: Dict[str, float], 
                                   sketch2: Dict[str, float]) -> float:
        """Weighted Jaccard similarity considering frequencies."""
        all_kmers = set(sketch1.keys()) | set(sketch2.keys())
        
        min_sum = sum(min(sketch1.get(k, 0), sketch2.get(k, 0)) for k in all_kmers)
        max_sum = sum(max(sketch1.get(k, 0), sketch2.get(k, 0)) for k in all_kmers)
        
        return min_sum / max_sum if max_sum > 0 else 0
    
    def _cosine_similarity(self, sketch1: Dict[str, float], 
                          sketch2: Dict[str, float]) -> float:
        """Cosine similarity between frequency vectors."""
        import math
        
        all_kmers = set(sketch1.keys()) | set(sketch2.keys())
        
        dot_product = sum(sketch1.get(k, 0) * sketch2.get(k, 0) for k in all_kmers)
        norm1 = math.sqrt(sum(sketch1.get(k, 0)**2 for k in all_kmers))
        norm2 = math.sqrt(sum(sketch2.get(k, 0)**2 for k in all_kmers))
        
        return dot_product / (norm1 * norm2) if norm1 > 0 and norm2 > 0 else 0


def process_sample_worker_sketch(args):
    """
    Worker function for multiprocessing sketch-based sample processing.
    
    Args:
        args: Tuple containing processing parameters for StreamingKmerProfiler
    
    Returns:
        Tuple containing (sample_name, sketch, error_message)
    """
    (filepath, sample_name, k, max_reads, min_kmer_freq, sketch_size, 
     sampling_method, use_fast_extraction, progress_dict, lock) = args
    
    try:
        # Create a temporary sketch profiler for this worker
        profiler = StreamingKmerProfiler(
            k=k, 
            max_reads=max_reads, 
            min_kmer_freq=min_kmer_freq,
            sketch_size=sketch_size,
            sampling_method=sampling_method,
            use_fast_extraction=use_fast_extraction
        )
        sketch = profiler.profile_sample(filepath, sample_name)
        
        # Update progress safely
        with lock:
            progress_dict['completed'] += 1
            progress_dict['current_sample'] = sample_name
        
        return (sample_name, sketch, None)
    
    except Exception as e:
        # Update progress and record error
        with lock:
            progress_dict['completed'] += 1
            progress_dict['failed'] += 1
            progress_dict['errors'][sample_name] = str(e)
        
        return (sample_name, None, str(e))


if __name__ == "__main__":
    # Quick test
    logging.basicConfig(level=logging.INFO)
    
    # Test with example data
    import tempfile
    from pathlib import Path
    
    # Create test FASTQ
    test_sequences = [
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    ]
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        for i, seq in enumerate(test_sequences * 100):  # 500 total sequences
            f.write(f"@read_{i}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write("~" * len(seq) + "\n")
        
        temp_file = f.name
    
    try:
        # Test different sampling methods
        for method in ['reservoir', 'frequency', 'adaptive']:
            print(f"\nTesting {method} sampling:")
            
            profiler = StreamingKmerProfiler(
                k=21, 
                sketch_size=1000, 
                sampling_method=method
            )
            
            sketch = profiler.profile_sample(temp_file, f"test_{method}")
            stats = profiler.get_sketch_statistics()
            memory_info = profiler.estimate_memory_usage()
            
            print(f"  Sketch size: {len(sketch)}")
            print(f"  Memory reduction: {memory_info['memory_reduction']:.1f}x")
            print(f"  Processing time: {stats[f'test_{method}']['processing_time']:.3f}s")
    
    finally:
        Path(temp_file).unlink()  # Clean up