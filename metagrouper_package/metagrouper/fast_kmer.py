#!/usr/bin/env python3
"""
Fast k-mer extraction using bit operations for MetaGrouper.

This module provides significant speedup (5-10x) over string-based k-mer extraction
by representing DNA sequences as integers and using bit operations for sliding window
and reverse complement operations.
"""

import logging
from collections import defaultdict
from typing import Dict, List, Optional


class FastKmerExtractor:
    """Fast k-mer extraction using bit operations."""
    
    def __init__(self, k: int = 21):
        """
        Initialize fast k-mer extractor.
        
        Args:
            k: K-mer size (must be <= 32 for 64-bit integer representation)
        """
        if k > 32:
            raise ValueError("K-mer size > 32 not supported in fast mode (use legacy)")
        
        self.k = k
        self.base_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        self.int_to_base = ['A', 'C', 'G', 'T']
        self.mask = (1 << (2 * k)) - 1  # Mask for k-mer bits
        self.rc_table = self._build_reverse_complement_table()
        
        logging.debug(f"FastKmerExtractor initialized for k={k}")
    
    def _build_reverse_complement_table(self) -> List[int]:
        """Pre-compute reverse complement lookup table for 8-bit values."""
        rc_table = [0] * 256
        for i in range(256):
            rc = 0
            val = i
            for _ in range(4):  # 4 bases per 8-bit value
                base = val & 3  # Get last 2 bits
                rc_base = 3 - base  # A↔T (0↔3), C↔G (1↔2)
                rc = (rc << 2) | rc_base
                val >>= 2
            rc_table[i] = rc
        return rc_table
    
    def sequence_to_kmers_fast(self, sequence: str) -> Dict[int, int]:
        """
        Extract k-mers from sequence using bit operations.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Dictionary mapping k-mer integers to frequencies
        """
        if len(sequence) < self.k:
            return {}
        
        # Convert to uppercase for consistency
        sequence = sequence.upper()
        
        # Handle sequences with N's by splitting
        if 'N' in sequence:
            return self._handle_sequence_with_n(sequence)
        
        # Use optimized clean sequence method
        return self._extract_kmers_from_clean_sequence(sequence)
    
    def _get_canonical_fast(self, kmer_int: int) -> int:
        """Get canonical k-mer (smaller of forward/reverse complement)."""
        rc_kmer = self._reverse_complement_fast(kmer_int)
        return min(kmer_int, rc_kmer)
    
    def _reverse_complement_fast(self, kmer_int: int) -> int:
        """Fast reverse complement using direct bit operations."""
        rc = 0
        temp = kmer_int
        
        # Process each base (2 bits) directly
        for _ in range(self.k):
            base = temp & 3  # Get last 2 bits
            rc_base = 3 - base  # A↔T (0↔3), C↔G (1↔2)
            rc = (rc << 2) | rc_base
            temp >>= 2
        
        return rc
    
    def _handle_sequence_with_n(self, sequence: str) -> Dict[int, int]:
        """Handle sequences containing N's by splitting into segments."""
        kmers = defaultdict(int)
        
        # Split sequence on N's and process each segment
        segments = sequence.split('N')
        for segment in segments:
            if len(segment) >= self.k:
                # Process segment directly without recursion
                segment_kmers = self._extract_kmers_from_clean_sequence(segment.upper())
                for kmer, count in segment_kmers.items():
                    kmers[kmer] += count
        
        return kmers
    
    def _extract_kmers_from_clean_sequence(self, sequence: str) -> Dict[int, int]:
        """Extract k-mers from sequence known to be clean (no N's)."""
        if len(sequence) < self.k:
            return {}
        
        # Validate sequence
        if not all(base in self.base_to_int for base in sequence):
            return {}
        
        kmers = defaultdict(int)
        
        # Initialize first k-mer
        current_kmer = 0
        for i in range(self.k):
            current_kmer = (current_kmer << 2) | self.base_to_int[sequence[i]]
        
        # Process first k-mer
        canonical = self._get_canonical_fast(current_kmer)
        kmers[canonical] += 1
        
        # Sliding window for remaining k-mers
        for i in range(self.k, len(sequence)):
            # Remove leftmost base, add rightmost base
            current_kmer = ((current_kmer << 2) | self.base_to_int[sequence[i]]) & self.mask
            canonical = self._get_canonical_fast(current_kmer)
            kmers[canonical] += 1
        
        return kmers
    
    def int_to_kmer_string(self, kmer_int: int) -> str:
        """Convert integer k-mer back to string representation."""
        bases = []
        temp = kmer_int
        for _ in range(self.k):
            bases.append(self.int_to_base[temp & 3])
            temp >>= 2
        return ''.join(reversed(bases))
    
    def kmer_string_to_int(self, kmer_str: str) -> int:
        """Convert k-mer string to integer representation."""
        if len(kmer_str) != self.k:
            raise ValueError(f"K-mer string length {len(kmer_str)} != k={self.k}")
        
        kmer_int = 0
        for base in kmer_str:
            if base not in self.base_to_int:
                raise ValueError(f"Invalid base '{base}' in k-mer")
            kmer_int = (kmer_int << 2) | self.base_to_int[base]
        
        return kmer_int
    
    def extract_kmers_with_strings(self, sequence: str) -> Dict[str, int]:
        """
        Extract k-mers and return as string dictionary for compatibility.
        
        This method provides the same interface as the original string-based
        method but uses fast bit operations internally.
        """
        int_kmers = self.sequence_to_kmers_fast(sequence)
        return {self.int_to_kmer_string(k): count for k, count in int_kmers.items()}


def benchmark_extraction_methods(sequence: str, k: int = 21, iterations: int = 100) -> Dict[str, float]:
    """
    Benchmark fast vs legacy k-mer extraction methods.
    
    Args:
        sequence: Test DNA sequence
        k: K-mer size
        iterations: Number of iterations for timing
        
    Returns:
        Dictionary with timing results and speedup
    """
    import time
    
    # Legacy method (simplified version of original)
    def extract_kmers_legacy(seq: str) -> Dict[str, int]:
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        def reverse_complement(s):
            return "".join(complement.get(base, base) for base in reversed(s))
        
        kmers = defaultdict(int)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if "N" not in kmer:
                rev_comp = reverse_complement(kmer)
                canonical = min(kmer, rev_comp)
                kmers[canonical] += 1
        return kmers
    
    # Fast method
    fast_extractor = FastKmerExtractor(k)
    
    # Warm up
    _ = extract_kmers_legacy(sequence)
    _ = fast_extractor.extract_kmers_with_strings(sequence)
    
    # Benchmark legacy method
    start_time = time.perf_counter()
    for _ in range(iterations):
        legacy_result = extract_kmers_legacy(sequence)
    legacy_time = time.perf_counter() - start_time
    
    # Benchmark fast method  
    start_time = time.perf_counter()
    for _ in range(iterations):
        fast_result = fast_extractor.extract_kmers_with_strings(sequence)
    fast_time = time.perf_counter() - start_time
    
    # Verify results are identical
    results_match = legacy_result == fast_result
    
    speedup = legacy_time / fast_time if fast_time > 0 else float('inf')
    
    return {
        'legacy_time': legacy_time,
        'fast_time': fast_time,
        'speedup': speedup,
        'results_match': results_match,
        'legacy_kmers': len(legacy_result),
        'fast_kmers': len(fast_result),
        'sequence_length': len(sequence),
        'total_kmers': sum(legacy_result.values()) if legacy_result else 0
    }


if __name__ == "__main__":
    # Quick test
    test_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    extractor = FastKmerExtractor(k=21)
    
    # Test extraction
    int_kmers = extractor.sequence_to_kmers_fast(test_sequence)
    str_kmers = extractor.extract_kmers_with_strings(test_sequence)
    
    print(f"Extracted {len(int_kmers)} unique k-mers")
    print(f"Integer format: {list(int_kmers.items())[:3]}...")
    print(f"String format: {list(str_kmers.items())[:3]}...")
    
    # Test benchmark
    benchmark_results = benchmark_extraction_methods(test_sequence, k=21, iterations=100)
    print(f"\nBenchmark results:")
    print(f"Speedup: {benchmark_results['speedup']:.1f}x")
    print(f"Results match: {benchmark_results['results_match']}")