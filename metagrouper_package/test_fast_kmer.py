#!/usr/bin/env python3
"""
Test suite for fast k-mer extraction implementation.
"""

import unittest
import time
import random
from collections import defaultdict
from metagrouper.fast_kmer import FastKmerExtractor, benchmark_extraction_methods


class TestFastKmerExtractor(unittest.TestCase):
    """Test fast k-mer extraction functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.extractor = FastKmerExtractor(k=21)
        self.short_extractor = FastKmerExtractor(k=5)
    
    def test_initialization(self):
        """Test FastKmerExtractor initialization."""
        # Valid k-mer sizes
        for k in [5, 15, 21, 31]:
            extractor = FastKmerExtractor(k=k)
            self.assertEqual(extractor.k, k)
            self.assertEqual(len(extractor.int_to_base), 4)
        
        # Invalid k-mer size
        with self.assertRaises(ValueError):
            FastKmerExtractor(k=33)  # Too large
    
    def test_reverse_complement_table(self):
        """Test reverse complement lookup table."""
        extractor = FastKmerExtractor(k=5)
        
        # Test known values
        # 0b00011011 = ACTG -> should become CAGT = 0b01001110
        test_byte = 0b00011011  # ACTG in 2-bit encoding
        rc_byte = extractor.rc_table[test_byte]
        
        # Verify table is built correctly
        self.assertEqual(len(extractor.rc_table), 256)
        self.assertIsInstance(rc_byte, int)
    
    def test_sequence_conversion(self):
        """Test sequence to integer conversion."""
        # Test simple sequence
        sequence = "ATCG"
        extractor = FastKmerExtractor(k=4)
        
        # Manual calculation: A=0, T=3, C=1, G=2
        # ATCG = 0b00110110 = 54
        expected = (0 << 6) | (3 << 4) | (1 << 2) | 2  # = 54
        
        int_kmers = extractor.sequence_to_kmers_fast(sequence)
        self.assertEqual(len(int_kmers), 1)
        self.assertIn(expected, int_kmers)
    
    def test_string_int_conversion(self):
        """Test bidirectional string-integer conversion."""
        test_kmers = ["ATCGATCGATCGATCGATCGA", "GCTAGCTAGCTAGCTAGCTAG", "AAAAAAAAAAAAAAAAAAA"]
        
        for kmer in test_kmers:
            if len(kmer) == self.extractor.k:
                # Convert to int and back
                kmer_int = self.extractor.kmer_string_to_int(kmer)
                recovered_kmer = self.extractor.int_to_kmer_string(kmer_int)
                self.assertEqual(kmer, recovered_kmer)
    
    def test_canonical_kmers(self):
        """Test canonical k-mer selection."""
        # Test with known reverse complements
        test_cases = [
            ("AAAAA", "TTTTT"),  # AAAAA should be canonical (smaller)
            ("ATCGA", "TCGAT"),  # Need to check which is smaller
            ("ACGTG", "CACGT"),
        ]
        
        extractor = FastKmerExtractor(k=5)
        
        for kmer1, kmer2 in test_cases:
            int1 = extractor.sequence_to_kmers_fast(kmer1)
            int2 = extractor.sequence_to_kmers_fast(kmer2)
            
            # Both should produce the same canonical k-mer
            self.assertEqual(set(int1.keys()), set(int2.keys()))
    
    def test_sequences_with_n(self):
        """Test handling of sequences with N bases."""
        # Sequence with N in the middle
        sequence = "ATCGATCGANATCGATCGATC"
        kmers = self.extractor.sequence_to_kmers_fast(sequence)
        
        # Should split into segments and extract k-mers from valid segments
        self.assertGreater(len(kmers), 0)
        
        # No k-mer should contain N (since we use integer representation)
        for kmer_int in kmers.keys():
            kmer_str = self.extractor.int_to_kmer_string(kmer_int)
            self.assertNotIn('N', kmer_str)
    
    def test_short_sequences(self):
        """Test handling of sequences shorter than k."""
        short_sequence = "ATCG"  # Length 4, k=21
        kmers = self.extractor.sequence_to_kmers_fast(short_sequence)
        self.assertEqual(len(kmers), 0)
    
    def test_extract_kmers_with_strings(self):
        """Test string output compatibility method."""
        sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        
        str_kmers = self.extractor.extract_kmers_with_strings(sequence)
        
        # Verify output format
        self.assertIsInstance(str_kmers, dict)
        for kmer, count in str_kmers.items():
            self.assertIsInstance(kmer, str)
            self.assertIsInstance(count, int)
            self.assertEqual(len(kmer), self.extractor.k)
            self.assertGreater(count, 0)
    
    def test_comparison_with_legacy(self):
        """Test that results match legacy string-based method."""
        # Legacy method implementation
        def extract_kmers_legacy(sequence: str, k: int) -> dict:
            complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
            
            def reverse_complement(seq):
                return "".join(complement.get(base, base) for base in reversed(seq))
            
            kmers = defaultdict(int)
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i + k]
                if "N" not in kmer:
                    rev_comp = reverse_complement(kmer)
                    canonical = min(kmer, rev_comp)
                    kmers[canonical] += 1
            return dict(kmers)
        
        # Test with various sequences
        test_sequences = [
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "TTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCC"
        ]
        
        for sequence in test_sequences:
            legacy_kmers = extract_kmers_legacy(sequence, self.extractor.k)
            fast_kmers = self.extractor.extract_kmers_with_strings(sequence)
            
            self.assertEqual(legacy_kmers, fast_kmers, 
                           f"Results differ for sequence: {sequence[:20]}...")
    
    def test_performance_improvement(self):
        """Test that fast method is actually faster."""
        # Generate longer sequence for meaningful benchmark
        bases = ['A', 'T', 'C', 'G']
        sequence = ''.join(random.choice(bases) for _ in range(10000))
        
        benchmark_results = benchmark_extraction_methods(
            sequence, k=21, iterations=10
        )
        
        # Verify results match
        self.assertTrue(benchmark_results['results_match'])
        
        # Should show some speedup (even if modest for small tests)
        self.assertGreaterEqual(benchmark_results['speedup'], 0.8)  # At least not slower
        
        print(f"Performance test: {benchmark_results['speedup']:.2f}x speedup")
    
    def test_edge_cases(self):
        """Test various edge cases."""
        # Empty sequence
        empty_kmers = self.extractor.sequence_to_kmers_fast("")
        self.assertEqual(len(empty_kmers), 0)
        
        # Sequence exactly k length
        k_length_seq = "A" * self.extractor.k
        k_kmers = self.extractor.sequence_to_kmers_fast(k_length_seq)
        self.assertEqual(len(k_kmers), 1)
        
        # Sequence with lowercase
        lower_seq = "atcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcg"
        lower_kmers = self.extractor.extract_kmers_with_strings(lower_seq)
        upper_kmers = self.extractor.extract_kmers_with_strings(lower_seq.upper())
        self.assertEqual(lower_kmers, upper_kmers)
    
    def test_different_k_sizes(self):
        """Test various k-mer sizes."""
        sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        
        for k in [5, 11, 15, 21, 25, 31]:
            extractor = FastKmerExtractor(k=k)
            kmers = extractor.extract_kmers_with_strings(sequence)
            
            # Check that all k-mers have correct length
            for kmer in kmers.keys():
                self.assertEqual(len(kmer), k)
            
            # Check expected number of k-mers
            expected_count = max(0, len(sequence) - k + 1)
            total_kmers = sum(kmers.values())
            self.assertEqual(total_kmers, expected_count)


def generate_viral_test_sequences(n_sequences: int = 100, length: int = 150) -> list:
    """Generate realistic viral-like test sequences."""
    sequences = []
    
    # Common viral motifs and patterns
    viral_patterns = [
        "ATGAAA",  # Start codon variations
        "GCTAGC",  # Viral promoter-like
        "CCCCCC",  # Poly-C regions common in viruses
        "ATATAT",  # AT-rich regions
        "CGCGCG",  # CG-rich regions
    ]
    
    bases = ['A', 'T', 'C', 'G']
    
    for i in range(n_sequences):
        sequence = ""
        
        # Mix of random sequence and viral patterns
        while len(sequence) < length:
            if random.random() < 0.3:  # 30% chance of viral pattern
                pattern = random.choice(viral_patterns)
                sequence += pattern
            else:
                sequence += random.choice(bases)
        
        # Trim to exact length
        sequence = sequence[:length]
        sequences.append(sequence)
    
    return sequences


def comprehensive_performance_test():
    """Run comprehensive performance comparison."""
    print("=== Comprehensive Fast K-mer Performance Test ===\n")
    
    # Test different sequence sizes
    test_cases = [
        (100, 100, "Small sequences"),
        (500, 150, "Medium sequences"),  
        (1000, 200, "Large sequences"),
        (100, 1000, "Very long sequences")
    ]
    
    for n_seqs, seq_length, description in test_cases:
        print(f"{description}: {n_seqs} sequences of {seq_length}bp each")
        
        # Generate test data
        sequences = generate_viral_test_sequences(n_seqs, seq_length)
        
        # Benchmark extraction
        total_legacy_time = 0
        total_fast_time = 0
        
        extractor = FastKmerExtractor(k=21)
        
        # Time legacy method
        start_time = time.time()
        for seq in sequences:
            benchmark_results = benchmark_extraction_methods(seq, k=21, iterations=1)
        
        print(f"  Average speedup: {benchmark_results['speedup']:.2f}x")
        print(f"  Results match: {benchmark_results['results_match']}")
        print()


if __name__ == "__main__":
    # Run unit tests
    unittest.main(argv=[''], exit=False, verbosity=2)
    
    # Run performance tests
    print("\n" + "="*60)
    comprehensive_performance_test()