#!/usr/bin/env python3
"""
Sparse similarity analyzer for memory-efficient distance computation.

This module enables analysis of 340+ samples by storing only significant
similarities in sparse matrices, dramatically reducing memory requirements.
"""

import logging
import numpy as np
import time
from typing import Dict, List, Tuple, Optional, Any, Set
from scipy.sparse import csr_matrix, lil_matrix
from collections import defaultdict

try:
    from datasketch import MinHashLSH, MinHash
    MINHASH_AVAILABLE = True
except ImportError:
    MINHASH_AVAILABLE = False
    logging.warning("datasketch not available - LSH approximation disabled")


class SparseSimilarityAnalyzer:
    """
    Memory-efficient similarity analysis using sparse matrices.
    
    This class enables analysis of large datasets (340+ samples) by storing
    only similarities above a threshold in sparse format.
    """
    
    def __init__(self, similarity_threshold: float = 0.1, 
                 max_comparisons: Optional[int] = None,
                 use_approximation: bool = False):
        """
        Initialize sparse similarity analyzer.
        
        Args:
            similarity_threshold: Minimum similarity to store (0.0 to 1.0)
            max_comparisons: Maximum pairwise comparisons before using approximation
            use_approximation: Force use of approximation methods (LSH)
        """
        self.similarity_threshold = similarity_threshold
        self.max_comparisons = max_comparisons or 1000000  # 1M comparisons default
        self.use_approximation = use_approximation
        
        # Storage for analysis results
        self.similarity_matrix = None
        self.sample_names = []
        self.computation_stats = {}
        
        logging.info(f"SparseSimilarityAnalyzer initialized: threshold={similarity_threshold}, "
                    f"max_comparisons={max_comparisons}")
    
    def compute_similarities(self, profiles: Dict[str, Dict[str, float]], 
                           method: str = 'jaccard') -> Tuple[csr_matrix, List[str]]:
        """
        Compute sparse similarity matrix from k-mer profiles/sketches.
        
        Args:
            profiles: Dictionary of sample_name -> k-mer_profile mappings
            method: Similarity method ('jaccard', 'weighted_jaccard', 'cosine')
            
        Returns:
            Tuple of (sparse_similarity_matrix, sample_names)
        """
        if not profiles:
            raise ValueError("No profiles provided")
        
        self.sample_names = list(profiles.keys())
        n_samples = len(self.sample_names)
        total_comparisons = n_samples * (n_samples - 1) // 2
        
        logging.info(f"Computing similarities for {n_samples} samples "
                    f"({total_comparisons:,} comparisons)")
        
        # Choose computation method based on dataset size
        if (self.use_approximation or 
            total_comparisons > self.max_comparisons or 
            (n_samples > 500 and MINHASH_AVAILABLE)):
            
            logging.info("Using LSH approximation for large dataset")
            similarity_matrix = self._compute_lsh_similarities(profiles, method)
        else:
            logging.info("Using exact sparse computation")
            similarity_matrix = self._compute_exact_sparse_similarities(profiles, method)
        
        self.similarity_matrix = similarity_matrix
        return similarity_matrix, self.sample_names
    
    def _compute_exact_sparse_similarities(self, profiles: Dict[str, Dict[str, float]], 
                                         method: str) -> csr_matrix:
        """Compute exact pairwise similarities, store only significant ones."""
        n_samples = len(self.sample_names)
        
        # Use LIL format for efficient construction
        similarity_matrix = lil_matrix((n_samples, n_samples))
        
        start_time = time.time()
        completed_comparisons = 0
        significant_pairs = 0
        
        # Progress tracking
        total_comparisons = n_samples * (n_samples - 1) // 2
        next_report = total_comparisons // 20  # Report every 5%
        
        for i in range(n_samples):
            profile_i = profiles[self.sample_names[i]]
            
            for j in range(i + 1, n_samples):
                profile_j = profiles[self.sample_names[j]]
                
                # Compute similarity
                if method == 'jaccard':
                    sim = self._jaccard_similarity(profile_i, profile_j)
                elif method == 'weighted_jaccard':
                    sim = self._weighted_jaccard_similarity(profile_i, profile_j)
                elif method == 'cosine':
                    sim = self._cosine_similarity(profile_i, profile_j)
                elif method == 'braycurtis':
                    # Convert Bray-Curtis distance to similarity
                    distance = self._braycurtis_distance(profile_i, profile_j)
                    sim = 1 - distance
                elif method == 'euclidean':
                    # Convert Euclidean distance to similarity
                    distance = self._euclidean_distance(profile_i, profile_j)
                    # Normalize to [0,1] similarity (assuming max distance of 2)
                    sim = max(0, 1 - distance / 2)
                else:
                    raise ValueError(f"Unknown similarity method: {method}. "
                                   f"Supported: jaccard, weighted_jaccard, cosine, braycurtis, euclidean")
                
                # Store if above threshold
                if sim > self.similarity_threshold:
                    similarity_matrix[i, j] = sim
                    similarity_matrix[j, i] = sim  # Symmetric
                    significant_pairs += 1
                
                completed_comparisons += 1
                
                # Progress reporting
                if completed_comparisons % next_report == 0:
                    progress = completed_comparisons / total_comparisons * 100
                    elapsed = time.time() - start_time
                    rate = completed_comparisons / elapsed
                    eta = (total_comparisons - completed_comparisons) / rate
                    
                    logging.info(f"Progress: {progress:.1f}% "
                                f"({completed_comparisons:,}/{total_comparisons:,}) "
                                f"- {significant_pairs:,} significant pairs - "
                                f"ETA: {eta:.1f}s")
        
        # Add diagonal (self-similarities = 1.0)
        for i in range(n_samples):
            similarity_matrix[i, i] = 1.0
        
        # Convert to CSR for efficient operations
        similarity_matrix = similarity_matrix.tocsr()
        
        # Compute statistics
        elapsed_time = time.time() - start_time
        sparsity = 1 - similarity_matrix.nnz / (n_samples * n_samples)
        
        self.computation_stats = {
            'method': 'exact_sparse',
            'total_comparisons': total_comparisons,
            'significant_pairs': significant_pairs,
            'sparsity': sparsity,
            'computation_time': elapsed_time,
            'pairs_per_second': total_comparisons / elapsed_time,
            'memory_saved_vs_dense': f"{sparsity:.1%}"
        }
        
        logging.info(f"Sparse computation complete: {sparsity:.1%} sparsity, "
                    f"{significant_pairs:,} significant pairs, "
                    f"{elapsed_time:.1f}s")
        
        return similarity_matrix
    
    def _compute_lsh_similarities(self, profiles: Dict[str, Dict[str, float]], 
                                method: str) -> csr_matrix:
        """Use Locality-Sensitive Hashing for approximate similarities."""
        if not MINHASH_AVAILABLE:
            logging.warning("LSH requested but datasketch not available, falling back to exact")
            return self._compute_exact_sparse_similarities(profiles, method)
        
        n_samples = len(self.sample_names)
        start_time = time.time()
        
        logging.info("Creating MinHash signatures...")
        
        # Create MinHash signatures for all samples
        minhashes = {}
        num_perm = 128  # Number of permutations for MinHash
        
        for sample_name, profile in profiles.items():
            m = MinHash(num_perm=num_perm)
            for kmer in profile.keys():
                m.update(str(kmer).encode('utf8'))
            minhashes[sample_name] = m
        
        # Build LSH index
        logging.info("Building LSH index...")
        lsh = MinHashLSH(threshold=self.similarity_threshold, num_perm=num_perm)
        
        for i, (sample_name, minhash) in enumerate(minhashes.items()):
            lsh.insert(f"{i}_{sample_name}", minhash)
        
        # Query for similar pairs
        logging.info("Querying for similar pairs...")
        similarity_matrix = lil_matrix((n_samples, n_samples))
        
        total_pairs_found = 0
        
        for i, (sample_name, minhash) in enumerate(minhashes.items()):
            similar_samples = lsh.query(minhash)
            
            for similar_sample in similar_samples:
                j = int(similar_sample.split('_')[0])
                if i != j:  # Skip self-similarity
                    # Estimate similarity using MinHash
                    sim = minhashes[self.sample_names[i]].jaccard(minhashes[self.sample_names[j]])
                    
                    if sim > self.similarity_threshold:
                        similarity_matrix[i, j] = sim
                        total_pairs_found += 1
        
        # Add diagonal (self-similarities = 1.0)
        for i in range(n_samples):
            similarity_matrix[i, i] = 1.0
        
        # Convert to CSR
        similarity_matrix = similarity_matrix.tocsr()
        
        # Compute statistics
        elapsed_time = time.time() - start_time
        sparsity = 1 - similarity_matrix.nnz / (n_samples * n_samples)
        
        self.computation_stats = {
            'method': 'lsh_approximation',
            'total_samples': n_samples,
            'significant_pairs': total_pairs_found,
            'sparsity': sparsity,
            'computation_time': elapsed_time,
            'num_permutations': num_perm,
            'memory_saved_vs_dense': f"{sparsity:.1%}"
        }
        
        logging.info(f"LSH computation complete: {sparsity:.1%} sparsity, "
                    f"{total_pairs_found:,} pairs found, "
                    f"{elapsed_time:.1f}s")
        
        return similarity_matrix
    
    def _jaccard_similarity(self, profile1: Dict[str, float], 
                           profile2: Dict[str, float]) -> float:
        """Jaccard similarity between k-mer profiles."""
        set1, set2 = set(profile1.keys()), set(profile2.keys())
        intersection = len(set1 & set2)
        union = len(set1 | set2)
        return intersection / union if union > 0 else 0
    
    def _weighted_jaccard_similarity(self, profile1: Dict[str, float], 
                                   profile2: Dict[str, float]) -> float:
        """Weighted Jaccard similarity considering frequencies."""
        all_kmers = set(profile1.keys()) | set(profile2.keys())
        
        min_sum = sum(min(profile1.get(k, 0), profile2.get(k, 0)) for k in all_kmers)
        max_sum = sum(max(profile1.get(k, 0), profile2.get(k, 0)) for k in all_kmers)
        
        return min_sum / max_sum if max_sum > 0 else 0
    
    def _cosine_similarity(self, profile1: Dict[str, float], 
                          profile2: Dict[str, float]) -> float:
        """Cosine similarity between frequency vectors."""
        import math
        
        all_kmers = set(profile1.keys()) | set(profile2.keys())
        
        dot_product = sum(profile1.get(k, 0) * profile2.get(k, 0) for k in all_kmers)
        norm1 = math.sqrt(sum(profile1.get(k, 0)**2 for k in all_kmers))
        norm2 = math.sqrt(sum(profile2.get(k, 0)**2 for k in all_kmers))
        
        return dot_product / (norm1 * norm2) if norm1 > 0 and norm2 > 0 else 0
    
    def _braycurtis_distance(self, profile1: Dict[str, float], 
                            profile2: Dict[str, float]) -> float:
        """Bray-Curtis distance between frequency vectors."""
        all_kmers = set(profile1.keys()) | set(profile2.keys())
        
        numerator = sum(abs(profile1.get(k, 0) - profile2.get(k, 0)) for k in all_kmers)
        denominator = sum(profile1.get(k, 0) + profile2.get(k, 0) for k in all_kmers)
        
        return numerator / denominator if denominator > 0 else 0
    
    def _euclidean_distance(self, profile1: Dict[str, float], 
                           profile2: Dict[str, float]) -> float:
        """Euclidean distance between frequency vectors."""
        import math
        
        all_kmers = set(profile1.keys()) | set(profile2.keys())
        
        squared_diffs = sum((profile1.get(k, 0) - profile2.get(k, 0))**2 for k in all_kmers)
        
        return math.sqrt(squared_diffs)
    
    def get_most_similar_samples(self, sample_name: str, 
                               top_n: int = 10) -> List[Tuple[str, float]]:
        """
        Get most similar samples to a given sample.
        
        Args:
            sample_name: Name of the query sample
            top_n: Number of most similar samples to return
            
        Returns:
            List of (sample_name, similarity_score) tuples
        """
        if self.similarity_matrix is None:
            raise ValueError("Must compute similarities first")
        
        if sample_name not in self.sample_names:
            raise ValueError(f"Sample {sample_name} not found")
        
        sample_idx = self.sample_names.index(sample_name)
        similarities = self.similarity_matrix[sample_idx, :].toarray().flatten()
        
        # Get indices of most similar samples (excluding self)
        similar_indices = np.argsort(similarities)[::-1][1:top_n+1]
        
        results = [
            (self.sample_names[idx], similarities[idx])
            for idx in similar_indices
            if similarities[idx] > 0  # Only return non-zero similarities
        ]
        
        return results
    
    def find_clusters(self, min_cluster_size: int = 2, 
                     min_similarity: float = None) -> Dict[int, List[str]]:
        """
        Find clusters of similar samples using connected components.
        
        Args:
            min_cluster_size: Minimum number of samples per cluster
            min_similarity: Minimum similarity for clustering (uses threshold if None)
            
        Returns:
            Dictionary mapping cluster_id -> list of sample names
        """
        if self.similarity_matrix is None:
            raise ValueError("Must compute similarities first")
        
        threshold = min_similarity or self.similarity_threshold
        
        # Create adjacency matrix for clustering
        adjacency = self.similarity_matrix >= threshold
        
        # Find connected components
        from scipy.sparse.csgraph import connected_components
        n_components, labels = connected_components(adjacency, directed=False)
        
        # Group samples by cluster
        clusters = defaultdict(list)
        for i, cluster_id in enumerate(labels):
            clusters[cluster_id].append(self.sample_names[i])
        
        # Filter by minimum cluster size
        filtered_clusters = {
            cluster_id: samples 
            for cluster_id, samples in clusters.items()
            if len(samples) >= min_cluster_size
        }
        
        logging.info(f"Found {len(filtered_clusters)} clusters "
                    f"(min_size={min_cluster_size}, threshold={threshold:.3f})")
        
        return filtered_clusters
    
    def compute_summary_statistics(self) -> Dict[str, Any]:
        """Compute summary statistics for the similarity analysis."""
        if self.similarity_matrix is None:
            raise ValueError("Must compute similarities first")
        
        # Extract non-zero similarities (excluding diagonal)
        similarities = self.similarity_matrix.data
        non_diagonal_similarities = similarities[similarities < 1.0]
        
        n_samples = len(self.sample_names)
        n_significant_pairs = len(non_diagonal_similarities)
        total_possible_pairs = n_samples * (n_samples - 1) // 2
        
        stats = {
            'n_samples': n_samples,
            'n_significant_pairs': n_significant_pairs,
            'total_possible_pairs': total_possible_pairs,
            'sparsity': 1 - self.similarity_matrix.nnz / (n_samples * n_samples),
            'mean_similarity': np.mean(non_diagonal_similarities) if len(non_diagonal_similarities) > 0 else 0,
            'std_similarity': np.std(non_diagonal_similarities) if len(non_diagonal_similarities) > 0 else 0,
            'min_similarity': np.min(non_diagonal_similarities) if len(non_diagonal_similarities) > 0 else 0,
            'max_similarity': np.max(non_diagonal_similarities) if len(non_diagonal_similarities) > 0 else 0,
            'median_similarity': np.median(non_diagonal_similarities) if len(non_diagonal_similarities) > 0 else 0,
            'similarity_threshold': self.similarity_threshold,
            'computation_stats': self.computation_stats
        }
        
        logging.info(f"Similarity statistics: {n_significant_pairs:,} significant pairs, "
                    f"mean similarity: {stats['mean_similarity']:.3f}")
        
        return stats
    
    def save_similarity_matrix(self, output_path: str, format: str = 'csv'):
        """
        Save similarity matrix to file.
        
        Args:
            output_path: Path to save the matrix
            format: Output format ('csv', 'npz', 'mtx')
        """
        if self.similarity_matrix is None:
            raise ValueError("Must compute similarities first")
        
        if format == 'csv':
            # Convert to dense for CSV (only for small matrices)
            if len(self.sample_names) > 1000:
                logging.warning("Large matrix - CSV export may use significant memory")
            
            dense_matrix = self.similarity_matrix.toarray()
            import pandas as pd
            df = pd.DataFrame(dense_matrix, index=self.sample_names, columns=self.sample_names)
            df.to_csv(output_path)
            
        elif format == 'npz':
            # Save as sparse matrix
            from scipy.sparse import save_npz
            save_npz(output_path, self.similarity_matrix)
            
            # Also save sample names
            sample_names_path = output_path.replace('.npz', '_samples.txt')
            with open(sample_names_path, 'w') as f:
                for name in self.sample_names:
                    f.write(f"{name}\n")
                    
        elif format == 'mtx':
            # Matrix Market format
            from scipy.io import mmwrite
            mmwrite(output_path, self.similarity_matrix)
            
        else:
            raise ValueError(f"Unknown format: {format}")
        
        logging.info(f"Similarity matrix saved to {output_path} (format: {format})")


def benchmark_sparse_vs_dense(profiles: Dict[str, Dict[str, float]], 
                             similarity_threshold: float = 0.1) -> Dict[str, Any]:
    """
    Benchmark sparse vs dense similarity computation.
    
    Args:
        profiles: Sample profiles for testing
        similarity_threshold: Threshold for sparse computation
        
    Returns:
        Dictionary with benchmark results
    """
    import time
    import psutil
    import os
    
    def get_memory_usage():
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024  # MB
    
    n_samples = len(profiles)
    results = {}
    
    logging.info(f"Benchmarking with {n_samples} samples, threshold={similarity_threshold}")
    
    # Test sparse method
    logging.info("Testing sparse similarity computation...")
    start_memory = get_memory_usage()
    start_time = time.perf_counter()
    
    sparse_analyzer = SparseSimilarityAnalyzer(similarity_threshold=similarity_threshold)
    sparse_matrix, sample_names = sparse_analyzer.compute_similarities(profiles, method='jaccard')
    
    sparse_time = time.perf_counter() - start_time
    sparse_memory = get_memory_usage() - start_memory
    sparse_stats = sparse_analyzer.compute_summary_statistics()
    
    results['sparse'] = {
        'time': sparse_time,
        'memory_mb': sparse_memory,
        'sparsity': sparse_stats['sparsity'],
        'significant_pairs': sparse_stats['n_significant_pairs'],
        'matrix_size': sparse_matrix.nnz
    }
    
    # Estimate dense method (only compute if small enough)
    if n_samples <= 100:  # Only test dense for small datasets
        logging.info("Testing dense similarity computation (small dataset)...")
        
        start_memory = get_memory_usage()
        start_time = time.perf_counter()
        
        # Simulate dense computation
        dense_matrix = np.zeros((n_samples, n_samples))
        for i in range(n_samples):
            for j in range(i, n_samples):
                if i == j:
                    dense_matrix[i, j] = 1.0
                else:
                    # Use same similarity computation
                    sim = sparse_analyzer._jaccard_similarity(
                        profiles[sample_names[i]], 
                        profiles[sample_names[j]]
                    )
                    dense_matrix[i, j] = sim
                    dense_matrix[j, i] = sim
        
        dense_time = time.perf_counter() - start_time
        dense_memory = get_memory_usage() - start_memory
        
        results['dense'] = {
            'time': dense_time,
            'memory_mb': dense_memory,
            'sparsity': 0.0,  # Dense matrix
            'significant_pairs': n_samples * (n_samples - 1) // 2,
            'matrix_size': n_samples * n_samples * 8  # Estimated bytes
        }
        
        # Compute speedup and memory savings
        results['comparison'] = {
            'time_speedup': dense_time / sparse_time if sparse_time > 0 else float('inf'),
            'memory_reduction': dense_memory / sparse_memory if sparse_memory > 0 else float('inf'),
            'sparsity_benefit': sparse_stats['sparsity']
        }
    
    else:
        # Estimate dense requirements for large datasets
        estimated_dense_memory = (n_samples * n_samples * 8) / (1024 * 1024)  # MB
        estimated_dense_time = sparse_time * (1 / sparse_stats['sparsity']) if sparse_stats['sparsity'] < 1 else sparse_time * 100
        
        results['dense_estimated'] = {
            'time_estimated': estimated_dense_time,
            'memory_mb_estimated': estimated_dense_memory,
            'feasible': estimated_dense_memory < 8000  # Less than 8GB
        }
        
        results['comparison'] = {
            'time_speedup_estimated': estimated_dense_time / sparse_time,
            'memory_reduction_estimated': estimated_dense_memory / sparse_memory,
            'sparsity_benefit': sparse_stats['sparsity']
        }
    
    return results


if __name__ == "__main__":
    # Quick test
    logging.basicConfig(level=logging.INFO)
    
    # Create test profiles
    test_profiles = {}
    for i in range(10):
        # Create overlapping k-mer profiles
        profile = {}
        
        # Common k-mers (high similarity)
        for j in range(50):
            kmer = f"ATCG{'A' * j:0>17}"  # Pad to 21 chars
            profile[kmer] = 0.01 + (j % 3) * 0.005
        
        # Sample-specific k-mers
        for j in range(20):
            kmer = f"GCTA{i:02d}{'T' * j:0>13}"  # Sample-specific
            profile[kmer] = 0.005 + j * 0.001
        
        test_profiles[f"sample_{i:02d}"] = profile
    
    # Test sparse analyzer
    analyzer = SparseSimilarityAnalyzer(similarity_threshold=0.2)
    similarity_matrix, sample_names = analyzer.compute_similarities(test_profiles)
    
    print(f"Created sparse matrix: {similarity_matrix.shape}")
    print(f"Sparsity: {1 - similarity_matrix.nnz / (similarity_matrix.shape[0] * similarity_matrix.shape[1]):.1%}")
    
    # Test clustering
    clusters = analyzer.find_clusters(min_cluster_size=2)
    print(f"Found {len(clusters)} clusters")
    
    # Test statistics
    stats = analyzer.compute_summary_statistics()
    print(f"Mean similarity: {stats['mean_similarity']:.3f}")
    print(f"Significant pairs: {stats['n_significant_pairs']}")