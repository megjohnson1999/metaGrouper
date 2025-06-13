"""
Similarity analysis functionality for k-mer profiles.

This module contains the SimilarityAnalyzer class for computing distance matrices
and performing dimensionality reduction on k-mer profiles.
"""

import logging
import numpy as np
from collections import defaultdict
from typing import Dict, Tuple
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.metrics.pairwise import pairwise_distances


class SimilarityAnalyzer:
    """Analyze similarities between k-mer profiles with memory optimization."""

    def __init__(self, profiles: Dict[str, Dict[str, float]], memory_efficient: bool = True):
        self.profiles = profiles
        self.sample_names = list(profiles.keys())
        self.distance_matrix = None
        self.similarity_matrix = None
        self.memory_efficient = memory_efficient

    def compute_distance_matrix(self, metric: str = "braycurtis") -> np.ndarray:
        """Compute pairwise distance matrix between samples with memory optimization."""
        logging.info(f"Computing distance matrix using {metric} metric")
        
        n_samples = len(self.sample_names)
        
        if self.memory_efficient and n_samples > 50:
            # Memory-efficient computation for large datasets
            logging.info("Using memory-efficient distance computation")
            return self._compute_distance_matrix_efficient(metric)
        else:
            # Standard computation for smaller datasets
            return self._compute_distance_matrix_standard(metric)

    def _compute_distance_matrix_standard(self, metric: str) -> np.ndarray:
        """Standard distance matrix computation (loads all data into memory)."""
        # Get all unique k-mers across all samples
        all_kmers = set()
        for profile in self.profiles.values():
            all_kmers.update(profile.keys())
        all_kmers = sorted(list(all_kmers))
        
        logging.info(f"Using {len(all_kmers)} unique k-mers across {len(self.sample_names)} samples")

        # Create feature matrix
        feature_matrix = np.zeros((len(self.sample_names), len(all_kmers)))
        for i, sample in enumerate(self.sample_names):
            for j, kmer in enumerate(all_kmers):
                feature_matrix[i, j] = self.profiles[sample].get(kmer, 0)

        # Compute pairwise distances
        self.distance_matrix = pairwise_distances(feature_matrix, metric=metric)
        self.similarity_matrix = 1 - self.distance_matrix

        return self.distance_matrix

    def _compute_distance_matrix_efficient(self, metric: str) -> np.ndarray:
        """Memory-efficient distance matrix computation for large datasets."""
        from scipy.spatial.distance import pdist, squareform
        from scipy.sparse import csr_matrix
        
        n_samples = len(self.sample_names)
        
        # Get common k-mers (present in multiple samples) to reduce dimensionality
        kmer_counts = defaultdict(int)
        for profile in self.profiles.values():
            for kmer in profile.keys():
                kmer_counts[kmer] += 1
        
        # Keep k-mers present in at least 2 samples or with high frequency
        min_samples = max(2, int(0.1 * n_samples))  # At least 10% of samples
        common_kmers = [kmer for kmer, count in kmer_counts.items() if count >= min_samples]
        
        if not common_kmers:
            # Fallback to all k-mers if no common ones
            common_kmers = list(kmer_counts.keys())
        
        logging.info(f"Using {len(common_kmers)} common k-mers (from {len(kmer_counts)} total)")
        
        # Build sparse feature matrix
        data, row_indices, col_indices = [], [], []
        for i, sample in enumerate(self.sample_names):
            profile = self.profiles[sample]
            for j, kmer in enumerate(common_kmers):
                if kmer in profile and profile[kmer] > 0:
                    data.append(profile[kmer])
                    row_indices.append(i)
                    col_indices.append(j)
        
        sparse_matrix = csr_matrix((data, (row_indices, col_indices)), 
                                 shape=(n_samples, len(common_kmers)))
        
        # Convert to dense for distance computation (only if manageable size)
        if sparse_matrix.nnz < 1000000:  # Less than 1M non-zero elements
            feature_matrix = sparse_matrix.toarray()
            self.distance_matrix = pairwise_distances(feature_matrix, metric=metric)
        else:
            # For very large datasets, compute distances chunk by chunk
            logging.info("Computing distances in chunks for very large dataset")
            self.distance_matrix = np.zeros((n_samples, n_samples))
            
            chunk_size = 10
            for i in range(0, n_samples, chunk_size):
                end_i = min(i + chunk_size, n_samples)
                chunk_i = sparse_matrix[i:end_i].toarray()
                
                for j in range(i, n_samples, chunk_size):
                    end_j = min(j + chunk_size, n_samples)
                    chunk_j = sparse_matrix[j:end_j].toarray()
                    
                    # Compute distances for this chunk
                    chunk_distances = pairwise_distances(chunk_i, chunk_j, metric=metric)
                    
                    # Store in matrix
                    self.distance_matrix[i:end_i, j:end_j] = chunk_distances
                    if i != j:  # Fill symmetric part
                        self.distance_matrix[j:end_j, i:end_i] = chunk_distances.T
        
        self.similarity_matrix = 1 - self.distance_matrix
        return self.distance_matrix

    def perform_pca(self, n_components: int = 2) -> Tuple[np.ndarray, PCA]:
        """Perform PCA on k-mer profiles."""
        logging.info("Performing PCA analysis")

        # Get all unique k-mers
        all_kmers = set()
        for profile in self.profiles.values():
            all_kmers.update(profile.keys())
        all_kmers = sorted(list(all_kmers))

        # Create feature matrix
        feature_matrix = np.zeros((len(self.sample_names), len(all_kmers)))
        for i, sample in enumerate(self.sample_names):
            for j, kmer in enumerate(all_kmers):
                feature_matrix[i, j] = self.profiles[sample].get(kmer, 0)

        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(feature_matrix)

        logging.info(f"PCA explained variance ratio: {pca.explained_variance_ratio_}")
        return pca_result, pca

    def perform_mds(self, n_components: int = 2) -> np.ndarray:
        """Perform MDS on distance matrix."""
        logging.info("Performing MDS analysis")

        if self.distance_matrix is None:
            self.compute_distance_matrix()

        mds = MDS(
            n_components=n_components, dissimilarity="precomputed", random_state=42
        )
        mds_result = mds.fit_transform(self.distance_matrix)

        return mds_result