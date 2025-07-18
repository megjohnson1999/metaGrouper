"""
Similarity analysis functionality for k-mer profiles.

This module contains the SimilarityAnalyzer class for computing distance matrices
and performing dimensionality reduction on k-mer profiles.
"""

import logging
import numpy as np
from collections import defaultdict
from typing import Dict, Tuple, Optional
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.manifold import MDS
from sklearn.metrics.pairwise import pairwise_distances


class SimilarityAnalyzer:
    """Analyze similarities between k-mer profiles with memory optimization."""

    def __init__(self, profiles: Dict[str, Dict[str, float]], memory_efficient: bool = True, 
                 prevalence_threshold: float = 0.1):
        self.profiles = profiles
        self.sample_names = list(profiles.keys())
        self.distance_matrix = None
        self.similarity_matrix = None
        self.memory_efficient = memory_efficient
        self.prevalence_threshold = prevalence_threshold
        
        # Auto-adjust prevalence threshold for very large datasets
        n_samples = len(self.sample_names)
        if n_samples > 500 and prevalence_threshold < 0.2:
            self.prevalence_threshold = 0.2  # 20% for large datasets
            logging.info(f"Auto-adjusted prevalence threshold to {self.prevalence_threshold:.1%} for large dataset ({n_samples} samples)")
        elif n_samples > 1000 and prevalence_threshold < 0.3:
            self.prevalence_threshold = 0.3  # 30% for very large datasets
            logging.info(f"Auto-adjusted prevalence threshold to {self.prevalence_threshold:.1%} for very large dataset ({n_samples} samples)")

    def compute_distance_matrix(self, metric: str = "jaccard") -> np.ndarray:
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
        
        # Keep k-mers present in at least prevalence_threshold of samples (default 10% for memory efficiency)
        min_samples = max(2, int(self.prevalence_threshold * n_samples))
        common_kmers = [kmer for kmer, count in kmer_counts.items() if count >= min_samples]
        
        if not common_kmers:
            # Fallback to all k-mers if no common ones
            common_kmers = list(kmer_counts.keys())
        
        reduction_percent = 100 * (1 - len(common_kmers) / len(kmer_counts)) if kmer_counts else 0
        logging.info(f"K-mer filtering: {len(kmer_counts):,} total → {len(common_kmers):,} retained "
                    f"({reduction_percent:.1f}% reduction, min_samples={min_samples})")
        
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

    def perform_pca(self, n_components: int = 2, use_sparse: bool = None) -> Tuple[np.ndarray, PCA]:
        """Perform PCA on k-mer profiles with automatic sparse detection."""
        logging.info("Performing PCA analysis")

        # Get all unique k-mers
        all_kmers = set()
        for profile in self.profiles.values():
            all_kmers.update(profile.keys())
        all_kmers = sorted(list(all_kmers))

        n_samples = len(self.sample_names)
        n_features = len(all_kmers)
        
        # Auto-detect whether to use sparse (TruncatedSVD) based on size
        if use_sparse is None:
            use_sparse = n_features > 10000 or (n_samples * n_features) > 1e7
            
        if use_sparse:
            logging.info(f"Using TruncatedSVD for sparse PCA ({n_features:,} features)")
            return self._perform_sparse_pca(n_components, all_kmers)
        else:
            # Create dense feature matrix
            feature_matrix = np.zeros((n_samples, n_features))
            for i, sample in enumerate(self.sample_names):
                for j, kmer in enumerate(all_kmers):
                    feature_matrix[i, j] = self.profiles[sample].get(kmer, 0)

            pca = PCA(n_components=n_components)
            pca_result = pca.fit_transform(feature_matrix)

            logging.info(f"PCA explained variance ratio: {pca.explained_variance_ratio_}")
            return pca_result, pca
    
    def _perform_sparse_pca(self, n_components: int, all_kmers: list) -> Tuple[np.ndarray, TruncatedSVD]:
        """Perform PCA using TruncatedSVD on sparse matrix."""
        from scipy.sparse import csr_matrix
        
        # Build sparse feature matrix
        data, row_indices, col_indices = [], [], []
        kmer_to_idx = {kmer: idx for idx, kmer in enumerate(all_kmers)}
        
        for i, sample in enumerate(self.sample_names):
            profile = self.profiles[sample]
            for kmer, count in profile.items():
                if kmer in kmer_to_idx and count > 0:
                    data.append(count)
                    row_indices.append(i)
                    col_indices.append(kmer_to_idx[kmer])
        
        sparse_matrix = csr_matrix((data, (row_indices, col_indices)), 
                                 shape=(len(self.sample_names), len(all_kmers)))
        
        # Use TruncatedSVD for sparse matrices
        svd = TruncatedSVD(n_components=n_components, random_state=42)
        pca_result = svd.fit_transform(sparse_matrix)
        
        # Make it compatible with PCA interface
        svd.explained_variance_ratio_ = svd.explained_variance_ratio_
        
        logging.info(f"TruncatedSVD explained variance ratio: {svd.explained_variance_ratio_}")
        logging.info(f"Sparse matrix: {sparse_matrix.nnz:,} non-zero elements out of {sparse_matrix.shape[0] * sparse_matrix.shape[1]:,}")
        
        return pca_result, svd

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