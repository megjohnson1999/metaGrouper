#!/usr/bin/env python3
"""
Compression-based similarity analyzer for MetaGrouper.

This module provides analysis tools that work with compression-based
distance matrices instead of k-mer profiles.
"""

import numpy as np
import pandas as pd
import logging
from typing import Dict, List, Optional, Tuple, Any
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


class CompressionAnalyzer:
    """
    Analyze compression-based similarity matrices and perform dimensionality reduction.
    
    This class replaces the k-mer profile based SimilarityAnalyzer with methods
    that work directly on compression distance matrices.
    """
    
    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str]):
        """
        Initialize analyzer with compression-based distance matrix.
        
        Args:
            distance_matrix: NCD distance matrix of shape (n_samples, n_samples)
            sample_names: List of sample names corresponding to matrix rows/columns
        """
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names
        self.n_samples = len(sample_names)
        
        # Validate inputs
        if distance_matrix.shape != (self.n_samples, self.n_samples):
            raise ValueError(f"Distance matrix shape {distance_matrix.shape} doesn't match "
                           f"number of samples {self.n_samples}")
        
        # Storage for analysis results
        self.pca_result = None
        self.mds_result = None
        self.pca_model = None
        self.mds_model = None
        
        logging.info(f"Initialized CompressionAnalyzer with {self.n_samples} samples")
    
    def perform_pca(self, n_components: int = 2) -> Tuple[np.ndarray, Any]:
        """
        Perform PCA on the distance matrix.
        
        Note: For distance matrices, we convert to similarity and then to coordinates.
        
        Args:
            n_components: Number of principal components
            
        Returns:
            Tuple of (pca_coordinates, pca_model)
        """
        logging.info(f"Performing PCA with {n_components} components")
        
        # Convert distance matrix to similarity matrix
        # Use double centering (classical MDS approach) to get coordinates
        max_dist = np.max(self.distance_matrix)
        similarity_matrix = max_dist - self.distance_matrix
        
        # Double centering
        n = similarity_matrix.shape[0]
        H = np.eye(n) - np.ones((n, n)) / n
        B = -0.5 * H @ (self.distance_matrix ** 2) @ H
        
        # PCA on the centered matrix
        self.pca_model = PCA(n_components=n_components)
        self.pca_result = self.pca_model.fit_transform(B)
        
        logging.info(f"PCA explained variance ratio: {self.pca_model.explained_variance_ratio_}")
        
        return self.pca_result, self.pca_model
    
    def perform_mds(self, n_components: int = 2, metric: bool = True) -> np.ndarray:
        """
        Perform Multidimensional Scaling on the distance matrix.
        
        Args:
            n_components: Number of MDS dimensions
            metric: Whether to use metric MDS (default) or non-metric
            
        Returns:
            np.ndarray: MDS coordinates
        """
        logging.info(f"Performing {'metric' if metric else 'non-metric'} MDS "
                    f"with {n_components} components")
        
        self.mds_model = MDS(
            n_components=n_components,
            dissimilarity='precomputed',
            metric=metric,
            random_state=42
        )
        
        self.mds_result = self.mds_model.fit_transform(self.distance_matrix)
        
        if hasattr(self.mds_model, 'stress_'):
            logging.info(f"MDS stress: {self.mds_model.stress_:.4f}")
        
        return self.mds_result
    
    def perform_clustering(self, 
                          n_clusters_range: Tuple[int, int] = (2, 8),
                          methods: List[str] = ['kmeans', 'hierarchical']) -> Dict[str, Any]:
        """
        Perform clustering on the distance matrix using multiple methods.
        
        Args:
            n_clusters_range: Range of cluster numbers to test
            methods: Clustering methods to use
            
        Returns:
            Dict with clustering results
        """
        logging.info(f"Performing clustering with methods: {methods}")
        
        results = {}
        min_clusters, max_clusters = n_clusters_range
        
        # For clustering, we need coordinate representation
        if self.mds_result is None:
            self.perform_mds()
        
        coordinates = self.mds_result
        
        for method in methods:
            results[method] = {}
            best_score = -1
            best_k = None
            
            for k in range(min_clusters, max_clusters + 1):
                if k >= self.n_samples:
                    continue
                
                try:
                    if method == 'kmeans':
                        clusterer = KMeans(n_clusters=k, random_state=42, n_init=10)
                        labels = clusterer.fit_predict(coordinates)
                    
                    elif method == 'hierarchical':
                        # Use distance matrix directly for hierarchical clustering
                        condensed_dist = squareform(self.distance_matrix)
                        linkage_matrix = linkage(condensed_dist, method='average')
                        labels = fcluster(linkage_matrix, k, criterion='maxclust') - 1
                    
                    else:
                        logging.warning(f"Unknown clustering method: {method}")
                        continue
                    
                    # Calculate silhouette score
                    if len(np.unique(labels)) > 1:
                        score = silhouette_score(coordinates, labels)
                        
                        results[method][k] = {
                            'labels': labels,
                            'silhouette_score': score,
                            'n_clusters': k
                        }
                        
                        if score > best_score:
                            best_score = score
                            best_k = k
                
                except Exception as e:
                    logging.warning(f"Clustering failed for {method}, k={k}: {e}")
                    continue
            
            # Store best result
            if best_k is not None:
                results[method]['optimal'] = {
                    'n_clusters': best_k,
                    'silhouette_score': best_score,
                    'labels': results[method][best_k]['labels']
                }
                
                logging.info(f"{method.capitalize()} clustering: "
                           f"optimal k={best_k}, silhouette={best_score:.3f}")
        
        return results
    
    def get_sample_similarities(self, sample_name: str, top_n: int = 5) -> List[Tuple[str, float]]:
        """
        Get most similar samples to a given sample.
        
        Args:
            sample_name: Name of the query sample
            top_n: Number of most similar samples to return
            
        Returns:
            List of (sample_name, similarity_score) tuples
        """
        if sample_name not in self.sample_names:
            raise ValueError(f"Sample {sample_name} not found in dataset")
        
        sample_idx = self.sample_names.index(sample_name)
        distances = self.distance_matrix[sample_idx, :]
        
        # Convert distances to similarities (1 - distance)
        similarities = 1 - distances
        
        # Get indices of most similar samples (excluding self)
        similar_indices = np.argsort(similarities)[::-1][1:top_n+1]
        
        results = [
            (self.sample_names[idx], similarities[idx])
            for idx in similar_indices
        ]
        
        return results
    
    def compute_summary_statistics(self) -> Dict[str, float]:
        """
        Compute summary statistics for the distance matrix.
        
        Returns:
            Dict with summary statistics
        """
        # Extract upper triangle (excluding diagonal)
        upper_triangle = self.distance_matrix[np.triu_indices_from(self.distance_matrix, k=1)]
        
        stats = {
            'mean_distance': np.mean(upper_triangle),
            'median_distance': np.median(upper_triangle),
            'std_distance': np.std(upper_triangle),
            'min_distance': np.min(upper_triangle),
            'max_distance': np.max(upper_triangle),
            'n_samples': self.n_samples,
            'n_comparisons': len(upper_triangle)
        }
        
        logging.info(f"Distance matrix statistics: "
                    f"mean={stats['mean_distance']:.3f}, "
                    f"std={stats['std_distance']:.3f}, "
                    f"range=[{stats['min_distance']:.3f}, {stats['max_distance']:.3f}]")
        
        return stats
    
    def identify_outliers(self, method: str = 'isolation', contamination: float = 0.1) -> List[str]:
        """
        Identify outlier samples based on compression distances.
        
        Args:
            method: Outlier detection method ('isolation', 'distance_based')
            contamination: Expected proportion of outliers
            
        Returns:
            List of outlier sample names
        """
        logging.info(f"Identifying outliers using {method} method")
        
        if method == 'distance_based':
            # Samples with consistently high distances to all others
            mean_distances = np.mean(self.distance_matrix, axis=1)
            threshold = np.percentile(mean_distances, (1 - contamination) * 100)
            outlier_indices = np.where(mean_distances > threshold)[0]
        
        elif method == 'isolation':
            # Use sklearn's IsolationForest on MDS coordinates
            if self.mds_result is None:
                self.perform_mds()
            
            from sklearn.ensemble import IsolationForest
            
            iso_forest = IsolationForest(contamination=contamination, random_state=42)
            outlier_labels = iso_forest.fit_predict(self.mds_result)
            outlier_indices = np.where(outlier_labels == -1)[0]
        
        else:
            raise ValueError(f"Unknown outlier detection method: {method}")
        
        outlier_samples = [self.sample_names[idx] for idx in outlier_indices]
        
        logging.info(f"Identified {len(outlier_samples)} outliers: {outlier_samples}")
        
        return outlier_samples
    
    def export_distance_matrix(self, output_path: str):
        """
        Export distance matrix to CSV file.
        
        Args:
            output_path: Path to save the distance matrix
        """
        df = pd.DataFrame(
            self.distance_matrix,
            index=self.sample_names,
            columns=self.sample_names
        )
        
        df.to_csv(output_path)
        logging.info(f"Distance matrix exported to {output_path}")