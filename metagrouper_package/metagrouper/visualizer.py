"""
Visualization functionality for metagenomic analysis results.

This module contains the Visualizer class for generating plots and visualizations
of sample relationships, distance matrices, and dimensionality reduction results.
"""

import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List
from sklearn.decomposition import PCA


class Visualizer:
    """Generate visualizations for sample relationships."""

    def __init__(self, sample_names: List[str]):
        self.sample_names = sample_names

    def plot_distance_heatmap(self, distance_matrix: np.ndarray, output_path: str):
        """Plot distance matrix as heatmap."""
        plt.figure(figsize=(10, 8))

        # Create DataFrame for better labeling
        df = pd.DataFrame(
            distance_matrix, index=self.sample_names, columns=self.sample_names
        )

        sns.heatmap(df, annot=True, cmap="viridis", fmt=".3f")
        plt.title("Sample Distance Matrix")
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"Distance heatmap saved to {output_path}")

    def plot_pca(self, pca_result: np.ndarray, pca: PCA, output_path: str):
        """Plot PCA results."""
        plt.figure(figsize=(10, 8))

        plt.scatter(pca_result[:, 0], pca_result[:, 1], s=100, alpha=0.7)

        # Add sample labels
        for i, sample in enumerate(self.sample_names):
            plt.annotate(
                sample,
                (pca_result[i, 0], pca_result[i, 1]),
                xytext=(5, 5),
                textcoords="offset points",
            )

        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
        plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
        plt.title("PCA of K-mer Profiles")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"PCA plot saved to {output_path}")

    def plot_mds(self, mds_result: np.ndarray, output_path: str):
        """Plot MDS results."""
        plt.figure(figsize=(10, 8))

        plt.scatter(mds_result[:, 0], mds_result[:, 1], s=100, alpha=0.7)

        # Add sample labels
        for i, sample in enumerate(self.sample_names):
            plt.annotate(
                sample,
                (mds_result[i, 0], mds_result[i, 1]),
                xytext=(5, 5),
                textcoords="offset points",
            )

        plt.xlabel("MDS Dimension 1")
        plt.ylabel("MDS Dimension 2")
        plt.title("MDS of Sample Distances")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"MDS plot saved to {output_path}")

    def plot_sample_overview(self, distance_matrix: np.ndarray, pca_result: np.ndarray, 
                           mds_result: np.ndarray, pca: PCA, output_path: str):
        """Create a comprehensive overview plot with multiple visualizations."""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Distance matrix heatmap
        df = pd.DataFrame(distance_matrix, index=self.sample_names, columns=self.sample_names)
        sns.heatmap(df, annot=True, cmap="viridis", fmt=".3f", ax=axes[0, 0], cbar_kws={'shrink': 0.8})
        axes[0, 0].set_title("Sample Distance Matrix")
        
        # PCA plot
        axes[0, 1].scatter(pca_result[:, 0], pca_result[:, 1], s=100, alpha=0.7)
        for i, sample in enumerate(self.sample_names):
            axes[0, 1].annotate(sample, (pca_result[i, 0], pca_result[i, 1]), 
                              xytext=(5, 5), textcoords="offset points")
        axes[0, 1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
        axes[0, 1].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
        axes[0, 1].set_title("PCA of K-mer Profiles")
        axes[0, 1].grid(True, alpha=0.3)
        
        # MDS plot
        axes[1, 0].scatter(mds_result[:, 0], mds_result[:, 1], s=100, alpha=0.7)
        for i, sample in enumerate(self.sample_names):
            axes[1, 0].annotate(sample, (mds_result[i, 0], mds_result[i, 1]), 
                              xytext=(5, 5), textcoords="offset points")
        axes[1, 0].set_xlabel("MDS Dimension 1")
        axes[1, 0].set_ylabel("MDS Dimension 2")
        axes[1, 0].set_title("MDS of Sample Distances")
        axes[1, 0].grid(True, alpha=0.3)
        
        # Sample statistics
        axes[1, 1].axis('off')
        stats_text = f"""Sample Analysis Overview
        
Number of samples: {len(self.sample_names)}
Distance matrix shape: {distance_matrix.shape}
Mean distance: {np.mean(distance_matrix):.3f}
Std distance: {np.std(distance_matrix):.3f}
Min distance: {np.min(distance_matrix[distance_matrix > 0]):.3f}
Max distance: {np.max(distance_matrix):.3f}

PCA Variance Explained:
PC1: {pca.explained_variance_ratio_[0]:.1%}
PC2: {pca.explained_variance_ratio_[1]:.1%}
Total: {sum(pca.explained_variance_ratio_):.1%}
        """
        axes[1, 1].text(0.1, 0.9, stats_text, transform=axes[1, 1].transAxes, 
                        fontsize=10, verticalalignment='top', fontfamily='monospace')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        
        logging.info(f"Sample overview plot saved to {output_path}")

    def plot_clustering_dendrogram(self, distance_matrix: np.ndarray, output_path: str):
        """Plot hierarchical clustering dendrogram."""
        from scipy.cluster.hierarchy import dendrogram, linkage
        
        plt.figure(figsize=(12, 8))
        
        # Perform hierarchical clustering
        linkage_matrix = linkage(distance_matrix, method='ward')
        
        # Create dendrogram
        dendrogram(linkage_matrix, labels=self.sample_names, orientation='top')
        plt.title("Hierarchical Clustering of Samples")
        plt.xlabel("Samples")
        plt.ylabel("Distance")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        
        logging.info(f"Clustering dendrogram saved to {output_path}")

    def plot_distance_distribution(self, distance_matrix: np.ndarray, output_path: str):
        """Plot distribution of pairwise distances."""
        plt.figure(figsize=(10, 6))
        
        # Extract upper triangle (excluding diagonal) for pairwise distances
        distances = distance_matrix[np.triu_indices_from(distance_matrix, k=1)]
        
        plt.hist(distances, bins=20, alpha=0.7, edgecolor='black')
        plt.axvline(np.mean(distances), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(distances):.3f}')
        plt.axvline(np.median(distances), color='green', linestyle='--', 
                   label=f'Median: {np.median(distances):.3f}')
        
        plt.xlabel("Pairwise Distance")
        plt.ylabel("Frequency")
        plt.title("Distribution of Pairwise Sample Distances")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        
        logging.info(f"Distance distribution plot saved to {output_path}")