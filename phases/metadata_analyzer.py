#!/usr/bin/env python3
"""
MetaGrouper Phase 2: Metadata Variable Analysis

This module adds metadata analysis capabilities including PERMANOVA,
variable importance ranking, and clustering analysis.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
import logging
from pathlib import Path
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import LabelEncoder
from sklearn.cluster import KMeans, AgglomerativeClustering, DBSCAN
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
import warnings

warnings.filterwarnings("ignore")


class PermanovaAnalyzer:
    """PERMANOVA (Permutational Multivariate Analysis of Variance) implementation."""

    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str]):
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names
        self.n_samples = len(sample_names)

    def _calculate_sum_of_squares(
        self, distance_matrix: np.ndarray, groups: np.ndarray
    ) -> Tuple[float, float]:
        """Calculate within-group and total sum of squares."""
        n = distance_matrix.shape[0]

        # Total sum of squares
        total_ss = np.sum(distance_matrix**2) / n

        # Within-group sum of squares
        within_ss = 0
        unique_groups = np.unique(groups)

        for group in unique_groups:
            group_indices = np.where(groups == group)[0]
            if len(group_indices) > 1:
                group_distances = distance_matrix[np.ix_(group_indices, group_indices)]
                within_ss += np.sum(group_distances**2) / len(group_indices)

        return within_ss, total_ss

    def permanova_test(
        self, metadata_variable: np.ndarray, n_permutations: int = 999
    ) -> Dict[str, float]:
        """Perform PERMANOVA test for a single metadata variable."""

        # Remove samples with missing metadata
        valid_indices = ~pd.isna(metadata_variable)
        if not np.any(valid_indices):
            return {"f_statistic": np.nan, "p_value": np.nan, "r_squared": np.nan}

        valid_distance_matrix = self.distance_matrix[
            np.ix_(valid_indices, valid_indices)
        ]
        valid_groups = metadata_variable[valid_indices]

        # Calculate observed F-statistic
        within_ss, total_ss = self._calculate_sum_of_squares(
            valid_distance_matrix, valid_groups
        )

        if total_ss == 0:
            return {"f_statistic": np.nan, "p_value": np.nan, "r_squared": np.nan}

        between_ss = total_ss - within_ss

        # Degrees of freedom
        n_groups = len(np.unique(valid_groups))
        n_samples = len(valid_groups)
        df_between = n_groups - 1
        df_within = n_samples - n_groups

        if df_between == 0 or df_within == 0:
            return {"f_statistic": np.nan, "p_value": np.nan, "r_squared": np.nan}

        # F-statistic
        f_observed = (between_ss / df_between) / (within_ss / df_within)

        # R-squared (proportion of variation explained)
        r_squared = between_ss / total_ss

        # Permutation test
        f_permuted = []
        for _ in range(n_permutations):
            # Permute group labels
            permuted_groups = np.random.permutation(valid_groups)
            perm_within_ss, perm_total_ss = self._calculate_sum_of_squares(
                valid_distance_matrix, permuted_groups
            )
            perm_between_ss = perm_total_ss - perm_within_ss

            if perm_within_ss > 0 and df_within > 0:
                f_perm = (perm_between_ss / df_between) / (perm_within_ss / df_within)
                f_permuted.append(f_perm)

        # P-value
        if f_permuted:
            p_value = (np.sum(np.array(f_permuted) >= f_observed) + 1) / (
                len(f_permuted) + 1
            )
        else:
            p_value = np.nan

        return {
            "f_statistic": f_observed,
            "p_value": p_value,
            "r_squared": r_squared,
            "n_samples": n_samples,
            "n_groups": n_groups,
        }


class MetadataAnalyzer:
    """Analyze metadata variables and their relationship to sample similarities."""

    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str]):
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names
        self.metadata = None
        self.permanova_results = {}
        self.cluster_results = {}

    def load_metadata(self, metadata_file: str, sample_id_column: str = "sample_id"):
        """Load metadata from CSV/TSV file."""
        logging.info(f"Loading metadata from {metadata_file}")

        # Detect file format
        if metadata_file.endswith(".tsv") or metadata_file.endswith(".txt"):
            self.metadata = pd.read_csv(metadata_file, sep="\t")
        else:
            self.metadata = pd.read_csv(metadata_file)

        # Ensure sample_id column exists
        if sample_id_column not in self.metadata.columns:
            raise ValueError(
                f"Sample ID column '{sample_id_column}' not found in metadata"
            )

        # Align metadata with sample names
        self.metadata = self.metadata.set_index(sample_id_column)
        self.metadata = self.metadata.reindex(self.sample_names)

        logging.info(
            f"Loaded metadata for {len(self.metadata)} samples with "
            f"{len(self.metadata.columns)} variables"
        )

        # Report missing data
        missing_samples = self.metadata.index[self.metadata.isnull().all(axis=1)]
        if len(missing_samples) > 0:
            logging.warning(f"Missing metadata for samples: {list(missing_samples)}")

    def validate_sample_size(self, groups):
        """Validate that groups have sufficient sample size for reliable PERMANOVA results."""
        for group_name, group_data in groups.items():
            if len(group_data) < 10:
                warnings.warn(f"Group '{group_name}' has only {len(group_data)} samples. "
                             f"PERMANOVA requires ≥10 samples per group for reliable results.")

    def analyze_variables(
        self, variables: Optional[List[str]] = None, n_permutations: int = 999
    ) -> pd.DataFrame:
        """Analyze all or specified metadata variables using PERMANOVA."""
        if self.metadata is None:
            raise ValueError("Metadata not loaded. Call load_metadata() first.")

        if variables is None:
            variables = list(self.metadata.columns)

        logging.info(f"Analyzing {len(variables)} metadata variables")

        permanova = PermanovaAnalyzer(self.distance_matrix, self.sample_names)
        results = []

        for variable in variables:
            logging.info(f"Analyzing variable: {variable}")

            if variable not in self.metadata.columns:
                logging.warning(f"Variable '{variable}' not found in metadata")
                continue

            # Prepare variable data
            var_data = self.metadata[variable].copy()

            # Validate sample size for categorical variables
            if var_data.dtype == "object":
                # Check group sizes before analysis
                groups = var_data.dropna().groupby(var_data.dropna()).apply(list).to_dict()
                self.validate_sample_size(groups)
                
                # Categorical variable
                var_data = var_data.astype("category")
                var_encoded = LabelEncoder().fit_transform(var_data.dropna())
                var_array = np.full(len(var_data), np.nan)
                var_array[~var_data.isna()] = var_encoded
            else:
                # Numerical variable - bin into categories for PERMANOVA
                var_array = var_data.values
                if not np.all(np.isnan(var_array)):
                    # Create quantile-based bins
                    valid_values = var_array[~np.isnan(var_array)]
                    if (
                        len(np.unique(valid_values)) > 10
                    ):  # Only bin if many unique values
                        quantiles = np.percentile(valid_values, [33, 67])
                        var_binned = np.full_like(var_array, np.nan)
                        var_binned[~np.isnan(var_array)] = np.digitize(
                            valid_values, quantiles
                        )
                        var_array = var_binned

            # Run PERMANOVA
            result = permanova.permanova_test(var_array, n_permutations)
            result["variable"] = variable
            result["variable_type"] = (
                "categorical"
                if self.metadata[variable].dtype == "object"
                else "numerical"
            )
            result["missing_count"] = self.metadata[variable].isna().sum()

            results.append(result)
            self.permanova_results[variable] = result

        # Create results DataFrame
        results_df = pd.DataFrame(results)
        if not results_df.empty:
            results_df = results_df.sort_values("r_squared", ascending=False)

        return results_df

    def identify_clusters(
        self,
        n_clusters_range: Tuple[int, int] = (2, 8),
        methods: List[str] = ["kmeans", "hierarchical"],
    ) -> Dict[str, Any]:
        """Identify natural clusters in the data using multiple methods."""
        logging.info("Identifying natural clusters in sample data")

        results = {}

        # Convert distance matrix to feature space using MDS
        from sklearn.manifold import MDS

        mds = MDS(
            n_components=min(10, len(self.sample_names) - 1),
            dissimilarity="precomputed",
            random_state=42,
        )
        X = mds.fit_transform(self.distance_matrix)

        for method in methods:
            logging.info(f"Clustering with {method}")
            method_results = {}

            for n_clusters in range(n_clusters_range[0], n_clusters_range[1] + 1):
                if n_clusters >= len(self.sample_names):
                    continue

                # Apply clustering method
                if method == "kmeans":
                    clusterer = KMeans(
                        n_clusters=n_clusters, random_state=42, n_init=10
                    )
                    labels = clusterer.fit_predict(X)
                elif method == "hierarchical":
                    clusterer = AgglomerativeClustering(
                        n_clusters=n_clusters, metric="precomputed", linkage="average"
                    )
                    labels = clusterer.fit_predict(self.distance_matrix)
                else:
                    continue

                # Calculate clustering metrics
                silhouette = silhouette_score(X, labels)

                method_results[n_clusters] = {
                    "labels": labels,
                    "silhouette_score": silhouette,
                    "n_clusters": n_clusters,
                }

            # Find optimal number of clusters
            if method_results:
                best_k = max(
                    method_results.keys(),
                    key=lambda k: method_results[k]["silhouette_score"],
                )
                method_results["optimal"] = method_results[best_k]

            results[method] = method_results

        self.cluster_results = results
        return results

    def compare_clustering_with_metadata(
        self, clustering_labels: np.ndarray, variables: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """Compare clustering results with metadata variables."""
        if self.metadata is None:
            raise ValueError("Metadata not loaded")

        if variables is None:
            variables = list(self.metadata.columns)

        results = []

        for variable in variables:
            if variable not in self.metadata.columns:
                continue

            var_data = self.metadata[variable].dropna()
            if len(var_data) == 0:
                continue

            # Get corresponding cluster labels
            valid_indices = [
                i for i, name in enumerate(self.sample_names) if name in var_data.index
            ]

            if len(valid_indices) < 2:
                continue

            cluster_subset = clustering_labels[valid_indices]

            # Calculate agreement metrics
            if var_data.dtype == "object":
                # Categorical variable
                var_encoded = LabelEncoder().fit_transform(var_data.values)
                ari = adjusted_rand_score(var_encoded, cluster_subset)
            else:
                # Numerical variable - use correlation with cluster centroids
                cluster_means = []
                for cluster_id in np.unique(cluster_subset):
                    cluster_mask = cluster_subset == cluster_id
                    if np.any(cluster_mask):
                        cluster_mean = var_data.iloc[cluster_mask].mean()
                        cluster_means.extend([cluster_mean] * np.sum(cluster_mask))

                if len(cluster_means) == len(var_data):
                    ari = stats.pearsonr(var_data.values, cluster_means)[0] ** 2
                else:
                    ari = np.nan

            results.append(
                {
                    "variable": variable,
                    "adjusted_rand_index": ari,
                    "variable_type": (
                        "categorical" if var_data.dtype == "object" else "numerical"
                    ),
                    "n_valid_samples": len(valid_indices),
                }
            )

        return pd.DataFrame(results).sort_values("adjusted_rand_index", ascending=False)


class MetadataVisualizer:
    """Generate visualizations for metadata analysis."""

    def __init__(self, sample_names: List[str], metadata: pd.DataFrame):
        self.sample_names = sample_names
        self.metadata = metadata

    def plot_variable_importance(self, results_df: pd.DataFrame, output_path: str):
        """Plot variable importance (R-squared values) from PERMANOVA."""
        plt.figure(figsize=(12, 8))

        # Filter out variables with NaN R-squared
        valid_results = results_df.dropna(subset=["r_squared"])

        if valid_results.empty:
            plt.text(
                0.5,
                0.5,
                "No valid results to display",
                ha="center",
                va="center",
                transform=plt.gca().transAxes,
            )
            plt.title("Variable Importance (R-squared from PERMANOVA)")
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            plt.close()
            return

        # Create color map based on significance
        colors = [
            "red" if p < 0.05 else "orange" if p < 0.1 else "gray"
            for p in valid_results["p_value"]
        ]

        # Horizontal bar plot
        y_pos = np.arange(len(valid_results))
        plt.barh(y_pos, valid_results["r_squared"], color=colors, alpha=0.7)

        # Customize plot
        plt.yticks(y_pos, valid_results["variable"])
        plt.xlabel("R-squared (Proportion of Variation Explained)")
        plt.title("Variable Importance (PERMANOVA Analysis)")
        plt.grid(axis="x", alpha=0.3)

        # Add significance legend
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(facecolor="red", alpha=0.7, label="p < 0.05"),
            Patch(facecolor="orange", alpha=0.7, label="0.05 ≤ p < 0.1"),
            Patch(facecolor="gray", alpha=0.7, label="p ≥ 0.1"),
        ]
        plt.legend(handles=legend_elements, loc="lower right")

        # Add R-squared values as text
        for i, (_, row) in enumerate(valid_results.iterrows()):
            plt.text(
                row["r_squared"] + 0.01,
                i,
                f"{row['r_squared']:.3f}",
                va="center",
                fontsize=9,
            )

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"Variable importance plot saved to {output_path}")

    def plot_samples_by_variable(
        self,
        pca_result: np.ndarray,
        variable: str,
        output_path: str,
        pca_obj: Optional[PCA] = None,
    ):
        """Plot PCA colored by a specific metadata variable."""
        plt.figure(figsize=(10, 8))

        if variable not in self.metadata.columns:
            logging.warning(f"Variable '{variable}' not found in metadata")
            return

        var_data = self.metadata[variable]

        # Handle missing values
        valid_mask = ~var_data.isna()

        if not np.any(valid_mask):
            plt.text(
                0.5,
                0.5,
                f"No valid data for variable: {variable}",
                ha="center",
                va="center",
                transform=plt.gca().transAxes,
            )
            plt.title(f"PCA colored by {variable}")
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            plt.close()
            return

        # Plot samples
        if var_data.dtype == "object":
            # Categorical variable
            unique_values = var_data.dropna().unique()
            colors = plt.cm.Set1(np.linspace(0, 1, len(unique_values)))

            for i, value in enumerate(unique_values):
                mask = (var_data == value) & valid_mask
                if np.any(mask):
                    plt.scatter(
                        pca_result[mask, 0],
                        pca_result[mask, 1],
                        c=[colors[i]],
                        label=str(value),
                        s=100,
                        alpha=0.7,
                    )

            plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        else:
            # Numerical variable
            valid_data = var_data[valid_mask]
            scatter = plt.scatter(
                pca_result[valid_mask, 0],
                pca_result[valid_mask, 1],
                c=valid_data,
                cmap="viridis",
                s=100,
                alpha=0.7,
            )
            plt.colorbar(scatter, label=variable)

        # Add sample labels
        for i, sample in enumerate(self.sample_names):
            if valid_mask.iloc[i]:
                plt.annotate(
                    sample,
                    (pca_result[i, 0], pca_result[i, 1]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                    alpha=0.7,
                )

        # Labels and title
        if pca_obj is not None:
            plt.xlabel(f"PC1 ({pca_obj.explained_variance_ratio_[0]:.1%} variance)")
            plt.ylabel(f"PC2 ({pca_obj.explained_variance_ratio_[1]:.1%} variance)")
        else:
            plt.xlabel("PC1")
            plt.ylabel("PC2")

        plt.title(f"PCA colored by {variable}")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"PCA plot for {variable} saved to {output_path}")

    def plot_clustering_results(
        self,
        pca_result: np.ndarray,
        cluster_labels: np.ndarray,
        method: str,
        n_clusters: int,
        output_path: str,
    ):
        """Plot clustering results on PCA space."""
        plt.figure(figsize=(10, 8))

        # Plot clusters
        unique_clusters = np.unique(cluster_labels)
        colors = plt.cm.Set1(np.linspace(0, 1, len(unique_clusters)))

        for i, cluster_id in enumerate(unique_clusters):
            mask = cluster_labels == cluster_id
            plt.scatter(
                pca_result[mask, 0],
                pca_result[mask, 1],
                c=[colors[i]],
                label=f"Cluster {cluster_id}",
                s=100,
                alpha=0.7,
            )

        # Add sample labels
        for i, sample in enumerate(self.sample_names):
            plt.annotate(
                sample,
                (pca_result[i, 0], pca_result[i, 1]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
                alpha=0.7,
            )

        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title(f"{method.title()} Clustering (k={n_clusters})")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"Clustering plot saved to {output_path}")


def generate_summary_report(
    results_df: pd.DataFrame, cluster_results: Dict, output_path: str
):
    """Generate a summary report of the metadata analysis."""

    with open(output_path, "w") as f:
        f.write("# MetaGrouper Phase 2: Metadata Analysis Report\n\n")

        # Variable importance section
        f.write("## Variable Importance (PERMANOVA Results)\n\n")
        f.write(
            "Variables ranked by proportion of variation explained (R-squared):\n\n"
        )

        if not results_df.empty:
            valid_results = results_df.dropna(subset=["r_squared"])

            for _, row in valid_results.iterrows():
                significance = (
                    "***"
                    if row["p_value"] < 0.001
                    else (
                        "**"
                        if row["p_value"] < 0.01
                        else (
                            "*"
                            if row["p_value"] < 0.05
                            else "." if row["p_value"] < 0.1 else ""
                        )
                    )
                )

                f.write(
                    f"- **{row['variable']}**: R² = {row['r_squared']:.3f}, "
                    f"p = {row['p_value']:.3f}{significance}\n"
                )
                f.write(f"  - Type: {row['variable_type']}\n")
                f.write(f"  - Valid samples: {row['n_samples']}\n")
                f.write(f"  - Groups: {row['n_groups']}\n\n")
        else:
            f.write("No valid results found.\n\n")

        # Clustering results section
        f.write("## Clustering Analysis\n\n")

        for method, method_results in cluster_results.items():
            f.write(f"### {method.title()} Clustering\n\n")

            if "optimal" in method_results:
                optimal = method_results["optimal"]
                f.write(f"- **Optimal clusters**: {optimal['n_clusters']}\n")
                f.write(
                    f"- **Silhouette score**: {optimal['silhouette_score']:.3f}\n\n"
                )

                # Show silhouette scores for different k values
                f.write("Silhouette scores for different k values:\n")
                for k in sorted(
                    [k for k in method_results.keys() if isinstance(k, int)]
                ):
                    score = method_results[k]["silhouette_score"]
                    f.write(f"- k={k}: {score:.3f}\n")
                f.write("\n")

        # Recommendations section
        f.write("## Recommendations\n\n")

        if not results_df.empty:
            valid_results = results_df.dropna(subset=["r_squared"])
            if not valid_results.empty:
                top_variable = valid_results.iloc[0]
                f.write(f"### Primary Grouping Variable\n\n")
                f.write(
                    f"**{top_variable['variable']}** explains the most variation "
                    f"({top_variable['r_squared']:.1%}) in sample similarities.\n\n"
                )

                if top_variable["p_value"] < 0.05:
                    f.write(
                        "This variable shows statistically significant association "
                        "with sample composition (p < 0.05).\n\n"
                    )
                    f.write(
                        "**Recommendation**: Consider grouping samples by this variable "
                        "for co-assembly.\n\n"
                    )
                else:
                    f.write(
                        "This association is not statistically significant (p ≥ 0.05).\n\n"
                    )
                    f.write(
                        "**Recommendation**: Consider individual sample assembly or "
                        "global co-assembly.\n\n"
                    )

        # Add interpretation guide
        f.write("## Interpretation Guide\n\n")
        f.write(
            "- **R-squared**: Proportion of variation in sample composition explained by the variable\n"
        )
        f.write(
            "- **p-value**: Statistical significance (< 0.05 is typically significant)\n"
        )
        f.write(
            "- **Silhouette score**: Quality of clustering (higher is better, > 0.5 is good)\n"
        )
        f.write("- Significance codes: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1\n")

    logging.info(f"Summary report saved to {output_path}")
