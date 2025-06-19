#!/usr/bin/env python3
"""
Phase 2: Metadata Analysis for MetaGrouper

Implements PERMANOVA testing, clustering analysis, and statistical reporting
for understanding metadata associations with k-mer similarity patterns.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import logging
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr
import warnings

# Try to import advanced statistics packages
try:
    from skbio import DistanceMatrix
    from skbio.stats.ordination import pcoa
    from skbio.stats.distance import permanova
    SKBIO_AVAILABLE = True
except ImportError:
    SKBIO_AVAILABLE = False
    logging.warning("scikit-bio not available - using simplified PERMANOVA")

try:
    from sklearn.manifold import MDS
    SKLEARN_MDS_AVAILABLE = True
except ImportError:
    SKLEARN_MDS_AVAILABLE = False


class MetadataAnalyzer:
    """Analyzes metadata associations with k-mer similarity patterns."""
    
    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str]):
        """
        Initialize metadata analyzer.
        
        Args:
            distance_matrix: Pairwise distance matrix between samples
            sample_names: List of sample names corresponding to matrix
        """
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names
        self.metadata = None
        self.analysis_results = {}
        
        # Validate inputs
        if distance_matrix.shape[0] != distance_matrix.shape[1]:
            raise ValueError("Distance matrix must be square")
        if len(sample_names) != distance_matrix.shape[0]:
            raise ValueError("Sample names must match distance matrix dimensions")
    
    def load_metadata(self, metadata_file: str, sample_id_column: str = 'sample_id'):
        """
        Load metadata from CSV/TSV file.
        
        Args:
            metadata_file: Path to metadata file
            sample_id_column: Column name containing sample IDs
        """
        try:
            # Detect separator
            with open(metadata_file, 'r') as f:
                first_line = f.readline()
                separator = '\t' if '\t' in first_line else ','
            
            # Load metadata
            self.metadata = pd.read_csv(metadata_file, sep=separator)
            
            # Validate sample ID column
            if sample_id_column not in self.metadata.columns:
                raise ValueError(f"Sample ID column '{sample_id_column}' not found in metadata")
            
            # Set sample ID as index
            self.metadata = self.metadata.set_index(sample_id_column)
            
            # Filter to samples present in distance matrix
            available_samples = set(self.sample_names)
            metadata_samples = set(self.metadata.index)
            
            common_samples = available_samples & metadata_samples
            missing_samples = available_samples - metadata_samples
            
            if missing_samples:
                logging.warning(f"Missing metadata for {len(missing_samples)} samples: {list(missing_samples)[:5]}...")
            
            # Reorder metadata to match distance matrix
            self.metadata = self.metadata.reindex(self.sample_names)
            
            logging.info(f"Loaded metadata for {len(common_samples)}/{len(self.sample_names)} samples")
            logging.info(f"Available variables: {list(self.metadata.columns)}")
            
        except Exception as e:
            raise RuntimeError(f"Failed to load metadata: {e}")
    
    def _simple_permanova(self, distances: np.ndarray, groups: np.ndarray, 
                         n_permutations: int = 999) -> Dict[str, float]:
        """
        Simple PERMANOVA implementation when scikit-bio is not available.
        
        Args:
            distances: Distance matrix
            groups: Group assignments
            n_permutations: Number of permutations
            
        Returns:
            Dictionary with test statistics
        """
        # Remove NaN samples
        valid_mask = ~pd.isna(groups)
        if not valid_mask.any():
            return {'f_statistic': np.nan, 'p_value': np.nan, 'r_squared': np.nan}
        
        distances = distances[valid_mask][:, valid_mask]
        groups = groups[valid_mask]
        
        # Convert to numeric if categorical
        if not np.issubdtype(groups.dtype, np.number):
            unique_groups = pd.unique(groups)
            group_map = {group: i for i, group in enumerate(unique_groups)}
            groups = np.array([group_map[g] for g in groups])
        
        n_samples = len(groups)
        if n_samples < 3:
            return {'f_statistic': np.nan, 'p_value': np.nan, 'r_squared': np.nan}
        
        # Calculate F-statistic
        distances_sq = distances ** 2
        
        # Total sum of squares
        grand_mean = np.mean(distances_sq)
        ss_total = np.sum((distances_sq - grand_mean) ** 2)
        
        # Within-group sum of squares
        unique_groups = np.unique(groups)
        ss_within = 0
        df_within = 0
        
        for group in unique_groups:
            group_mask = groups == group
            group_distances = distances_sq[group_mask][:, group_mask]
            if len(group_distances) > 1:
                group_mean = np.mean(group_distances)
                ss_within += np.sum((group_distances - group_mean) ** 2)
                df_within += len(group_distances) - 1
        
        # Between-group sum of squares
        ss_between = ss_total - ss_within
        df_between = len(unique_groups) - 1
        
        if df_within == 0 or df_between == 0:
            return {'f_statistic': np.nan, 'p_value': np.nan, 'r_squared': np.nan}
        
        # F-statistic
        ms_between = ss_between / df_between
        ms_within = ss_within / df_within
        f_stat = ms_between / ms_within if ms_within > 0 else np.nan
        
        # R-squared
        r_squared = ss_between / ss_total if ss_total > 0 else 0
        
        # Permutation test for p-value
        if np.isnan(f_stat):
            p_value = np.nan
        else:
            permuted_f_stats = []
            for _ in range(n_permutations):
                perm_groups = np.random.permutation(groups)
                
                # Calculate permuted F-statistic
                ss_within_perm = 0
                for group in unique_groups:
                    group_mask = perm_groups == group
                    group_distances = distances_sq[group_mask][:, group_mask]
                    if len(group_distances) > 1:
                        group_mean = np.mean(group_distances)
                        ss_within_perm += np.sum((group_distances - group_mean) ** 2)
                
                ss_between_perm = ss_total - ss_within_perm
                ms_between_perm = ss_between_perm / df_between
                ms_within_perm = ss_within_perm / df_within
                f_stat_perm = ms_between_perm / ms_within_perm if ms_within_perm > 0 else 0
                permuted_f_stats.append(f_stat_perm)
            
            # Calculate p-value
            permuted_f_stats = np.array(permuted_f_stats)
            p_value = np.mean(permuted_f_stats >= f_stat)
        
        return {
            'f_statistic': f_stat,
            'p_value': p_value,
            'r_squared': r_squared
        }
    
    def analyze_variables(self, variables: Optional[List[str]] = None, 
                         n_permutations: int = 999) -> pd.DataFrame:
        """
        Analyze metadata variables using PERMANOVA.
        
        Args:
            variables: List of variables to analyze (default: all)
            n_permutations: Number of permutations for significance testing
            
        Returns:
            DataFrame with analysis results
        """
        if self.metadata is None:
            raise ValueError("No metadata loaded. Call load_metadata() first.")
        
        if variables is None:
            variables = list(self.metadata.columns)
        
        results = []
        
        logging.info(f"Analyzing {len(variables)} variables with {n_permutations} permutations")
        
        for var in variables:
            if var not in self.metadata.columns:
                logging.warning(f"Variable '{var}' not found in metadata")
                continue
            
            try:
                groups = self.metadata[var].values
                
                # Skip if all values are NaN or all the same
                valid_groups = groups[~pd.isna(groups)]
                if len(valid_groups) == 0 or len(np.unique(valid_groups)) < 2:
                    logging.warning(f"Variable '{var}' has insufficient variation")
                    continue
                
                # Run PERMANOVA
                if SKBIO_AVAILABLE and len(valid_groups) > 10:
                    try:
                        # Use scikit-bio for more robust analysis
                        valid_mask = ~pd.isna(groups)
                        valid_distances = self.distance_matrix[valid_mask][:, valid_mask]
                        valid_groups = groups[valid_mask]
                        
                        dm = DistanceMatrix(valid_distances)
                        result = permanova(dm, valid_groups, permutations=n_permutations)
                        
                        results.append({
                            'variable': var,
                            'f_statistic': result['test statistic'],
                            'p_value': result['p-value'],
                            'r_squared': result['test statistic'] / (result['test statistic'] + len(valid_groups) - len(np.unique(valid_groups))),
                            'n_samples': len(valid_groups),
                            'n_groups': len(np.unique(valid_groups))
                        })
                    except Exception as e:
                        logging.warning(f"scikit-bio PERMANOVA failed for {var}: {e}")
                        # Fall back to simple implementation
                        result = self._simple_permanova(self.distance_matrix, groups, n_permutations)
                        results.append({
                            'variable': var,
                            'f_statistic': result['f_statistic'],
                            'p_value': result['p_value'],
                            'r_squared': result['r_squared'],
                            'n_samples': len(valid_groups),
                            'n_groups': len(np.unique(valid_groups))
                        })
                else:
                    # Use simple implementation
                    result = self._simple_permanova(self.distance_matrix, groups, n_permutations)
                    results.append({
                        'variable': var,
                        'f_statistic': result['f_statistic'],
                        'p_value': result['p_value'],
                        'r_squared': result['r_squared'],
                        'n_samples': len(valid_groups),
                        'n_groups': len(np.unique(valid_groups))
                    })
                
            except Exception as e:
                logging.error(f"Failed to analyze variable '{var}': {e}")
                continue
        
        if not results:
            logging.warning("No variables could be analyzed")
            return pd.DataFrame()
        
        results_df = pd.DataFrame(results)
        
        # Sort by R-squared (effect size)  
        results_df = results_df.sort_values('r_squared', ascending=False)
        
        # Store results
        self.analysis_results['permanova'] = results_df
        
        logging.info(f"Successfully analyzed {len(results_df)} variables")
        
        return results_df
    
    def identify_clusters(self, n_clusters_range: Tuple[int, int] = (2, 8)) -> Dict[str, Any]:
        """
        Identify optimal clustering of samples.
        
        Args:
            n_clusters_range: Range of cluster numbers to test
            
        Returns:
            Dictionary with clustering results
        """
        min_clusters, max_clusters = n_clusters_range
        max_clusters = min(max_clusters, len(self.sample_names) - 1)
        
        if max_clusters <= min_clusters:
            logging.warning("Not enough samples for clustering analysis")
            return {}
        
        results = {}
        
        # K-means clustering
        kmeans_results = {'silhouette_scores': {}, 'labels': {}}
        best_kmeans = {'n_clusters': None, 'score': -1, 'labels': None}
        
        for n_clusters in range(min_clusters, max_clusters + 1):
            try:
                # Convert distance to similarity for k-means
                similarity_matrix = 1 - self.distance_matrix
                
                kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
                labels = kmeans.fit_predict(similarity_matrix)
                
                # Calculate silhouette score
                silhouette = silhouette_score(similarity_matrix, labels)
                
                kmeans_results['silhouette_scores'][n_clusters] = silhouette
                kmeans_results['labels'][n_clusters] = labels
                
                if silhouette > best_kmeans['score']:
                    best_kmeans = {
                        'n_clusters': n_clusters,
                        'score': silhouette,
                        'labels': labels
                    }
                    
            except Exception as e:
                logging.warning(f"K-means clustering failed for {n_clusters} clusters: {e}")
        
        if best_kmeans['n_clusters']:
            kmeans_results['optimal'] = best_kmeans
        
        results['kmeans'] = kmeans_results
        
        # Hierarchical clustering
        hierarchical_results = {'silhouette_scores': {}, 'labels': {}}
        best_hierarchical = {'n_clusters': None, 'score': -1, 'labels': None}
        
        try:
            # Create linkage matrix
            condensed_distances = squareform(self.distance_matrix, checks=False)
            linkage_matrix = linkage(condensed_distances, method='ward')
            
            for n_clusters in range(min_clusters, max_clusters + 1):
                try:
                    labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust') - 1
                    
                    # Calculate silhouette score
                    silhouette = silhouette_score(self.distance_matrix, labels, metric='precomputed')
                    
                    hierarchical_results['silhouette_scores'][n_clusters] = silhouette
                    hierarchical_results['labels'][n_clusters] = labels
                    
                    if silhouette > best_hierarchical['score']:
                        best_hierarchical = {
                            'n_clusters': n_clusters,
                            'score': silhouette,
                            'labels': labels
                        }
                        
                except Exception as e:
                    logging.warning(f"Hierarchical clustering failed for {n_clusters} clusters: {e}")
            
            if best_hierarchical['n_clusters']:
                hierarchical_results['optimal'] = best_hierarchical
                hierarchical_results['linkage_matrix'] = linkage_matrix
            
        except Exception as e:
            logging.warning(f"Hierarchical clustering failed: {e}")
        
        results['hierarchical'] = hierarchical_results
        
        # Store results
        self.analysis_results['clustering'] = results
        
        logging.info(f"Clustering analysis complete. Best k-means: {best_kmeans['n_clusters']} clusters "
                    f"(silhouette: {best_kmeans['score']:.3f})")
        
        return results


class MetadataVisualizer:
    """Creates visualizations for metadata analysis results."""
    
    def __init__(self, sample_names: List[str], metadata: pd.DataFrame):
        """
        Initialize visualizer.
        
        Args:
            sample_names: List of sample names
            metadata: Metadata DataFrame
        """
        self.sample_names = sample_names
        self.metadata = metadata
    
    def plot_variable_importance(self, results_df: pd.DataFrame, output_path: str):
        """Plot variable importance from PERMANOVA results."""
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            if results_df.empty:
                logging.warning("No results to plot")
                return
            
            # Filter out NaN values
            plot_data = results_df.dropna(subset=['r_squared'])
            
            if plot_data.empty:
                logging.warning("No valid results to plot")
                return
            
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # R-squared plot
            bars1 = ax1.barh(plot_data['variable'], plot_data['r_squared'])
            ax1.set_xlabel('R² (Effect Size)')
            ax1.set_title('Variable Importance (R²)')
            ax1.set_xlim(0, max(0.5, plot_data['r_squared'].max() * 1.1))
            
            # Color bars by significance
            for i, (_, row) in enumerate(plot_data.iterrows()):
                color = 'red' if row['p_value'] < 0.05 else 'lightblue'
                bars1[i].set_color(color)
            
            # P-value plot
            neg_log_p = -np.log10(plot_data['p_value'].clip(lower=1e-10))
            bars2 = ax2.barh(plot_data['variable'], neg_log_p)
            ax2.set_xlabel('-log₁₀(p-value)')
            ax2.set_title('Statistical Significance')
            ax2.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')
            ax2.legend()
            
            # Color bars by significance
            for i, p_val in enumerate(plot_data['p_value']):
                color = 'red' if p_val < 0.05 else 'lightblue'
                bars2[i].set_color(color)
            
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            logging.info(f"Variable importance plot saved to {output_path}")
            
        except ImportError:
            logging.warning("Matplotlib not available - skipping visualization")
        except Exception as e:
            logging.error(f"Failed to create variable importance plot: {e}")
    
    def plot_samples_by_variable(self, pca_result: np.ndarray, variable: str, 
                               output_path: str, pca_obj=None):
        """Plot PCA colored by metadata variable."""
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            if variable not in self.metadata.columns:
                logging.warning(f"Variable {variable} not found in metadata")
                return
            
            # Get variable values
            var_values = self.metadata[variable].values
            
            # Handle different variable types
            if pd.api.types.is_numeric_dtype(var_values):
                # Continuous variable
                scatter = plt.scatter(pca_result[:, 0], pca_result[:, 1], 
                                    c=var_values, cmap='viridis', alpha=0.7)
                plt.colorbar(scatter, label=variable)
            else:
                # Categorical variable
                unique_values = pd.unique(var_values[~pd.isna(var_values)])
                colors = plt.cm.Set1(np.linspace(0, 1, len(unique_values)))
                
                for i, value in enumerate(unique_values):
                    mask = var_values == value
                    plt.scatter(pca_result[mask, 0], pca_result[mask, 1], 
                              c=[colors[i]], label=str(value), alpha=0.7)
                
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Labels and formatting
            if pca_obj and hasattr(pca_obj, 'explained_variance_ratio_'):
                plt.xlabel(f'PC1 ({pca_obj.explained_variance_ratio_[0]:.1%} variance)')
                plt.ylabel(f'PC2 ({pca_obj.explained_variance_ratio_[1]:.1%} variance)')
            else:
                plt.xlabel('PC1')
                plt.ylabel('PC2')
            
            plt.title(f'PCA colored by {variable}')
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            logging.info(f"PCA plot by {variable} saved to {output_path}")
            
        except ImportError:
            logging.warning("Matplotlib not available - skipping visualization")
        except Exception as e:
            logging.error(f"Failed to create PCA plot for {variable}: {e}")
    
    def plot_clustering_results(self, pca_result: np.ndarray, cluster_labels: np.ndarray,
                              method: str, n_clusters: int, output_path: str):
        """Plot clustering results on PCA."""
        try:
            import matplotlib.pyplot as plt
            
            colors = plt.cm.Set1(np.linspace(0, 1, n_clusters))
            
            for cluster_id in range(n_clusters):
                mask = cluster_labels == cluster_id
                plt.scatter(pca_result[mask, 0], pca_result[mask, 1],
                          c=[colors[cluster_id]], label=f'Cluster {cluster_id + 1}',
                          alpha=0.7)
            
            plt.xlabel('PC1')
            plt.ylabel('PC2')
            plt.title(f'{method.title()} Clustering ({n_clusters} clusters)')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            logging.info(f"Clustering plot saved to {output_path}")
            
        except ImportError:
            logging.warning("Matplotlib not available - skipping visualization")
        except Exception as e:
            logging.error(f"Failed to create clustering plot: {e}")


def generate_summary_report(metadata_results: pd.DataFrame, cluster_results: Dict[str, Any],
                          output_path: str):
    """Generate comprehensive summary report."""
    try:
        report_lines = [
            "# MetaGrouper Phase 2: Metadata Analysis Report",
            "",
            "## Variable Analysis Results",
            ""
        ]
        
        if not metadata_results.empty:
            # Sort by R-squared
            sorted_results = metadata_results.sort_values('r_squared', ascending=False)
            
            report_lines.extend([
                "### Top Variables by Effect Size (R²)",
                ""
            ])
            
            for _, row in sorted_results.head(10).iterrows():
                significance = "**" if row['p_value'] < 0.05 else ""
                report_lines.append(
                    f"- {significance}{row['variable']}{significance}: "
                    f"R² = {row['r_squared']:.3f}, p = {row['p_value']:.3f} "
                    f"({row['n_groups']} groups, {row['n_samples']} samples)"
                )
            
            report_lines.extend([
                "",
                "### Statistical Significance",
                ""
            ])
            
            significant_vars = sorted_results[sorted_results['p_value'] < 0.05]
            if len(significant_vars) > 0:
                report_lines.append(f"**{len(significant_vars)} variables** show significant associations (p < 0.05):")
                for _, row in significant_vars.iterrows():
                    report_lines.append(f"- {row['variable']}: R² = {row['r_squared']:.3f}, p = {row['p_value']:.3f}")
            else:
                report_lines.append("**No variables** show significant associations (p < 0.05).")
        else:
            report_lines.append("No metadata variables could be analyzed.")
        
        # Clustering results
        if cluster_results:
            report_lines.extend([
                "",
                "## Clustering Analysis",
                ""
            ])
            
            for method, results in cluster_results.items():
                if 'optimal' in results:
                    optimal = results['optimal']
                    report_lines.extend([
                        f"### {method.title()} Clustering",
                        f"- **Optimal clusters:** {optimal['n_clusters']}",
                        f"- **Silhouette score:** {optimal['score']:.3f}",
                        ""
                    ])
        
        # Interpretation guidelines
        report_lines.extend([
            "",
            "## Interpretation Guidelines",
            "",
            "### Effect Size (R²)",
            "- **R² > 0.20:** Strong association (>20% variance explained)",
            "- **R² = 0.10-0.20:** Moderate association",
            "- **R² < 0.10:** Weak association",
            "",
            "### Statistical Significance",
            "- **p < 0.05:** Statistically significant",
            "- **p < 0.01:** Highly significant",
            "- **p < 0.001:** Very highly significant",
            "",
            "### Clustering Quality",
            "- **Silhouette > 0.70:** Strong clustering",
            "- **Silhouette = 0.50-0.70:** Reasonable clustering", 
            "- **Silhouette < 0.50:** Weak clustering",
            ""
        ])
        
        # Write report
        with open(output_path, 'w') as f:
            f.write('\n'.join(report_lines))
        
        logging.info(f"Summary report saved to {output_path}")
        
    except Exception as e:
        logging.error(f"Failed to generate summary report: {e}")