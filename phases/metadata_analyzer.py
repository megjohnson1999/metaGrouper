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
import re
import os

warnings.filterwarnings("ignore")


def detect_file_paths(values: pd.Series, threshold: float = 0.8) -> bool:
    """
    Detect if a column contains file paths.
    
    Args:
        values: Series of column values (non-null)
        threshold: Minimum fraction of values that should look like paths
        
    Returns:
        True if column appears to contain file paths
    """
    if len(values) == 0:
        return False
    
    path_indicators = [
        r'/',  # Unix paths
        r'\\',  # Windows paths  
        r'\.[a-zA-Z0-9]{2,4}$',  # File extensions
        r'^[A-Z]:\\',  # Windows drive letters
        r'/[a-zA-Z0-9_.-]+/',  # Directory structure
    ]
    
    path_count = 0
    for value in values.astype(str):
        if any(re.search(pattern, value) for pattern in path_indicators):
            path_count += 1
    
    return (path_count / len(values)) >= threshold


def detect_sequential_ids(values: pd.Series, threshold: float = 0.8) -> bool:
    """
    Detect if a column contains sequential IDs or row numbers.
    
    Args:
        values: Series of column values (non-null)
        threshold: Minimum fraction that should be sequential
        
    Returns:
        True if column appears to be sequential IDs
    """
    if len(values) < 3:
        return False
    
    try:
        # Try to convert to numeric
        numeric_values = pd.to_numeric(values, errors='coerce')
        if numeric_values.isna().sum() > len(values) * 0.1:  # Too many non-numeric
            return False
        
        # Check if all values are unique integers (like database row IDs)
        unique_values = numeric_values.dropna().unique()
        if len(unique_values) == len(values) and all(v == int(v) for v in unique_values):
            # Check if the range is suspiciously sequential (like 1,2,3...n)
            min_val, max_val = min(unique_values), max(unique_values)
            expected_range = max_val - min_val + 1
            if expected_range == len(unique_values) and min_val <= 10:
                return True  # Likely database row IDs starting from low numbers
        
        # Check if values are mostly sequential
        sorted_values = sorted(numeric_values.dropna())
        sequential_count = 0
        
        for i in range(1, len(sorted_values)):
            if abs(sorted_values[i] - sorted_values[i-1]) <= 1:
                sequential_count += 1
        
        return (sequential_count / (len(sorted_values) - 1)) >= threshold
        
    except:
        return False


def detect_technical_hashes(values: pd.Series, threshold: float = 0.7) -> bool:
    """
    Detect if a column contains UUIDs, hashes, or other technical identifiers.
    
    Args:
        values: Series of column values (non-null)
        threshold: Minimum fraction that should look like technical IDs
        
    Returns:
        True if column appears to contain technical identifiers
    """
    if len(values) == 0:
        return False
    
    technical_patterns = [
        r'^[A-Fa-f0-9]{8}-[A-Fa-f0-9]{4}-[A-Fa-f0-9]{4}-[A-Fa-f0-9]{4}-[A-Fa-f0-9]{12}$',  # UUID
        r'^[A-Fa-f0-9]{32}$',  # MD5
        r'^[A-Fa-f0-9]{40}$',  # SHA1
        r'^[A-Za-z0-9_-]{20,}$',  # Long alphanumeric strings (20+ chars)
        r'^\w+_I\d+_\d+_',  # Sequencing IDs like M595_I7606_7544_...
        r'^\w+\.\w+\.\w+\.\w+',  # Multi-dot separated technical IDs
        r'^[A-Z]\d+_[A-Z]\d+_',  # Pattern like M595_I7606_
        r'\.fastq',  # FASTQ file names
        r'\.gz$',  # Compressed files
    ]
    
    technical_count = 0
    for value in values.astype(str):
        if any(re.search(pattern, value) for pattern in technical_patterns):
            technical_count += 1
    
    return (technical_count / len(values)) >= threshold


def calculate_information_content(values: pd.Series) -> float:
    """
    Calculate the information content (entropy) of a column.
    
    Args:
        values: Series of column values (non-null)
        
    Returns:
        Entropy value (higher = more informative)
    """
    if len(values) == 0:
        return 0.0
    
    # Calculate value frequencies
    value_counts = values.value_counts(normalize=True)
    
    # Calculate entropy
    entropy = -sum(p * np.log2(p) for p in value_counts if p > 0)
    
    return entropy


def detect_constant_values(values: pd.Series, max_unique_ratio: float = 0.005) -> bool:
    """
    Detect if a column has essentially constant values.
    Only excludes if there's truly no meaningful variation.
    
    Args:
        values: Series of column values (non-null)
        max_unique_ratio: Maximum ratio of unique values to consider constant
        
    Returns:
        True if column is essentially constant
    """
    if len(values) == 0:
        return True
    
    unique_count = len(values.unique())
    
    # Never exclude binary variables (2 unique values) - often biologically important
    if unique_count == 2:
        return False
    
    # Don't exclude if there are 3+ unique values
    if unique_count >= 3:
        return False
    
    # Only exclude single-value columns
    return unique_count <= 1


def is_biological_id(column_name: str) -> bool:
    """
    Check if an ID column represents biological grouping rather than technical identifiers.
    
    Args:
        column_name: Name of the metadata column
        
    Returns:
        True if the column represents biological grouping, False if technical
    """
    biological_patterns = [
        'PID', 'Patient', 'Subject', 'Participant', 'Individual',
        'GEMM', 'Study', 'Cohort', 'Group', 'Batch', 'Site',
        'Family', 'Twin', 'Sibling', 'Parent', 'Child'
    ]
    column_lower = column_name.lower()
    return any(pattern.lower() in column_lower for pattern in biological_patterns)


def filter_metadata_variables(
    metadata: pd.DataFrame, 
    auto_filter: bool = True,
    exclude_variables: Optional[List[str]] = None,
    include_variables: Optional[List[str]] = None,
    max_unique_ratio: float = 0.2,  # Reduced from 0.5 to 0.2
    min_unique_count: int = 2,
    max_missing_ratio: float = 0.2
) -> Tuple[List[str], Dict[str, str]]:
    """
    Intelligently filter metadata variables to focus on biologically relevant columns.
    
    Args:
        metadata: DataFrame with metadata
        auto_filter: Whether to apply automatic filtering
        exclude_variables: List of variables to explicitly exclude
        include_variables: List of variables to explicitly include (overrides auto-filtering)
        max_unique_ratio: Maximum ratio of unique values to total samples (0.5 = 50%)
        min_unique_count: Minimum number of unique values required
        max_missing_ratio: Maximum ratio of missing values allowed (0.2 = 20%)
        
    Returns:
        Tuple of (filtered_variables, exclusion_reasons)
    """
    exclude_variables = exclude_variables or []
    include_variables = include_variables or []
    exclusion_reasons = {}
    
    # If specific variables are requested, use only those
    if include_variables:
        valid_vars = [var for var in include_variables if var in metadata.columns]
        excluded_vars = [var for var in include_variables if var not in metadata.columns]
        for var in excluded_vars:
            exclusion_reasons[var] = "Variable not found in metadata"
        return valid_vars, exclusion_reasons
    
    if not auto_filter:
        # Return all columns except explicitly excluded ones
        all_vars = [col for col in metadata.columns if col not in exclude_variables]
        for var in exclude_variables:
            if var in metadata.columns:
                exclusion_reasons[var] = "Manually excluded"
        return all_vars, exclusion_reasons
    
    # Auto-filtering logic
    filtered_variables = []
    
    for column_name in metadata.columns:
        # Skip if manually excluded
        if column_name in exclude_variables:
            exclusion_reasons[column_name] = "Manually excluded"
            continue
            
        # Get column values
        values = metadata[column_name].dropna()
        
        # Content-based detection of technical columns
        if detect_file_paths(values):
            exclusion_reasons[column_name] = "Contains file paths"
            continue
            
        if detect_sequential_ids(values):
            exclusion_reasons[column_name] = "Sequential row numbers/IDs"
            continue
            
        if detect_technical_hashes(values):
            exclusion_reasons[column_name] = "Technical hashes/UUIDs"
            continue
            
        if detect_constant_values(values):
            exclusion_reasons[column_name] = "Essentially constant values"
            continue
        
        # Statistical filters
        if len(values) == 0:
            exclusion_reasons[column_name] = "No valid values"
            continue
            
        unique_count = len(values.unique())
        total_count = len(metadata)
        
        # Too few unique values (no variation)
        if unique_count < min_unique_count:
            exclusion_reasons[column_name] = f"Too few unique values ({unique_count})"
            continue
            
        # Too many unique values (likely continuous ID or noise)
        unique_ratio = unique_count / total_count
        if unique_ratio > max_unique_ratio:
            # Check if this is an individual-level identifier (patient IDs, etc.)
            if unique_ratio > 0.8:  # >80% unique values
                exclusion_reasons[column_name] = f"Individual-level identifier ({unique_ratio:.1%} unique values). Not suitable for grouping analysis"
                continue
            
            # Check if high uniqueness is informative vs noise
            information_content = calculate_information_content(values)
            
            # Allow high uniqueness if it has high information content
            # and doesn't look like technical noise
            if information_content < 2.0:  # Low information content
                exclusion_reasons[column_name] = f"Too many unique values with low information content ({unique_ratio:.1%})"
                continue
                
        # Too many missing values
        missing_ratio = (total_count - len(values)) / total_count
        if missing_ratio > max_missing_ratio:
            exclusion_reasons[column_name] = f"Too many missing values ({missing_ratio:.1%})"
            continue
            
        # If we get here, include the variable
        filtered_variables.append(column_name)
    
    # Prioritize biological variables by sorting
    def biological_priority(var_name):
        """Sort biological variables first."""
        biological_keywords = [
            'disease', 'diagnosis', 'status', 'group', 'case', 'control',
            'sex', 'gender', 'age', 'delivery', 'birth', 'hla', 'genetic',
            'country', 'location', 'site', 'time', 'month', 'year', 'onset'
        ]
        var_lower = var_name.lower()
        bio_score = sum(1 for keyword in biological_keywords if keyword in var_lower)
        return (-bio_score, var_name)  # Negative for descending order
    
    filtered_variables.sort(key=biological_priority)
    
    logging.info(f"Auto-filtered metadata: {len(filtered_variables)}/{len(metadata.columns)} variables retained")
    
    return filtered_variables, exclusion_reasons


class PermanovaAnalyzer:
    """PERMANOVA (Permutational Multivariate Analysis of Variance) implementation."""

    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str]):
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names
        self.n_samples = len(sample_names)
        
        # Validate distance matrix
        self._validate_distance_matrix()
    
    def _validate_distance_matrix(self):
        """
        Validate that the distance matrix meets PERMANOVA requirements.
        
        Checks:
        - Square matrix
        - Symmetric
        - Non-negative values
        - Zero diagonal
        - Finite values (no NaN or inf)
        """
        # Check if square
        if self.distance_matrix.shape[0] != self.distance_matrix.shape[1]:
            raise ValueError(f"Distance matrix must be square. Got shape {self.distance_matrix.shape}")
        
        # Check size matches sample names
        if self.distance_matrix.shape[0] != len(self.sample_names):
            raise ValueError(f"Distance matrix size ({self.distance_matrix.shape[0]}) doesn't match "
                           f"number of samples ({len(self.sample_names)})")
        
        # Check for finite values
        if not np.all(np.isfinite(self.distance_matrix)):
            raise ValueError("Distance matrix contains NaN or infinite values")
        
        # Check symmetry (with tolerance for floating point errors)
        if not np.allclose(self.distance_matrix, self.distance_matrix.T, rtol=1e-10):
            # Make it symmetric by averaging with transpose
            logging.warning("Distance matrix was not perfectly symmetric. Forcing symmetry.")
            self.distance_matrix = (self.distance_matrix + self.distance_matrix.T) / 2
        
        # Check diagonal is zero
        if not np.allclose(np.diag(self.distance_matrix), 0, atol=1e-10):
            logging.warning("Distance matrix diagonal contains non-zero values. Setting to zero.")
            np.fill_diagonal(self.distance_matrix, 0)
        
        # Check non-negative
        if np.any(self.distance_matrix < 0):
            raise ValueError("Distance matrix contains negative values")
        
        # Check for variation
        if np.allclose(self.distance_matrix, 0):
            raise ValueError("Distance matrix contains all zeros (no variation)")

    def _calculate_sum_of_squares(
        self, distance_matrix: np.ndarray, groups: np.ndarray
    ) -> Tuple[float, float]:
        """
        Calculate within-group and total sum of squares for PERMANOVA.
        
        Uses the correct formulation from Anderson (2001):
        - Total SS = (1/n) * sum of squared distances / 2
        - Within SS = sum over groups of (1/n_g) * sum of squared distances within group / 2
        """
        n = distance_matrix.shape[0]

        # Total sum of squares
        # SS_T = (1/n) * sum(d_ij^2) / 2
        total_ss = np.sum(distance_matrix**2) / (2 * n)

        # Within-group sum of squares
        within_ss = 0
        unique_groups = np.unique(groups)

        for group in unique_groups:
            group_indices = np.where(groups == group)[0]
            n_group = len(group_indices)
            
            if n_group > 1:
                # Extract submatrix for this group
                group_distances = distance_matrix[np.ix_(group_indices, group_indices)]
                # SS_W,g = (1/n_g) * sum(d_ij^2 within group) / 2
                within_ss += np.sum(group_distances**2) / (2 * n_group)

        return within_ss, total_ss
    
    def permdisp_test(self, distance_matrix: np.ndarray, groups: np.ndarray, n_permutations: int = 999) -> Dict[str, Any]:
        """
        PERMDISP test for homogeneity of multivariate dispersions.
        
        Tests the null hypothesis that groups have the same dispersion (variance).
        This is an important assumption check for PERMANOVA.
        
        Based on Anderson (2006) "Distance-based tests for homogeneity of multivariate dispersions"
        """
        n = distance_matrix.shape[0]
        unique_groups = np.unique(groups)
        n_groups = len(unique_groups)
        
        if n_groups < 2:
            return {"f_statistic": np.nan, "p_value": np.nan, "warning": None}
        
        # Calculate spatial median (centroid) for each group using PCoA
        # First, convert distance matrix to coordinates using classical MDS
        from sklearn.manifold import MDS
        
        # Ensure symmetric matrix
        dist_symmetric = (distance_matrix + distance_matrix.T) / 2
        np.fill_diagonal(dist_symmetric, 0)
        
        # Use MDS to get coordinates (using as many dimensions as needed)
        n_components = min(n - 1, 10)  # Use up to 10 dimensions
        mds = MDS(n_components=n_components, dissimilarity='precomputed', random_state=42)
        coords = mds.fit_transform(dist_symmetric)
        
        # Calculate centroid for each group and distances to centroid
        group_dispersions = {}
        all_distances_to_centroid = np.zeros(n)
        
        for group in unique_groups:
            group_mask = groups == group
            group_coords = coords[group_mask]
            
            # Calculate spatial median (more robust than mean)
            group_centroid = np.median(group_coords, axis=0)
            
            # Calculate distances from each point to group centroid
            distances = np.sqrt(np.sum((group_coords - group_centroid)**2, axis=1))
            group_dispersions[group] = distances
            all_distances_to_centroid[group_mask] = distances
        
        # Perform ANOVA on distances to centroid
        # Calculate within-group and between-group sum of squares
        grand_mean = np.mean(all_distances_to_centroid)
        
        ss_total = np.sum((all_distances_to_centroid - grand_mean)**2)
        ss_within = 0
        
        for group in unique_groups:
            group_distances = group_dispersions[group]
            group_mean = np.mean(group_distances)
            ss_within += np.sum((group_distances - group_mean)**2)
        
        ss_between = ss_total - ss_within
        
        # Degrees of freedom
        df_between = n_groups - 1
        df_within = n - n_groups
        
        if df_within <= 0:
            return {"f_statistic": np.nan, "p_value": np.nan, "warning": "Insufficient degrees of freedom"}
        
        # F-statistic
        ms_between = ss_between / df_between
        ms_within = ss_within / df_within
        f_observed = ms_between / ms_within if ms_within > 0 else np.nan
        
        # Permutation test
        f_permuted = []
        for _ in range(n_permutations):
            # Permute group labels
            perm_groups = np.random.permutation(groups)
            
            # Recalculate F-statistic with permuted groups
            perm_distances = np.zeros(n)
            for group in unique_groups:
                group_mask = perm_groups == group
                if np.sum(group_mask) > 0:
                    group_coords = coords[group_mask]
                    group_centroid = np.median(group_coords, axis=0)
                    distances = np.sqrt(np.sum((group_coords - group_centroid)**2, axis=1))
                    perm_distances[group_mask] = distances
            
            # Calculate F for permuted data
            perm_mean = np.mean(perm_distances)
            perm_ss_total = np.sum((perm_distances - perm_mean)**2)
            perm_ss_within = 0
            
            for group in unique_groups:
                group_mask = perm_groups == group
                if np.sum(group_mask) > 0:
                    group_distances = perm_distances[group_mask]
                    group_mean = np.mean(group_distances)
                    perm_ss_within += np.sum((group_distances - group_mean)**2)
            
            perm_ss_between = perm_ss_total - perm_ss_within
            perm_ms_between = perm_ss_between / df_between
            perm_ms_within = perm_ss_within / df_within
            
            if perm_ms_within > 0:
                f_perm = perm_ms_between / perm_ms_within
                f_permuted.append(f_perm)
        
        # Calculate p-value
        if f_permuted and not np.isnan(f_observed):
            p_value = (np.sum(np.array(f_permuted) >= f_observed) + 1) / (len(f_permuted) + 1)
        else:
            p_value = np.nan
        
        # Determine warning
        warning = None
        if p_value < 0.05:
            warning = "Groups have significantly different dispersions (heterogeneous variances). PERMANOVA results may be confounded."
        
        return {
            "f_statistic": f_observed,
            "p_value": p_value,
            "warning": warning
        }

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
        
        # Adjusted R-squared (penalizes for degrees of freedom)
        # Formula: R²_adj = 1 - (1-R²) * (n-1)/(n-k-1)
        # where n = number of samples, k = number of groups - 1
        if n_samples > df_between + 1:
            r_squared_adj = 1 - (1 - r_squared) * (n_samples - 1) / (n_samples - df_between - 1)
        else:
            r_squared_adj = r_squared  # Fallback when adjustment impossible
        
        # Omega squared (ω²) - less biased effect size for ANOVA-type designs
        # Formula: ω² = (SS_between - df_between * MS_within) / (SS_total + MS_within)
        # This is more conservative than η² (which is equivalent to R²)
        ms_within = within_ss / df_within if df_within > 0 else 0
        omega_squared = max(0, (between_ss - df_between * ms_within) / (total_ss + ms_within))
        
        # Partial omega squared for comparison
        # Formula: ω²_p = (SS_between - df_between * MS_within) / (SS_between + (n - df_between) * MS_within)
        partial_omega_squared = max(0, (between_ss - df_between * ms_within) / 
                                   (between_ss + (n_samples - df_between) * ms_within)) if ms_within > 0 else 0

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

        # Effect size categorization based on omega squared (more conservative)
        # Using Cohen's f guidelines converted to ω²: small=0.01, medium=0.06, large=0.14
        if omega_squared < 0.01:
            effect_size = "negligible"
        elif omega_squared < 0.06:
            effect_size = "small"
        elif omega_squared < 0.14:
            effect_size = "medium"
        else:
            effect_size = "large"
        
        # Power analysis warning
        power_warning = None
        if n_groups > 2:
            expected_samples_per_group = n_samples / n_groups
            if expected_samples_per_group < 10:
                power_warning = f"Low statistical power: {expected_samples_per_group:.1f} samples per group (recommended: ≥10)"
        
        # Run PERMDISP test for homogeneity of dispersions
        permdisp_result = self.permdisp_test(valid_distance_matrix, valid_groups, n_permutations=min(99, n_permutations))
        
        return {
            "f_statistic": f_observed,
            "p_value": p_value,
            "r_squared": r_squared,
            "r_squared_adj": r_squared_adj,
            "omega_squared": omega_squared,
            "partial_omega_squared": partial_omega_squared,
            "effect_size": effect_size,
            "power_warning": power_warning,
            "n_samples": n_samples,
            "n_groups": n_groups,
            "permdisp_f": permdisp_result["f_statistic"],
            "permdisp_p": permdisp_result["p_value"],
            "variance_warning": permdisp_result["warning"]
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
        
        # Check for duplicate sample IDs and handle them
        if self.metadata.index.duplicated().any():
            duplicates = self.metadata.index[self.metadata.index.duplicated(keep=False)]
            logging.warning(f"Found {len(duplicates)} duplicate sample IDs in metadata: {list(duplicates.unique())}")
            
            # Remove duplicates, keeping the first occurrence
            self.metadata = self.metadata[~self.metadata.index.duplicated(keep='first')]
            logging.info(f"Removed duplicates, kept first occurrence for each sample ID")
        
        # Convert sample names to match metadata index type
        # If metadata index is numeric and sample names are strings, convert sample names to int
        # If metadata index is string and sample names are numeric, convert to string
        try:
            if pd.api.types.is_integer_dtype(self.metadata.index):
                # Try to convert sample names to integers
                sample_names_converted = [int(name) for name in self.sample_names]
                self.metadata = self.metadata.reindex(sample_names_converted)
            elif pd.api.types.is_string_dtype(self.metadata.index):
                # Convert sample names to strings (they usually already are)
                sample_names_converted = [str(name) for name in self.sample_names]
                self.metadata = self.metadata.reindex(sample_names_converted)
            else:
                # Default behavior - use as is
                self.metadata = self.metadata.reindex(self.sample_names)
        except (ValueError, TypeError) as e:
            # If conversion fails, try converting metadata index to string to match sample names
            logging.warning(f"Could not convert sample names to match metadata index type: {e}")
            logging.info("Converting metadata index to string type for matching")
            self.metadata.index = self.metadata.index.astype(str)
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
    
    def validate_categorical_variable(self, var_data, variable_name, min_group_size=10, max_group_ratio=0.2):
        """
        Validate categorical variable suitability for PERMANOVA analysis.
        
        Args:
            var_data: Series of categorical data
            variable_name: Name of the variable
            min_group_size: Minimum samples per group (default: 10)
            max_group_ratio: Maximum ratio of groups to samples (default: 0.2)
            
        Returns:
            dict: Validation results with 'valid', 'reason', and 'stats' keys
        """
        valid_data = var_data.dropna()
        n_samples = len(valid_data)
        
        if n_samples < min_group_size * 2:
            return {
                'valid': False,
                'reason': f'Too few samples ({n_samples}) for meaningful grouping analysis',
                'stats': {'n_samples': n_samples, 'n_groups': 0}
            }
        
        # Get group counts
        group_counts = valid_data.value_counts()
        n_groups = len(group_counts)
        
        # Check group count ratio
        group_ratio = n_groups / n_samples
        if group_ratio > max_group_ratio:
            return {
                'valid': False,
                'reason': f'Too many groups ({n_groups}) relative to samples ({n_samples}). '
                         f'Ratio {group_ratio:.2f} exceeds maximum {max_group_ratio}',
                'stats': {'n_samples': n_samples, 'n_groups': n_groups, 'group_ratio': group_ratio}
            }
        
        # Check minimum group sizes
        small_groups = group_counts[group_counts < min_group_size]
        if len(small_groups) > 0:
            small_group_ratio = len(small_groups) / n_groups
            if small_group_ratio > 0.5:  # More than half the groups are too small
                return {
                    'valid': False,
                    'reason': f'{len(small_groups)} of {n_groups} groups have <{min_group_size} samples. '
                             f'Groups with small sizes: {dict(small_groups)}',
                    'stats': {'n_samples': n_samples, 'n_groups': n_groups, 'small_groups': len(small_groups)}
                }
        
        # Check for individual-level variables (like patient IDs)
        if n_groups > n_samples * 0.8:  # >80% of samples have unique values
            return {
                'valid': False,
                'reason': f'Variable appears to be individual-level identifier '
                         f'({n_groups} groups for {n_samples} samples). '
                         f'Not suitable for grouping analysis',
                'stats': {'n_samples': n_samples, 'n_groups': n_groups, 'uniqueness': n_groups/n_samples}
            }
        
        # Check for common individual-level identifier patterns
        variable_lower = variable_name.lower()
        individual_patterns = [
            'patient_id', 'subject_id', 'participant_id', 'sample_id', 'cambridge_patient_id',
            'id_', '_id', 'uuid', 'barcode', 'accession', 'identifier'
        ]
        
        if any(pattern in variable_lower for pattern in individual_patterns):
            if n_groups > n_samples * 0.5:  # >50% unique for ID-like variables
                return {
                    'valid': False,
                    'reason': f'Variable name suggests individual-level identifier '
                             f'({n_groups} groups for {n_samples} samples). '
                             f'Not suitable for grouping analysis',
                    'stats': {'n_samples': n_samples, 'n_groups': n_groups, 'uniqueness': n_groups/n_samples}
                }
        
        return {
            'valid': True,
            'reason': 'Variable passes validation checks',
            'stats': {
                'n_samples': n_samples,
                'n_groups': n_groups,
                'group_ratio': group_ratio,
                'min_group_size': group_counts.min(),
                'max_group_size': group_counts.max(),
                'mean_group_size': group_counts.mean()
            }
        }
    
    def _apply_multiple_testing_correction(self, results_df: pd.DataFrame, method: str = 'fdr_bh') -> pd.DataFrame:
        """
        Apply multiple testing correction to p-values.
        
        Args:
            results_df: DataFrame with p_value column
            method: Correction method ('fdr_bh', 'fdr_by', 'bonferroni', 'holm')
                   - 'fdr_bh': Benjamini-Hochberg FDR (default, less conservative)
                   - 'fdr_by': Benjamini-Yekutieli FDR (more conservative) 
                   - 'bonferroni': Bonferroni correction (most conservative)
                   - 'holm': Holm-Bonferroni (step-down method)
        
        Returns:
            DataFrame with added p_adjusted column
        """
        # Get valid p-values (exclude NaN)
        valid_mask = results_df['p_value'].notna()
        valid_p_values = results_df.loc[valid_mask, 'p_value'].values
        
        if len(valid_p_values) == 0:
            results_df['p_adjusted'] = np.nan
            return results_df
        
        # Apply correction using statsmodels if available, otherwise implement simple methods
        try:
            from statsmodels.stats.multitest import multipletests
            
            # Apply correction
            rejected, adjusted_p, alpha_sidak, alpha_bonf = multipletests(
                valid_p_values, alpha=0.05, method=method
            )
            
            # Add to dataframe
            results_df.loc[valid_mask, 'p_adjusted'] = adjusted_p
            results_df.loc[valid_mask, 'significant_after_correction'] = rejected
            
        except ImportError:
            # Fallback to manual implementation
            logging.warning("statsmodels not available. Using manual correction methods.")
            
            if method == 'bonferroni':
                # Bonferroni: p_adj = p * n_tests
                n_tests = len(valid_p_values)
                adjusted_p = np.minimum(valid_p_values * n_tests, 1.0)
                results_df.loc[valid_mask, 'p_adjusted'] = adjusted_p
                
            elif method == 'holm':
                # Holm-Bonferroni: step-down procedure
                n_tests = len(valid_p_values)
                sorted_indices = np.argsort(valid_p_values)
                sorted_p = valid_p_values[sorted_indices]
                
                # Apply Holm correction
                adjusted_p = np.zeros_like(sorted_p)
                for i in range(len(sorted_p)):
                    adjusted_p[i] = min(sorted_p[i] * (n_tests - i), 1.0)
                    if i > 0:
                        adjusted_p[i] = max(adjusted_p[i], adjusted_p[i-1])
                
                # Restore original order
                results_df.loc[valid_mask, 'p_adjusted'] = adjusted_p[np.argsort(sorted_indices)]
                
            else:
                # Default to Bonferroni if method not implemented
                logging.warning(f"Method '{method}' not implemented. Using Bonferroni correction.")
                n_tests = len(valid_p_values)
                adjusted_p = np.minimum(valid_p_values * n_tests, 1.0)
                results_df.loc[valid_mask, 'p_adjusted'] = adjusted_p
            
            # Mark significance
            results_df.loc[valid_mask, 'significant_after_correction'] = results_df.loc[valid_mask, 'p_adjusted'] < 0.05
        
        # Fill NaN for invalid p-values
        results_df.loc[~valid_mask, 'p_adjusted'] = np.nan
        results_df.loc[~valid_mask, 'significant_after_correction'] = False
        
        # Add correction method info
        results_df['correction_method'] = method
        
        return results_df
    
    def _adaptive_binning(self, values, min_bin_size=30, max_bins=5):
        """
        Create adaptive bins for continuous variables that:
        1. Ensure minimum sample size per bin for statistical power
        2. Maximize between-group variance
        3. Respect data distribution
        
        Args:
            values: Array of continuous values (no NaN)
            min_bin_size: Minimum samples per bin (default: 30)
            max_bins: Maximum number of bins (default: 5)
            
        Returns:
            Array of bin assignments (0, 1, 2, ...)
        """
        n_samples = len(values)
        
        # Determine optimal number of bins
        max_possible_bins = min(max_bins, n_samples // min_bin_size)
        if max_possible_bins < 2:
            # Not enough samples to bin meaningfully
            return np.zeros(n_samples, dtype=int)
        
        # Try different numbers of bins and pick the best
        best_bins = 2
        best_score = -1
        
        for n_bins in range(2, max_possible_bins + 1):
            # Create quantile-based bins
            quantiles = np.linspace(0, 100, n_bins + 1)
            bin_edges = np.percentile(values, quantiles)
            
            # Handle edge case where values are identical
            if len(np.unique(bin_edges)) != len(bin_edges):
                continue
                
            # Assign bins
            bin_assignments = np.digitize(values, bin_edges[1:-1])
            
            # Check minimum bin size requirement
            bin_counts = np.bincount(bin_assignments)
            if np.any(bin_counts < min_bin_size):
                continue
            
            # Calculate between-group variance score
            # Higher score = better separation
            total_var = np.var(values)
            between_group_var = 0
            
            for bin_idx in range(n_bins):
                mask = bin_assignments == bin_idx
                if np.sum(mask) > 0:
                    bin_mean = np.mean(values[mask])
                    bin_size = np.sum(mask)
                    between_group_var += bin_size * (bin_mean - np.mean(values)) ** 2
            
            between_group_var /= n_samples
            
            if total_var > 0:
                variance_ratio = between_group_var / total_var
                # Penalize too many bins to avoid overfitting
                score = variance_ratio - 0.05 * (n_bins - 2)
                
                if score > best_score:
                    best_score = score
                    best_bins = n_bins
        
        # Create final binning with optimal number of bins
        quantiles = np.linspace(0, 100, best_bins + 1)
        bin_edges = np.percentile(values, quantiles)
        
        # Ensure unique edges
        bin_edges = np.unique(bin_edges)
        
        # Final bin assignment
        bin_assignments = np.digitize(values, bin_edges[1:-1])
        
        return bin_assignments

    def generate_filtering_report(
        self, 
        exclusion_reasons: Dict[str, str], 
        included_variables: List[str],
        output_path: Optional[str] = None
    ) -> str:
        """Generate a detailed report of variable filtering decisions."""
        report_lines = [
            "# MetaGrouper Variable Filtering Report",
            "",
            f"**Analysis Date:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"**Total Variables in Metadata:** {len(self.metadata.columns)}",
            f"**Variables Included in Analysis:** {len(included_variables)}",
            f"**Variables Excluded:** {len(exclusion_reasons)}",
            "",
            "## Included Variables (Biologically Relevant)",
            ""
        ]
        
        for var in included_variables:
            if var in self.metadata.columns:
                values = self.metadata[var].dropna()
                unique_count = len(values.unique())
                missing_count = len(self.metadata) - len(values)
                
                # Determine variable type
                if pd.api.types.is_numeric_dtype(values):
                    var_type = "numerical"
                else:
                    var_type = "categorical"
                
                report_lines.append(f"- **{var}** ({var_type}): {unique_count} unique values, {missing_count} missing")
        
        if exclusion_reasons:
            report_lines.extend([
                "",
                "## Excluded Variables (Technical/Low Quality)",
                ""
            ])
            
            # Group exclusions by reason
            exclusion_groups = {}
            for var, reason in exclusion_reasons.items():
                if reason not in exclusion_groups:
                    exclusion_groups[reason] = []
                exclusion_groups[reason].append(var)
            
            for reason, vars_list in exclusion_groups.items():
                report_lines.append(f"### {reason}")
                for var in vars_list:
                    report_lines.append(f"- {var}")
                report_lines.append("")
        
        report_lines.extend([
            "## Recommendations",
            "",
            "- **High Priority Variables**: Focus on variables with biological significance (disease, demographics, genetics)",
            "- **Patient Grouping**: Consider including patient/subject IDs (PID, GEMM) for assembly grouping",
            "- **Temporal Analysis**: Include time-related variables (month, age) for longitudinal studies",
            "- **Technical Variables**: Exclude lab processing variables (Plate, Well, Barcode) unless needed for batch correction",
            ""
        ])
        
        report_text = "\n".join(report_lines)
        
        if output_path:
            with open(output_path, 'w') as f:
                f.write(report_text)
            logging.info(f"Variable filtering report saved to {output_path}")
        
        return report_text

    def analyze_variables(
        self, 
        variables: Optional[List[str]] = None, 
        n_permutations: int = 999,
        auto_filter: bool = False,
        exclude_variables: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """Analyze all or specified metadata variables using PERMANOVA."""
        if self.metadata is None:
            raise ValueError("Metadata not loaded. Call load_metadata() first.")

        # Apply smart filtering if requested
        exclusion_reasons = {}
        if auto_filter or variables is None:
            if variables is None:
                # No variables specified - use auto-filtering
                filtered_vars, exclusion_reasons = filter_metadata_variables(
                    self.metadata, 
                    auto_filter=True,
                    exclude_variables=exclude_variables
                )
                variables = filtered_vars
                
                # Log filtering results
                if exclusion_reasons:
                    logging.info(f"Auto-filtering excluded {len(exclusion_reasons)} variables:")
                    for var, reason in exclusion_reasons.items():
                        logging.debug(f"  Excluded '{var}': {reason}")
            else:
                # Variables specified but auto-filtering requested - validate them
                filtered_vars, exclusion_reasons = filter_metadata_variables(
                    self.metadata,
                    auto_filter=False,  # Don't auto-filter, just validate
                    include_variables=variables,
                    exclude_variables=exclude_variables
                )
                variables = filtered_vars
        else:
            # Use specified variables as-is
            if exclude_variables:
                excluded_vars = [v for v in variables if v in exclude_variables]
                variables = [v for v in variables if v not in exclude_variables]
                for var in excluded_vars:
                    exclusion_reasons[var] = "Manually excluded"

        logging.info(f"Analyzing {len(variables)} metadata variables")
        
        # Store filtering info for later use
        self.filtering_report_data = {
            'exclusion_reasons': exclusion_reasons,
            'included_variables': variables
        }

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
                # Validate categorical variable suitability
                validation = self.validate_categorical_variable(var_data, variable)
                
                if not validation['valid']:
                    logging.warning(f"Skipping variable '{variable}': {validation['reason']}")
                    # Still add to results but with warning flags
                    result = {
                        "variable": variable,
                        "variable_type": "categorical",
                        "f_statistic": np.nan,
                        "p_value": np.nan,
                        "r_squared": np.nan,
                        "missing_count": var_data.isna().sum(),
                        "validation_warning": validation['reason'],
                        "n_groups": validation['stats'].get('n_groups', 0),
                        "valid_samples": validation['stats'].get('n_samples', 0),
                        "min_group_size": validation['stats'].get('min_group_size', 0),
                        "max_group_size": validation['stats'].get('max_group_size', 0),
                        "mean_group_size": validation['stats'].get('mean_group_size', 0)
                    }
                    results.append(result)
                    self.permanova_results[variable] = result
                    continue
                
                # Check group sizes before analysis (existing validation)
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
                    # Use adaptive binning strategy
                    valid_values = var_array[~np.isnan(var_array)]
                    if (
                        len(np.unique(valid_values)) > 10
                    ):  # Only bin if many unique values
                        var_binned = self._adaptive_binning(valid_values, min_bin_size=30)
                        var_array = np.full_like(var_array, np.nan)
                        var_array[~np.isnan(var_array)] = var_binned
                        
                        # Log binning results
                        unique_bins = np.unique(var_binned[~np.isnan(var_binned)])
                        bin_counts = [np.sum(var_binned == b) for b in unique_bins]
                        logging.info(f"Binned {variable} into {len(unique_bins)} groups: {bin_counts} samples each")

            # Run PERMANOVA
            result = permanova.permanova_test(var_array, n_permutations)
            result["variable"] = variable
            result["variable_type"] = (
                "categorical"
                if self.metadata[variable].dtype == "object"
                else "numerical"
            )
            result["missing_count"] = self.metadata[variable].isna().sum()
            
            # Add group information
            valid_groups = var_array[~np.isnan(var_array)]
            if len(valid_groups) > 0:
                unique_groups = np.unique(valid_groups)
                result["n_groups"] = len(unique_groups)
                result["valid_samples"] = len(valid_groups)
                
                # Add group size distribution for better interpretation
                group_sizes = [np.sum(valid_groups == g) for g in unique_groups]
                result["min_group_size"] = min(group_sizes)
                result["max_group_size"] = max(group_sizes)
                result["mean_group_size"] = np.mean(group_sizes)
            else:
                result["n_groups"] = 0
                result["valid_samples"] = 0
                result["min_group_size"] = 0
                result["max_group_size"] = 0
                result["mean_group_size"] = 0

            results.append(result)
            self.permanova_results[variable] = result

        # Create results DataFrame
        results_df = pd.DataFrame(results)
        if not results_df.empty:
            # Apply multiple testing correction
            results_df = self._apply_multiple_testing_correction(results_df)
            
            # Sort by adjusted R² if available, otherwise by regular R²
            sort_column = "r_squared_adj" if "r_squared_adj" in results_df.columns else "r_squared"
            results_df = results_df.sort_values(sort_column, ascending=False, na_position='last')

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

        # Ensure distance matrix is symmetric (fix floating point differences)
        symmetric_distance_matrix = (self.distance_matrix + self.distance_matrix.T) / 2
        np.fill_diagonal(symmetric_distance_matrix, 0)  # Ensure diagonal is exactly zero

        mds = MDS(
            n_components=min(10, len(self.sample_names) - 1),
            dissimilarity="precomputed",
            random_state=42,
        )
        X = mds.fit_transform(symmetric_distance_matrix)

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
                    labels = clusterer.fit_predict(symmetric_distance_matrix)
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
            "Variables ranked by adjusted R-squared (penalized for degrees of freedom):\n\n"
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

                # Show effect sizes
                r_squared_text = f"R² = {row['r_squared']:.3f}"
                if 'r_squared_adj' in row and not pd.isna(row['r_squared_adj']):
                    r_squared_text += f", R²_adj = {row['r_squared_adj']:.3f}"
                if 'omega_squared' in row and not pd.isna(row['omega_squared']):
                    r_squared_text += f", ω² = {row['omega_squared']:.3f}"
                
                # Add effect size
                effect_size_text = ""
                if 'effect_size' in row and not pd.isna(row['effect_size']):
                    effect_size_text = f" ({row['effect_size']} effect)"
                
                # Add adjusted p-value if available
                p_value_text = f"p = {row['p_value']:.3f}"
                if 'p_adjusted' in row and not pd.isna(row['p_adjusted']):
                    p_value_text += f" (p_adj = {row['p_adjusted']:.3f})"
                
                f.write(
                    f"- **{row['variable']}**: {r_squared_text}, "
                    f"{p_value_text}{significance}{effect_size_text}\n"
                )
                f.write(f"  - Type: {row['variable_type']}\n")
                f.write(f"  - Valid samples: {row.get('valid_samples', row['n_samples'])}\n")
                f.write(f"  - Groups: {row.get('n_groups', 'N/A')}\n")
                
                # Add group size information if available
                if 'min_group_size' in row and row['min_group_size'] > 0:
                    f.write(f"  - Group sizes: {row['min_group_size']}-{row['max_group_size']} "
                           f"(mean: {row['mean_group_size']:.1f})\n")
                
                # Add power warning if present
                if 'power_warning' in row and row['power_warning'] is not None:
                    f.write(f"  - ⚠️  {row['power_warning']}\n")
                
                # Add PERMDISP results
                if 'permdisp_p' in row and not pd.isna(row['permdisp_p']):
                    f.write(f"  - PERMDISP: F={row['permdisp_f']:.3f}, p={row['permdisp_p']:.3f}\n")
                
                # Add variance warning if present
                if 'variance_warning' in row and row['variance_warning'] is not None:
                    f.write(f"  - ⚠️  {row['variance_warning']}\n")
                
                f.write("\n")
            
            # Add warnings section for invalid variables
            warning_results = results_df[results_df['validation_warning'].notna()] if 'validation_warning' in results_df.columns else pd.DataFrame()
            if not warning_results.empty:
                f.write("## Variables with Statistical Concerns\n\n")
                f.write("The following variables were flagged with statistical concerns:\n\n")
                for _, row in warning_results.iterrows():
                    f.write(f"- **{row['variable']}**: {row['validation_warning']}\n")
                f.write("\n")
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
            "- **R-squared_adj**: Adjusted R-squared that penalizes for high numbers of groups\n"
        )
        f.write(
            "- **Omega squared (ω²)**: Less biased effect size measure (preferred over R²)\n"
        )
        f.write(
            "- **Effect size**: Based on ω² (negligible <0.01, small <0.06, medium <0.14, large ≥0.14)\n"
        )
        f.write(
            "- **p-value**: Statistical significance (< 0.05 is typically significant)\n"
        )
        f.write(
            "- **p_adjusted**: P-value after multiple testing correction (FDR/Bonferroni)\n"
        )
        f.write(
            "- **Group sizes**: Range of samples per group (minimum 10 recommended for reliable results)\n"
        )
        f.write(
            "- **PERMDISP**: Test for homogeneity of dispersions (p<0.05 indicates heterogeneous variances)\n"
        )
        f.write(
            "- **Silhouette score**: Quality of clustering (higher is better, > 0.5 is good)\n"
        )
        f.write("- Significance codes: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1\n")
        f.write("\n**Note**: Variables with >20% of samples as unique groups may show inflated R² values.\n")
        f.write("**Power**: Variables with <10 samples per group may have unreliable results.\n")
        f.write("**Assumptions**: PERMANOVA assumes homogeneous dispersions. Check PERMDISP results.\n")

    logging.info(f"Summary report saved to {output_path}")
