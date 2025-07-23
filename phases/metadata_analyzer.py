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
    max_unique_ratio: float = 0.5,
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
            # Check if high uniqueness is informative (e.g., patient IDs) vs noise
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

    def permdisp_test(
        self, metadata_variable: np.ndarray, n_permutations: int = 999
    ) -> Dict[str, float]:
        """
        Perform PERMDISP test for homogeneity of dispersions.
        Tests if groups have equal variances in multivariate space.
        """
        # Remove samples with missing metadata
        valid_indices = ~pd.isna(metadata_variable)
        if not np.any(valid_indices):
            return {"f_statistic": np.nan, "p_value": np.nan}
            
        valid_distance_matrix = self.distance_matrix[
            np.ix_(valid_indices, valid_indices)
        ]
        
        # Ensure distance matrix is symmetric (fix floating point differences)
        valid_distance_matrix = (valid_distance_matrix + valid_distance_matrix.T) / 2
        np.fill_diagonal(valid_distance_matrix, 0)  # Ensure diagonal is exactly zero
        
        valid_metadata = metadata_variable[valid_indices]
        groups = np.unique(valid_metadata)
        
        if len(groups) < 2:
            return {"f_statistic": np.nan, "p_value": np.nan}
            
        # Calculate group centroids using PCoA
        from sklearn.decomposition import PCA
        from sklearn.manifold import MDS
        from scipy.spatial.distance import squareform, pdist
        
        # Convert distance matrix to embedding
        n_components = min(len(valid_metadata) - 1, 10)
        mds = MDS(n_components=n_components, dissimilarity='precomputed', random_state=42)
        embedding = mds.fit_transform(valid_distance_matrix)
        
        # Calculate distances to group centroids
        group_dispersions = {}
        for group in groups:
            group_mask = valid_metadata == group
            if np.sum(group_mask) < 2:
                continue
                
            # Calculate centroid
            group_embedding = embedding[group_mask]
            centroid = np.mean(group_embedding, axis=0)
            
            # Calculate distances to centroid
            distances = np.sqrt(np.sum((group_embedding - centroid) ** 2, axis=1))
            group_dispersions[group] = distances
            
        # Perform F-test on dispersions
        all_dispersions = []
        group_labels = []
        for group, dispersions in group_dispersions.items():
            all_dispersions.extend(dispersions)
            group_labels.extend([group] * len(dispersions))
            
        all_dispersions = np.array(all_dispersions)
        group_labels = np.array(group_labels)
        
        # Calculate F-statistic
        observed_f = self._calculate_f_statistic_dispersions(all_dispersions, group_labels)
        
        # Permutation test
        f_permuted = []
        for _ in range(n_permutations):
            perm_labels = np.random.permutation(group_labels)
            f_perm = self._calculate_f_statistic_dispersions(all_dispersions, perm_labels)
            f_permuted.append(f_perm)
            
        f_permuted = np.array(f_permuted)
        p_value = np.sum(f_permuted >= observed_f) / n_permutations
        
        return {"f_statistic": observed_f, "p_value": p_value}
    
    def _calculate_f_statistic_dispersions(self, dispersions: np.ndarray, groups: np.ndarray) -> float:
        """Calculate F-statistic for testing equality of dispersions."""
        unique_groups = np.unique(groups)
        
        # Calculate group means and overall mean
        group_means = {}
        for group in unique_groups:
            group_mask = groups == group
            group_means[group] = np.mean(dispersions[group_mask])
            
        overall_mean = np.mean(dispersions)
        
        # Calculate between-group and within-group sum of squares
        ss_between = 0
        ss_within = 0
        
        for group in unique_groups:
            group_mask = groups == group
            n_group = np.sum(group_mask)
            
            # Between-group SS
            ss_between += n_group * (group_means[group] - overall_mean) ** 2
            
            # Within-group SS
            ss_within += np.sum((dispersions[group_mask] - group_means[group]) ** 2)
            
        # Calculate degrees of freedom
        df_between = len(unique_groups) - 1
        df_within = len(dispersions) - len(unique_groups)
        
        # Calculate F-statistic
        if df_within > 0 and ss_within > 0:
            ms_between = ss_between / df_between
            ms_within = ss_within / df_within
            f_statistic = ms_between / ms_within
        else:
            f_statistic = np.nan
            
        return f_statistic

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
        
        # Ensure distance matrix is symmetric (fix floating point differences)
        valid_distance_matrix = (valid_distance_matrix + valid_distance_matrix.T) / 2
        np.fill_diagonal(valid_distance_matrix, 0)  # Ensure diagonal is exactly zero
        
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
        """Load metadata from CSV/TSV file with auto-detection of sample ID column."""
        logging.info(f"Loading metadata from {metadata_file}")

        # Detect file format
        if metadata_file.endswith(".tsv") or metadata_file.endswith(".txt"):
            self.metadata = pd.read_csv(metadata_file, sep="\t")
        else:
            self.metadata = pd.read_csv(metadata_file)

        # Auto-detect sample ID column if the specified one doesn't exist or doesn't match well
        actual_sample_id_column = self._auto_detect_sample_id_column(sample_id_column)
        
        # Align metadata with sample names
        self.metadata = self.metadata.set_index(actual_sample_id_column)
        
        # Ensure index is string type to match sample_names
        self.metadata.index = self.metadata.index.astype(str)
        
        # Check for duplicate sample IDs and handle them
        if self.metadata.index.duplicated().any():
            duplicates = self.metadata.index[self.metadata.index.duplicated(keep=False)]
            logging.warning(f"Found {len(duplicates)} duplicate sample IDs in metadata: {list(duplicates.unique())}")
            
            # Remove duplicates, keeping the first occurrence
            self.metadata = self.metadata[~self.metadata.index.duplicated(keep='first')]
            logging.info(f"Removed duplicates, kept first occurrence for each sample ID")
        
        # Ensure sample_names are also strings
        sample_names_str = [str(name) for name in self.sample_names]
        self.metadata = self.metadata.reindex(sample_names_str)

        logging.info(
            f"Loaded metadata for {len(self.metadata)} samples with "
            f"{len(self.metadata.columns)} variables"
        )

        # Report missing data
        missing_samples = self.metadata.index[self.metadata.isnull().all(axis=1)]
        if len(missing_samples) > 0:
            logging.warning(f"Missing metadata for samples: {list(missing_samples)}")
    
    def _auto_detect_sample_id_column(self, preferred_column: str = "sample_id"):
        """Auto-detect the best sample ID column using the same logic as Interactive Report."""
        logging.info(f"Auto-detecting sample ID column (preferred: {preferred_column})")
        
        # First, try the preferred column if it exists and has good matches
        if preferred_column in self.metadata.columns:
            col_values = self.metadata[preferred_column].astype(str).values
            sample_names_str = [str(name) for name in self.sample_names]
            exact_matches = sum(1 for name in sample_names_str if name in col_values)
            
            if exact_matches > len(self.sample_names) * 0.5:  # If >50% match
                logging.info(f"Using preferred column '{preferred_column}' with {exact_matches}/{len(self.sample_names)} matches")
                return preferred_column
        
        # Try common column names
        possible_columns = ['sample_id', 'sample', 'accession', 'run_id', 'srr', 'sample_name', 'id', 'sra_accession', 'database_ID']
        
        logging.info(f"Checking common column names: {possible_columns}")
        for col in possible_columns:
            if col in self.metadata.columns:
                # Convert both to string and check for matches
                col_values_str = self.metadata[col].astype(str).values
                sample_names_str = [str(name) for name in self.sample_names]
                
                # Try exact matches first
                exact_matches = sum(1 for name in sample_names_str if name in col_values_str)
                logging.info(f"   {col}: {exact_matches}/{len(self.sample_names)} exact matches")
                
                if exact_matches > 0:  # Any match is good enough for common columns
                    logging.info(f"✅ Using column '{col}' as sample identifier")
                    return col
        
        # Try all columns as a last resort
        logging.info("Checking all columns for any matches...")
        for col in self.metadata.columns:
            col_values = self.metadata[col].astype(str).values
            sample_names_str = [str(name) for name in self.sample_names]
            matches = sum(1 for name in sample_names_str if name in col_values)
            
            if matches > len(self.sample_names) * 0.5:  # If >50% match
                logging.info(f"   {col}: {matches} potential matches found")
                logging.info(f"✅ Using column '{col}' as sample identifier")
                return col
        
        # If no good matches found, raise an error
        raise ValueError(
            f"Could not find a suitable sample ID column. Checked columns: {list(self.metadata.columns)}. "
            f"Sample names: {self.sample_names[:5]}..."
        )

    def validate_sample_size(self, groups):
        """Validate that groups have sufficient sample size for reliable PERMANOVA results."""
        for group_name, group_data in groups.items():
            if len(group_data) < 10:
                warnings.warn(f"Group '{group_name}' has only {len(group_data)} samples. "
                             f"PERMANOVA requires ≥10 samples per group for reliable results.")

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

            # Run PERMDISP first to check homogeneity assumption
            permdisp_result = permanova.permdisp_test(var_array, n_permutations)
            
            # Run PERMANOVA
            result = permanova.permanova_test(var_array, n_permutations)
            result["variable"] = variable
            result["variable_type"] = (
                "categorical"
                if self.metadata[variable].dtype == "object"
                else "numerical"
            )
            
            # Add PERMDISP results
            result["permdisp_p_value"] = permdisp_result["p_value"]
            result["permdisp_f_statistic"] = permdisp_result["f_statistic"]
            result["homogeneity_violated"] = permdisp_result["p_value"] < 0.05
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
