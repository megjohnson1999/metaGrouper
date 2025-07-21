#!/usr/bin/env python3
"""
MetaGrouper Phase 3: Assembly Strategy Recommendation Engine

This module analyzes k-mer similarities and metadata associations to recommend
optimal assembly strategies including individual assembly, group co-assembly,
and global co-assembly with specific grouping criteria and confidence scores.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Union
import logging
from pathlib import Path
import json
from dataclasses import dataclass, asdict
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score, calinski_harabasz_score
import matplotlib.pyplot as plt
import seaborn as sns


@dataclass
class AssemblyGroup:
    """Represents a group of samples recommended for co-assembly."""

    group_id: str
    sample_names: List[str]
    grouping_criterion: str
    criterion_value: Any
    avg_distance: float
    max_distance: float
    confidence_score: float
    expected_benefits: List[str]
    expected_challenges: List[str]


@dataclass
class ConfidenceBreakdown:
    """Detailed breakdown of confidence score components."""
    total_confidence: float
    isolation_score: Optional[float] = None  # For individual assembly
    clustering_score: Optional[float] = None  # For group assembly  
    stability_score: Optional[float] = None  # Group stability
    metadata_alignment: Optional[float] = None  # Metadata consistency
    statistical_significance: Optional[float] = None  # P-values, effect sizes
    sample_size_adequacy: Optional[float] = None  # Sufficient samples for reliable grouping
    explanation: str = ""


@dataclass
class AssemblyRecommendation:
    """Complete assembly strategy recommendation."""

    strategy: str  # 'individual', 'grouped', 'global'
    groups: List[AssemblyGroup]
    overall_confidence: float
    confidence_breakdown: Optional[ConfidenceBreakdown]
    primary_criterion: Optional[str]
    decision_rationale: str
    assembly_commands: Dict[str, List[str]]
    performance_predictions: Dict[str, Any]


class ConfidenceCalculator:
    """Data-driven confidence calculation for assembly recommendations."""
    
    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str], metadata: Optional[pd.DataFrame] = None):
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names
        self.metadata = metadata
        self.n_samples = len(sample_names)
        
    def calculate_individual_confidence(self, sample_idx: int) -> ConfidenceBreakdown:
        """Calculate confidence for individual assembly of a specific sample."""
        # Isolation score: How different this sample is from all others
        sample_distances = self.distance_matrix[sample_idx, :]
        other_distances = np.delete(sample_distances, sample_idx)
        
        if len(other_distances) == 0:
            isolation_score = 1.0  # Only sample, definitely individual
        else:
            # Higher mean distance = more isolated = higher confidence for individual assembly
            mean_distance_to_others = np.mean(other_distances)
            min_distance_to_others = np.min(other_distances)
            
            # Calculate isolation score based on dataset distribution
            # Use dataset statistics rather than arbitrary thresholds
            all_distances = self.distance_matrix[np.triu_indices(self.n_samples, k=1)]
            median_distance = np.median(all_distances)
            percentile_75 = np.percentile(all_distances, 75)
            
            # Score based on how much more isolated this sample is than typical
            # Higher than median = some isolation, higher than 75th percentile = strong isolation
            if median_distance > 0:
                isolation_score = min(mean_distance_to_others / percentile_75, 1.0)
            else:
                isolation_score = 0.5  # All samples identical
            
            # Additional boost for samples with no close neighbors
            if min_distance_to_others > median_distance:
                isolation_boost = min_distance_to_others / median_distance
                isolation_score = min(isolation_score * (1.0 + isolation_boost * 0.2), 1.0)
        
        # Metadata uniqueness (if available)
        metadata_alignment = None
        if self.metadata is not None and sample_idx < len(self.metadata):
            # Count how many metadata features are unique or rare for this sample
            sample_row = self.metadata.iloc[sample_idx]
            uniqueness_score = 0
            total_features = 0
            
            for col in self.metadata.columns:
                if col != 'sample_id' and not pd.isna(sample_row[col]):
                    value_counts = self.metadata[col].value_counts()
                    if len(value_counts) > 1:  # Skip constant columns
                        frequency = value_counts.get(sample_row[col], 0) / len(self.metadata)
                        uniqueness_score += (1 - frequency)  # Rare values get higher scores
                        total_features += 1
            
            metadata_alignment = uniqueness_score / max(total_features, 1)
        
        # Combine scores
        if metadata_alignment is not None:
            total_confidence = (isolation_score * 0.7 + metadata_alignment * 0.3)
        else:
            total_confidence = isolation_score
            
        explanation = f"Sample is {'well-isolated' if isolation_score > 0.6 else 'moderately isolated'} from other samples (isolation={isolation_score:.2f})"
        if metadata_alignment is not None:
            explanation += f", with {'unique' if metadata_alignment > 0.5 else 'common'} metadata profile (uniqueness={metadata_alignment:.2f})"
            
        return ConfidenceBreakdown(
            total_confidence=total_confidence,
            isolation_score=isolation_score,
            metadata_alignment=metadata_alignment,
            explanation=explanation
        )
    
    def calculate_group_confidence(self, group_indices: List[int], group_label: str = "") -> ConfidenceBreakdown:
        """Calculate confidence for co-assembly of a group of samples."""
        if len(group_indices) < 2:
            return ConfidenceBreakdown(total_confidence=0.0, explanation="Group too small for co-assembly")
            
        # Extract group distance matrix
        group_distances = self.distance_matrix[np.ix_(group_indices, group_indices)]
        
        # Clustering tightness - use silhouette-like score
        if len(group_indices) >= 2:
            # Compare intra-group distances to distances to rest of samples
            intra_group_distances = []
            inter_group_distances = []
            
            for i, idx1 in enumerate(group_indices):
                for j, idx2 in enumerate(group_indices):
                    if i != j:
                        intra_group_distances.append(self.distance_matrix[idx1, idx2])
                        
                # Distances to samples NOT in this group
                for other_idx in range(self.n_samples):
                    if other_idx not in group_indices:
                        inter_group_distances.append(self.distance_matrix[idx1, other_idx])
            
            if intra_group_distances and inter_group_distances:
                mean_intra = np.mean(intra_group_distances)
                mean_inter = np.mean(inter_group_distances)
                
                # Silhouette-like score: want low intra-group, high inter-group distances
                if mean_inter > 0:
                    clustering_score = (mean_inter - mean_intra) / max(mean_inter, mean_intra)
                    clustering_score = max(0, clustering_score)  # Clamp to [0,1]
                else:
                    clustering_score = 0.0
            else:
                clustering_score = 0.0
        else:
            clustering_score = 0.0
        
        # Sample size adequacy - adaptive scoring based on dataset size and assembly theory
        # Optimal group size scales with total dataset size and computational constraints
        total_samples = self.n_samples
        
        # Calculate adaptive optimal range based on dataset size
        if total_samples <= 20:
            # Small datasets: be more permissive with group sizes
            optimal_min, optimal_max = 2, max(6, total_samples // 3)
        elif total_samples <= 100:
            # Medium datasets: classic co-assembly range
            optimal_min, optimal_max = 3, min(15, total_samples // 5)
        else:
            # Large datasets: larger groups become more beneficial
            optimal_min, optimal_max = 5, min(25, total_samples // 8)
        
        group_size = len(group_indices)
        
        # Adaptive scoring function - bell curve centered on optimal range
        if group_size < optimal_min:
            # Too small - limited diversity
            size_score = 0.6 + 0.3 * (group_size / optimal_min)
        elif optimal_min <= group_size <= optimal_max:
            # Optimal range - highest scores
            size_score = 1.0
        else:
            # Too large - diminishing returns due to complexity
            excess = group_size - optimal_max
            max_excess = max(10, optimal_max)  # Allow some flexibility
            size_score = max(0.5, 1.0 - (excess / max_excess) * 0.4)
        
        # Metadata consistency (if available)
        metadata_alignment = None
        if self.metadata is not None and len(group_indices) >= 2:
            group_metadata = self.metadata.iloc[group_indices]
            consistency_scores = []
            
            for col in group_metadata.columns:
                if col != 'sample_id' and not group_metadata[col].isna().all():
                    # For categorical: higher consistency = more samples with same value
                    if group_metadata[col].dtype == 'object':
                        mode_count = group_metadata[col].value_counts().iloc[0] if len(group_metadata[col].value_counts()) > 0 else 0
                        consistency = mode_count / len(group_indices)
                        consistency_scores.append(consistency)
                    else:
                        # For numerical: lower coefficient of variation = higher consistency  
                        if group_metadata[col].std() > 0:
                            cv = group_metadata[col].std() / group_metadata[col].mean() if group_metadata[col].mean() != 0 else 1
                            consistency = max(0, 1 - min(cv, 1))  # Convert CV to consistency score
                            consistency_scores.append(consistency)
            
            metadata_alignment = np.mean(consistency_scores) if consistency_scores else 0.5
        
        # Combine scores
        scores = [clustering_score, size_score]
        weights = [0.6, 0.4]
        
        if metadata_alignment is not None:
            scores.append(metadata_alignment)
            weights = [0.5, 0.3, 0.2]  # Reweight
        
        total_confidence = np.average(scores, weights=weights)
        
        # Generate explanation
        explanation = f"Group of {len(group_indices)} samples with "
        explanation += f"{'tight' if clustering_score > 0.6 else 'loose' if clustering_score > 0.3 else 'weak'} clustering (score={clustering_score:.2f}), "
        explanation += f"{'optimal' if size_score >= 0.9 else 'adequate' if size_score >= 0.7 else 'suboptimal'} size (score={size_score:.2f})"
        if metadata_alignment is not None:
            explanation += f", {'consistent' if metadata_alignment > 0.6 else 'mixed'} metadata (score={metadata_alignment:.2f})"
        
        return ConfidenceBreakdown(
            total_confidence=total_confidence,
            clustering_score=clustering_score,
            sample_size_adequacy=size_score,
            metadata_alignment=metadata_alignment,
            explanation=explanation
        )
    
    def calculate_overall_strategy_confidence(self, strategy: str, groups: List[List[int]]) -> ConfidenceBreakdown:
        """Calculate overall confidence for the chosen strategy."""
        if strategy == "individual":
            # For individual strategy, confidence is based on how well samples are separated
            if self.n_samples <= 1:
                return ConfidenceBreakdown(total_confidence=1.0, explanation="Single sample - individual assembly certain")
            
            # Calculate average isolation of all samples
            individual_scores = []
            for i in range(self.n_samples):
                individual_conf = self.calculate_individual_confidence(i)
                individual_scores.append(individual_conf.total_confidence)
            
            overall_score = np.mean(individual_scores)
            explanation = f"Individual assembly recommended: average sample isolation = {overall_score:.2f}"
            
            return ConfidenceBreakdown(
                total_confidence=overall_score,
                isolation_score=overall_score,
                explanation=explanation
            )
            
        elif strategy == "grouped":
            # For grouped strategy, confidence is average of group confidences
            group_confidences = []
            for group_indices in groups:
                group_conf = self.calculate_group_confidence(group_indices)
                group_confidences.append(group_conf.total_confidence)
            
            if group_confidences:
                overall_score = np.mean(group_confidences)
                explanation = f"Grouped assembly: {len(groups)} groups with average confidence {overall_score:.2f}"
            else:
                overall_score = 0.0
                explanation = "No valid groups found for grouped assembly"
                
            return ConfidenceBreakdown(
                total_confidence=overall_score,
                clustering_score=overall_score,
                explanation=explanation
            )
            
        elif strategy == "global":
            # For global strategy, treat all samples as one big group
            all_indices = list(range(self.n_samples))
            global_conf = self.calculate_group_confidence(all_indices, "global")
            global_conf.explanation = f"Global co-assembly: all {self.n_samples} samples together - " + global_conf.explanation
            return global_conf
            
        else:
            return ConfidenceBreakdown(total_confidence=0.0, explanation=f"Unknown strategy: {strategy}")


class AssemblyStrategyEngine:
    """Core engine for determining optimal assembly strategies."""

    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str], metadata: Optional[pd.DataFrame] = None):
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names
        self.metadata = metadata
        self.n_samples = len(sample_names)

        # Configurable thresholds for assembly decisions
        self.similarity_threshold_high = 0.30   # Stringent grouping
        self.similarity_threshold_medium = 0.45  # Moderate grouping (default)
        self.significance_threshold = 0.05  # P-value threshold for metadata
        self.min_group_size = 2
        self.max_group_size = 20
        
        # Initialize confidence calculator
        self.confidence_calc = ConfidenceCalculator(distance_matrix, sample_names, metadata)

    def _calculate_group_statistics(
        self, group_indices: List[int]
    ) -> Tuple[float, float]:
        """Calculate within-group distance statistics."""
        if len(group_indices) < 2:
            return 0.0, 0.0

        group_distances = []
        for i in range(len(group_indices)):
            for j in range(i + 1, len(group_indices)):
                dist = self.distance_matrix[group_indices[i], group_indices[j]]
                group_distances.append(dist)

        avg_distance = np.mean(group_distances)
        max_distance = np.max(group_distances)

        return avg_distance, max_distance

    def _assess_group_quality(self, group_indices: List[int]) -> float:
        """Assess the quality of a proposed assembly group using data-driven confidence calculation."""
        confidence_breakdown = self.confidence_calc.calculate_group_confidence(group_indices)
        return confidence_breakdown.total_confidence

    def recommend_by_similarity(self) -> List[AssemblyGroup]:
        """Recommend assembly groups based purely on k-mer similarity."""
        logging.info("Generating similarity-based assembly recommendations")

        # Use hierarchical clustering to identify natural groups
        condensed_distances = squareform(self.distance_matrix, checks=False)
        linkage_matrix = linkage(condensed_distances, method="average")

        groups = []

        # Try different numbers of clusters
        for n_clusters in range(2, min(self.n_samples, 8)):
            cluster_labels = fcluster(linkage_matrix, n_clusters, criterion="maxclust")

            for cluster_id in range(1, n_clusters + 1):
                cluster_indices = [
                    i for i, label in enumerate(cluster_labels) if label == cluster_id
                ]

                if len(cluster_indices) < self.min_group_size:
                    continue

                avg_dist, max_dist = self._calculate_group_statistics(cluster_indices)

                # Only recommend groups with reasonable similarity
                if avg_dist <= self.similarity_threshold_medium:
                    confidence = self._assess_group_quality(cluster_indices)

                    group = AssemblyGroup(
                        group_id=f"similarity_cluster_{cluster_id}_{n_clusters}",
                        sample_names=[self.sample_names[i] for i in cluster_indices],
                        grouping_criterion="k-mer_similarity",
                        criterion_value=f"avg_distance_{avg_dist:.3f}",
                        avg_distance=avg_dist,
                        max_distance=max_dist,
                        confidence_score=confidence,
                        expected_benefits=[
                            "Improved assembly continuity",
                            "Better coverage of shared sequences",
                            "Enhanced detection of strain variants",
                        ],
                        expected_challenges=[
                            "Potential strain mixing",
                            "Increased computational requirements",
                            f"Assembly size may be {len(cluster_indices)}x larger",
                        ],
                    )
                    groups.append(group)

        # Sort by confidence and remove overlapping groups
        groups.sort(key=lambda x: x.confidence_score, reverse=True)

        # Remove overlapping groups (keep highest confidence)
        final_groups = []
        used_samples = set()

        for group in groups:
            group_samples = set(group.sample_names)
            if not group_samples.intersection(used_samples):
                final_groups.append(group)
                used_samples.update(group_samples)

        return final_groups

    def generate_metadata_grouping_recommendations(
        self, metadata_results: pd.DataFrame, metadata: pd.DataFrame
    ) -> List[Dict[str, Any]]:
        """Generate intelligent 'group by X' recommendations based on metadata analysis."""
        logging.info("Generating metadata grouping recommendations")
        
        if metadata_results.empty or metadata is None:
            return []
        
        recommendations = []
        
        # Focus on significant and explanatory variables
        significant_vars = metadata_results[
            (metadata_results["p_value"] < self.significance_threshold)
            | (metadata_results["r_squared"] > 0.15)  # Include moderately explanatory variables
        ].sort_values("r_squared", ascending=False)
        
        for _, row in significant_vars.iterrows():
            variable = row["variable"]
            
            if variable not in metadata.columns:
                continue
            
            # Analyze grouping potential for this variable
            var_data = metadata[variable].dropna()
            
            if len(var_data) < 2:
                continue
            
            unique_values = var_data.unique()
            value_counts = var_data.value_counts()
            
            # Calculate grouping statistics
            total_samples = len(var_data)
            num_groups = len(unique_values)
            group_sizes = value_counts.tolist()
            min_group_size = min(group_sizes)
            max_group_size = max(group_sizes)
            avg_group_size = total_samples / num_groups
            
            # Filter out groups that are too small
            viable_groups = [size for size in group_sizes if size >= self.min_group_size]
            num_viable_groups = len(viable_groups)
            samples_in_viable_groups = sum(viable_groups)
            
            if num_viable_groups < 2:  # Need at least 2 viable groups
                continue
            
            # Calculate confidence based on multiple factors
            r_squared = row["r_squared"]
            p_value = row.get("p_value", 1.0)
            
            # Confidence factors
            statistical_confidence = min(r_squared * 2, 1.0)  # R² contribution
            size_balance_factor = min(min_group_size / max_group_size, 1.0)  # Group size balance
            coverage_factor = samples_in_viable_groups / total_samples  # Sample coverage
            
            overall_confidence = (statistical_confidence * 0.5 + 
                                size_balance_factor * 0.2 + 
                                coverage_factor * 0.3)
            
            # Generate recommendation rationale
            rationale_parts = []
            
            if r_squared > 0.3:
                rationale_parts.append(f"Strong association (R² = {r_squared:.3f})")
            elif r_squared > 0.15:
                rationale_parts.append(f"Moderate association (R² = {r_squared:.3f})")
            
            if p_value < 0.001:
                rationale_parts.append("highly significant (p < 0.001)")
            elif p_value < 0.01:
                rationale_parts.append("very significant (p < 0.01)")
            elif p_value < 0.05:
                rationale_parts.append("significant (p < 0.05)")
            
            rationale_parts.append(f"creates {num_viable_groups} viable groups")
            rationale_parts.append(f"covers {samples_in_viable_groups}/{total_samples} samples")
            
            # Generate benefits and challenges
            benefits = []
            challenges = []
            
            if avg_group_size >= 3 and avg_group_size <= 10:
                benefits.append("Optimal group sizes for co-assembly")
            elif avg_group_size > 10:
                benefits.append("Large groups may improve assembly contiguity")
                challenges.append("Large groups may increase computational requirements")
            else:
                challenges.append("Small groups may limit co-assembly benefits")
            
            if r_squared > 0.3:
                benefits.append("Strong biological basis for grouping")
            
            if size_balance_factor > 0.7:
                benefits.append("Well-balanced group sizes")
            else:
                challenges.append("Uneven group size distribution")
            
            # Determine recommended strategy
            if avg_group_size <= 2:
                strategy = "individual"
                strategy_note = "Groups too small for effective co-assembly"
            elif avg_group_size <= 15 and overall_confidence > 0.6:
                strategy = "grouped_coassembly"
                strategy_note = "Recommended for group-wise co-assembly"
            elif num_viable_groups <= 3 and total_samples > 20:
                strategy = "grouped_coassembly"
                strategy_note = "Large groups suitable for co-assembly"
            else:
                strategy = "individual"
                strategy_note = "Consider individual assembly due to complexity"
            
            recommendation = {
                "variable": variable,
                "strategy": strategy,
                "confidence": overall_confidence,
                "r_squared": r_squared,
                "p_value": p_value,
                "num_groups": num_viable_groups,
                "total_samples": total_samples,
                "samples_in_viable_groups": samples_in_viable_groups,
                "group_sizes": viable_groups,
                "min_group_size": min_group_size,
                "max_group_size": max_group_size,
                "avg_group_size": avg_group_size,
                "rationale": " - ".join(rationale_parts),
                "strategy_note": strategy_note,
                "benefits": benefits,
                "challenges": challenges,
                "recommendation_text": f"Group by '{variable}': {strategy_note} "
                                     f"(confidence: {overall_confidence:.2f}, "
                                     f"R² = {r_squared:.3f}, "
                                     f"{num_viable_groups} groups of {min_group_size}-{max_group_size} samples)"
            }
            
            recommendations.append(recommendation)
        
        # Sort by confidence, then by R²
        recommendations.sort(key=lambda x: (x["confidence"], x["r_squared"]), reverse=True)
        
        return recommendations

    def recommend_by_metadata(
        self, metadata_results: pd.DataFrame, metadata: pd.DataFrame
    ) -> List[AssemblyGroup]:
        """Recommend assembly groups based on metadata analysis."""
        logging.info("Generating metadata-based assembly recommendations")

        if metadata_results.empty or metadata is None:
            return []

        groups = []

        # Focus on significant and high-explaining variables
        significant_vars = metadata_results[
            (metadata_results["p_value"] < self.significance_threshold)
            | (metadata_results["r_squared"] > 0.20)  # High explanatory power
        ].sort_values("r_squared", ascending=False)

        for _, row in significant_vars.iterrows():
            variable = row["variable"]

            if variable not in metadata.columns:
                continue

            # Group samples by this metadata variable
            var_data = metadata[variable].dropna()

            if len(var_data) < 2:
                continue

            unique_values = var_data.unique()

            for value in unique_values:
                value_samples = var_data[var_data == value]

                if len(value_samples) < self.min_group_size:
                    continue

                # Get indices for distance calculation
                sample_indices = [
                    self.sample_names.index(name)
                    for name in value_samples.index
                    if name in self.sample_names
                ]

                if len(sample_indices) < self.min_group_size:
                    continue

                avg_dist, max_dist = self._calculate_group_statistics(sample_indices)

                # Calculate confidence using new data-driven method
                group_conf_breakdown = self.confidence_calc.calculate_group_confidence(sample_indices)
                
                # Incorporate statistical significance from metadata analysis
                stat_confidence = min(1.0, row["r_squared"] * 2)
                p_significance = 1.0 if row["p_value"] < 0.01 else 0.9 if row["p_value"] < 0.05 else 0.7
                
                # Combine clustering-based confidence with statistical significance
                confidence = (group_conf_breakdown.total_confidence * 0.7 + 
                            stat_confidence * p_significance * 0.3)

                # Determine benefits and challenges
                benefits = [
                    f"Biologically meaningful grouping by {variable}",
                    "Reduced inter-sample contamination",
                    "Better representation of group-specific features",
                ]

                challenges = [
                    "May miss cross-group shared sequences",
                    f"Groups based on {variable} may have variable quality",
                ]

                if avg_dist > self.similarity_threshold_medium:
                    challenges.append(
                        f"High within-group diversity (avg dist: {avg_dist:.3f})"
                    )

                group = AssemblyGroup(
                    group_id=f"{variable}_{value}",
                    sample_names=list(value_samples.index),
                    grouping_criterion=variable,
                    criterion_value=value,
                    avg_distance=avg_dist,
                    max_distance=max_dist,
                    confidence_score=confidence,
                    expected_benefits=benefits,
                    expected_challenges=challenges,
                )
                groups.append(group)

        return groups

    def recommend_hybrid_strategy(
        self,
        similarity_groups: List[AssemblyGroup],
        metadata_groups: List[AssemblyGroup],
    ) -> List[AssemblyGroup]:
        """Combine similarity and metadata-based recommendations."""
        logging.info("Generating hybrid assembly recommendations")

        # Score and combine approaches
        all_groups = similarity_groups + metadata_groups

        # Groups now have properly calculated confidence scores based on data
        # No need for arbitrary bonuses - the ConfidenceCalculator already considers
        # distance metrics, metadata alignment, and clustering quality

        # Remove overlapping groups, keeping highest confidence
        final_groups = []
        used_samples = set()

        all_groups.sort(key=lambda x: x.confidence_score, reverse=True)

        for group in all_groups:
            group_samples = set(group.sample_names)
            if not group_samples.intersection(used_samples):
                final_groups.append(group)
                used_samples.update(group_samples)

        return final_groups


class AssemblyCommandGenerator:
    """Generate assembly commands for different tools and strategies."""

    def __init__(self, sample_paths: Dict[str, str] = None):
        self.sample_paths = sample_paths or {}

    def generate_megahit_commands(
        self, groups: List[AssemblyGroup]
    ) -> Dict[str, List[str]]:
        """Generate MEGAHIT co-assembly commands."""
        commands = {}

        for group in groups:
            if len(group.sample_names) == 1:
                # Individual assembly
                sample = group.sample_names[0]
                cmd = f"megahit -r {sample}.fastq -o {sample}_assembly --min-contig-len 500"
            else:
                # Co-assembly
                input_files = ",".join(
                    [f"{sample}.fastq" for sample in group.sample_names]
                )
                output_dir = f"{group.group_id}_coassembly"
                cmd = f"megahit -r {input_files} -o {output_dir} --min-contig-len 500 --k-list 21,29,39,59,79,99"

            commands[group.group_id] = [cmd]

        return commands

    def generate_spades_commands(
        self, groups: List[AssemblyGroup]
    ) -> Dict[str, List[str]]:
        """Generate SPAdes metagenomic assembly commands.
        
        Note: SPAdes performs individual assembly on combined reads, not true 
        co-assembly like MEGAHIT. For groups, reads are concatenated before 
        assembly, which may result in suboptimal performance compared to MEGAHIT.
        """
        commands = {}

        for group in groups:
            if len(group.sample_names) == 1:
                # Individual assembly
                sample = group.sample_names[0]
                cmd = f"spades.py --meta -s {sample}.fastq -o {sample}_spades_assembly"
            else:
                # Pseudo-co-assembly: combine reads then individual assembly
                # WARNING: This is not true co-assembly like MEGAHIT
                group_name = group.group_id
                # Cross-platform file combination using Python
                input_files = ' '.join([f'"{sample}.fastq"' for sample in group.sample_names])
                combine_cmd = f'python -c "import shutil; out=open(\'{group_name}_combined.fastq\', \'wb\'); [shutil.copyfileobj(open(f, \'rb\'), out) for f in [{input_files}]]; out.close()"'
                assembly_cmd = f"spades.py --meta -s {group_name}_combined.fastq -o {group_name}_spades_assembly"
                commands[group.group_id] = [combine_cmd, assembly_cmd]
                continue

            commands[group.group_id] = [cmd]

        return commands

    def generate_flye_commands(
        self, groups: List[AssemblyGroup]
    ) -> Dict[str, List[str]]:
        """Generate Flye assembly commands (for long reads)."""
        commands = {}

        for group in groups:
            if len(group.sample_names) == 1:
                sample = group.sample_names[0]
                cmd = f"flye --meta --nano-raw {sample}.fastq -o {sample}_flye_assembly"
            else:
                # Co-assembly
                group_name = group.group_id
                # Cross-platform file combination using Python
                input_files = ' '.join([f'"{sample}.fastq"' for sample in group.sample_names])
                combine_cmd = f'python -c "import shutil; out=open(\'{group_name}_combined.fastq\', \'wb\'); [shutil.copyfileobj(open(f, \'rb\'), out) for f in [{input_files}]]; out.close()"'
                assembly_cmd = f"flye --meta --nano-raw {group_name}_combined.fastq -o {group_name}_flye_assembly"
                commands[group.group_id] = [combine_cmd, assembly_cmd]
                continue

            commands[group.group_id] = [cmd]

        return commands


class PerformancePredictor:
    """Predict assembly performance based on grouping strategy."""

    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str]):
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names

    def predict_assembly_metrics(self, groups: List[AssemblyGroup]) -> Dict[str, Any]:
        """Predict assembly performance metrics."""
        predictions = {}

        total_samples = len(self.sample_names)
        grouped_samples = sum(len(group.sample_names) for group in groups)
        individual_samples = total_samples - grouped_samples

        # Overall strategy metrics
        predictions["strategy_summary"] = {
            "total_assemblies": len(groups) + individual_samples,
            "co_assemblies": len(groups),
            "individual_assemblies": individual_samples,
            "samples_in_coassembly": grouped_samples,
            "coassembly_percentage": (grouped_samples / total_samples) * 100,
        }

        # Per-group predictions
        group_predictions = {}
        for group in groups:
            n_samples = len(group.sample_names)
            avg_dist = group.avg_distance

            # Predict relative benefits/challenges
            expected_contiguity_improvement = max(
                1.0, 2.0 - avg_dist * 3
            )  # Lower distance = better contiguity
            expected_coverage_boost = min(n_samples * 0.8, 5.0)  # Diminishing returns
            computational_cost_multiplier = n_samples**1.5  # Non-linear scaling

            # Assembly quality predictions
            if avg_dist < 0.1:
                predicted_quality = "Very High"
                contamination_risk = "Low"
            elif avg_dist < 0.2:
                predicted_quality = "High"
                contamination_risk = "Low-Medium"
            elif avg_dist < 0.3:
                predicted_quality = "Medium"
                contamination_risk = "Medium"
            else:
                predicted_quality = "Low-Medium"
                contamination_risk = "High"

            group_predictions[group.group_id] = {
                "n_samples": n_samples,
                "avg_distance": avg_dist,
                "predicted_quality": predicted_quality,
                "contamination_risk": contamination_risk,
                "expected_contiguity_improvement": f"{expected_contiguity_improvement:.1f}x",
                "expected_coverage_boost": f"{expected_coverage_boost:.1f}x",
                "computational_cost_multiplier": f"{computational_cost_multiplier:.1f}x",
                "confidence": group.confidence_score,
            }

        predictions["group_predictions"] = group_predictions

        return predictions


class AssemblyRecommender:
    """Main class for generating comprehensive assembly recommendations."""

    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str], metadata: Optional[pd.DataFrame] = None):
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names
        self.metadata = metadata
        self.strategy_engine = AssemblyStrategyEngine(distance_matrix, sample_names, metadata)
        self.command_generator = AssemblyCommandGenerator()
        self.performance_predictor = PerformancePredictor(distance_matrix, sample_names)

    def generate_metadata_grouping_recommendations(
        self, metadata_results: pd.DataFrame, metadata: pd.DataFrame
    ) -> List[Dict[str, Any]]:
        """
        Generate intelligent 'group by X' recommendations based on metadata analysis.
        
        Args:
            metadata_results: PERMANOVA results DataFrame
            metadata: Metadata DataFrame
            
        Returns:
            List of grouping recommendations with confidence scores and explanations
        """
        return self.strategy_engine.generate_metadata_grouping_recommendations(
            metadata_results, metadata
        )

    def generate_recommendations(
        self,
        metadata_results: Optional[pd.DataFrame] = None,
        metadata: Optional[pd.DataFrame] = None,
    ) -> AssemblyRecommendation:
        """Generate comprehensive assembly recommendations."""
        logging.info("Generating comprehensive assembly recommendations")

        # Get recommendations from different approaches
        similarity_groups = self.strategy_engine.recommend_by_similarity()

        metadata_groups = []
        if metadata_results is not None and not metadata_results.empty:
            metadata_groups = self.strategy_engine.recommend_by_metadata(
                metadata_results, metadata
            )

        # Combine approaches
        if metadata_groups:
            final_groups = self.strategy_engine.recommend_hybrid_strategy(
                similarity_groups, metadata_groups
            )
        else:
            final_groups = similarity_groups

        # Determine overall strategy
        total_samples = len(self.sample_names)
        grouped_samples = sum(len(group.sample_names) for group in final_groups)

        # Calculate strategy confidence and store breakdown for reporting
        confidence_breakdown = None
        
        if not final_groups or grouped_samples < 2:
            strategy = "individual"
            rationale = (
                "No clear grouping patterns found. Individual assembly recommended."
            )
            # Calculate individual assembly confidence based on data
            confidence_breakdown = self.strategy_engine.confidence_calc.calculate_overall_strategy_confidence("individual", [])
            overall_confidence = confidence_breakdown.total_confidence
        elif grouped_samples == total_samples and len(final_groups) == 1:
            strategy = "global"
            rationale = (
                "All samples show strong similarity. Global co-assembly recommended."
            )
            # Calculate global strategy confidence using new method
            confidence_breakdown = self.strategy_engine.confidence_calc.calculate_overall_strategy_confidence("global", [])
            overall_confidence = confidence_breakdown.total_confidence
        else:
            strategy = "grouped"
            rationale = f"Mixed strategy: {len(final_groups)} co-assembly groups covering {grouped_samples}/{total_samples} samples."
            # Calculate grouped strategy confidence using new method
            group_indices_list = []
            for group in final_groups:
                group_sample_indices = [self.sample_names.index(name) for name in group.sample_names if name in self.sample_names]
                if group_sample_indices:
                    group_indices_list.append(group_sample_indices)
            confidence_breakdown = self.strategy_engine.confidence_calc.calculate_overall_strategy_confidence("grouped", group_indices_list)
            overall_confidence = confidence_breakdown.total_confidence

        # Determine primary criterion
        primary_criterion = None
        if final_groups:
            # Find most common grouping criterion
            criteria = [group.grouping_criterion for group in final_groups]
            primary_criterion = max(set(criteria), key=criteria.count)

        # Generate assembly commands
        assembly_commands = {
            "megahit": self.command_generator.generate_megahit_commands(final_groups),
            "spades": self.command_generator.generate_spades_commands(final_groups),
            "flye": self.command_generator.generate_flye_commands(final_groups),
        }

        # Predict performance
        performance_predictions = self.performance_predictor.predict_assembly_metrics(
            final_groups
        )

        recommendation = AssemblyRecommendation(
            strategy=strategy,
            groups=final_groups,
            overall_confidence=overall_confidence,
            confidence_breakdown=confidence_breakdown,
            primary_criterion=primary_criterion,
            decision_rationale=rationale,
            assembly_commands=assembly_commands,
            performance_predictions=performance_predictions,
        )

        return recommendation


def save_recommendations(recommendation: AssemblyRecommendation, output_path: str):
    """Save assembly recommendations to files."""
    output_dir = Path(output_path)
    output_dir.mkdir(exist_ok=True)

    # Save detailed recommendation as JSON
    rec_dict = asdict(recommendation)
    with open(output_dir / "assembly_recommendations.json", "w") as f:
        json.dump(rec_dict, f, indent=2, default=str)

    # Save human-readable summary
    with open(output_dir / "assembly_strategy.md", "w") as f:
        f.write("# MetaGrouper Assembly Strategy Recommendations\n\n")

        f.write(f"## Overall Strategy: {recommendation.strategy.title()}\n\n")
        f.write(f"**Confidence Score:** {recommendation.overall_confidence:.2f}\n\n")
        f.write(f"**Rationale:** {recommendation.decision_rationale}\n\n")

        if recommendation.primary_criterion:
            f.write(
                f"**Primary Grouping Criterion:** {recommendation.primary_criterion}\n\n"
            )

        # Strategy summary
        perf = recommendation.performance_predictions["strategy_summary"]
        f.write("## Strategy Summary\n\n")
        f.write(f"- **Total Assemblies:** {perf['total_assemblies']}\n")
        f.write(f"- **Co-assemblies:** {perf['co_assemblies']}\n")
        f.write(f"- **Individual Assemblies:** {perf['individual_assemblies']}\n")
        f.write(
            f"- **Samples in Co-assembly:** {perf['samples_in_coassembly']} ({perf['coassembly_percentage']:.1f}%)\n\n"
        )

        # Group details
        if recommendation.groups:
            f.write("## Recommended Assembly Groups\n\n")
            for i, group in enumerate(recommendation.groups, 1):
                f.write(f"### Group {i}: {group.group_id}\n\n")
                f.write(f"- **Samples:** {', '.join(group.sample_names)}\n")
                f.write(f"- **Grouping Criterion:** {group.grouping_criterion}\n")
                f.write(f"- **Criterion Value:** {group.criterion_value}\n")
                f.write(f"- **Average Distance:** {group.avg_distance:.3f}\n")
                f.write(f"- **Confidence Score:** {group.confidence_score:.2f}\n\n")


        # Assembly commands
        f.write("## Assembly Commands\n\n")
        for tool, commands in recommendation.assembly_commands.items():
            f.write(f"### {tool.upper()}\n\n")
            for group_id, cmd_list in commands.items():
                f.write(f"**{group_id}:**\n")
                for cmd in cmd_list:
                    f.write(f"```bash\n{cmd}\n```\n\n")

    # Save assembly commands as shell scripts
    for tool, commands in recommendation.assembly_commands.items():
        script_content = "#!/bin/bash\n\n"
        script_content += f"# MetaGrouper {tool.upper()} Assembly Commands\n"
        script_content += (
            f"# Generated assembly strategy: {recommendation.strategy}\n\n"
        )

        for group_id, cmd_list in commands.items():
            script_content += f"# Group: {group_id}\n"
            for cmd in cmd_list:
                script_content += f"{cmd}\n"
            script_content += "\n"

        with open(output_dir / f"run_{tool}_assemblies.sh", "w") as f:
            f.write(script_content)

    logging.info(f"Assembly recommendations saved to {output_path}")


def visualize_assembly_strategy(
    recommendation: AssemblyRecommendation,
    distance_matrix: np.ndarray,
    sample_names: List[str],
    output_path: str,
):
    """Create visualization of the recommended assembly strategy."""

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # 1. Distance matrix with group overlays
    im1 = ax1.imshow(distance_matrix, cmap="viridis", aspect="auto")
    ax1.set_title("Sample Distance Matrix with Assembly Groups")
    ax1.set_xticks(range(len(sample_names)))
    ax1.set_yticks(range(len(sample_names)))
    ax1.set_xticklabels(sample_names, rotation=45, ha="right")
    ax1.set_yticklabels(sample_names)

    # Overlay group boundaries
    colors = plt.cm.Set1(np.linspace(0, 1, len(recommendation.groups)))
    used_samples = set()

    for i, (group, color) in enumerate(zip(recommendation.groups, colors)):
        group_indices = [
            sample_names.index(name)
            for name in group.sample_names
            if name not in used_samples
        ]
        if len(group_indices) > 1:
            for idx in group_indices:
                ax1.axhline(y=idx - 0.5, color=color, linewidth=2, alpha=0.7)
                ax1.axvline(x=idx - 0.5, color=color, linewidth=2, alpha=0.7)
        used_samples.update(group.sample_names)

    plt.colorbar(im1, ax=ax1, label="Distance")

    # 2. Group confidence scores
    if recommendation.groups:
        group_names = [f"{group.group_id}" for group in recommendation.groups]
        confidences = [group.confidence_score for group in recommendation.groups]

        bars = ax2.bar(
            range(len(group_names)),
            confidences,
            color=colors[: len(group_names)],
            alpha=0.7,
        )
        ax2.set_title("Assembly Group Confidence Scores")
        ax2.set_xlabel("Assembly Groups")
        ax2.set_ylabel("Confidence Score")
        ax2.set_xticks(range(len(group_names)))
        ax2.set_xticklabels(group_names, rotation=45, ha="right")
        ax2.set_ylim(0, 1)
        ax2.grid(axis="y", alpha=0.3)

        # Add confidence values on bars
        for bar, conf in zip(bars, confidences):
            height = bar.get_height()
            ax2.text(
                bar.get_x() + bar.get_width() / 2.0,
                height + 0.01,
                f"{conf:.2f}",
                ha="center",
                va="bottom",
            )
    else:
        ax2.text(
            0.5,
            0.5,
            "No groups recommended\nIndividual assembly suggested",
            ha="center",
            va="center",
            transform=ax2.transAxes,
        )
        ax2.set_title("Assembly Group Confidence Scores")

    # 3. Strategy overview pie chart
    perf = recommendation.performance_predictions["strategy_summary"]

    labels = []
    sizes = []
    colors_pie = []

    if perf["co_assemblies"] > 0:
        labels.append(f"Co-assemblies ({perf['co_assemblies']})")
        sizes.append(perf["co_assemblies"])
        colors_pie.append("lightblue")

    if perf["individual_assemblies"] > 0:
        labels.append(f"Individual ({perf['individual_assemblies']})")
        sizes.append(perf["individual_assemblies"])
        colors_pie.append("lightcoral")

    if sizes:
        ax3.pie(
            sizes, labels=labels, colors=colors_pie, autopct="%1.1f%%", startangle=90
        )
    ax3.set_title(
        f"Assembly Strategy Distribution\n({recommendation.strategy.title()})"
    )

    # 4. Performance predictions
    if recommendation.groups:
        group_sizes = [len(group.sample_names) for group in recommendation.groups]
        avg_distances = [group.avg_distance for group in recommendation.groups]

        scatter = ax4.scatter(
            group_sizes,
            avg_distances,
            c=confidences,
            cmap="RdYlGn",
            s=100,
            alpha=0.7,
            edgecolors="black",
        )

        ax4.set_xlabel("Group Size (number of samples)")
        ax4.set_ylabel("Average Intra-group Distance")
        ax4.set_title("Group Size vs. Distance (colored by confidence)")
        ax4.grid(True, alpha=0.3)

        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax4)
        cbar.set_label("Confidence Score")

        # Annotate points
        for i, group in enumerate(recommendation.groups):
            ax4.annotate(
                f"G{i+1}",
                (group_sizes[i], avg_distances[i]),
                xytext=(5, 5),
                textcoords="offset points",
            )
    else:
        ax4.text(
            0.5,
            0.5,
            "No groups to analyze",
            ha="center",
            va="center",
            transform=ax4.transAxes,
        )
        ax4.set_title("Group Size vs. Distance Analysis")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    logging.info(f"Assembly strategy visualization saved to {output_path}")
