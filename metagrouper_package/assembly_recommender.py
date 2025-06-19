#!/usr/bin/env python3
"""
Phase 3: Assembly Strategy Recommendations for MetaGrouper

Generates optimal assembly strategies based on k-mer similarity patterns
and metadata associations, with support for multiple assemblers.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import logging
from dataclasses import dataclass
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import json


@dataclass
class AssemblyGroup:
    """Represents a group of samples for co-assembly."""
    group_id: int
    sample_names: List[str]
    strategy: str  # 'individual', 'grouped', 'global'
    confidence: float
    similarity_stats: Dict[str, float]
    metadata_support: Optional[Dict[str, Any]] = None


@dataclass
class AssemblyRecommendation:
    """Complete assembly strategy recommendation."""
    strategy: str  # Overall strategy: 'individual', 'grouped', 'global'
    groups: List[AssemblyGroup]
    overall_confidence: float
    summary: str
    assembly_commands: Dict[str, List[str]]  # Tool -> list of commands
    rationale: List[str]


class AssemblyStrategyEngine:
    """Core engine for determining optimal assembly strategies."""
    
    def __init__(self):
        # Strategy thresholds
        self.similarity_threshold_high = 0.85  # Very similar -> global assembly
        self.similarity_threshold_medium = 0.30  # Moderately similar -> grouped assembly
        self.min_group_size = 2
        self.max_group_size = 10
        
        # Confidence factors
        self.confidence_weights = {
            'similarity_consistency': 0.4,
            'group_size_optimality': 0.2,
            'metadata_support': 0.3,
            'clustering_quality': 0.1
        }
    
    def analyze_similarity_structure(self, distance_matrix: np.ndarray) -> Dict[str, Any]:
        """Analyze the overall structure of sample similarities."""
        # Convert distance to similarity
        similarity_matrix = 1 - distance_matrix
        
        # Basic statistics
        upper_tri_mask = np.triu(np.ones_like(similarity_matrix, dtype=bool), k=1)
        similarities = similarity_matrix[upper_tri_mask]
        
        stats = {
            'mean_similarity': np.mean(similarities),
            'median_similarity': np.median(similarities),
            'std_similarity': np.std(similarities),
            'min_similarity': np.min(similarities),
            'max_similarity': np.max(similarities),
            'q25_similarity': np.percentile(similarities, 25),
            'q75_similarity': np.percentile(similarities, 75)
        }
        
        # Determine overall similarity pattern
        if stats['median_similarity'] > self.similarity_threshold_high:
            pattern = 'very_high'
        elif stats['median_similarity'] > self.similarity_threshold_medium:
            pattern = 'moderate'
        elif stats['q75_similarity'] > self.similarity_threshold_medium:
            pattern = 'mixed'
        else:
            pattern = 'low'
        
        stats['similarity_pattern'] = pattern
        
        return stats
    
    def identify_assembly_groups(self, distance_matrix: np.ndarray, sample_names: List[str],
                               metadata_results: Optional[pd.DataFrame] = None) -> List[AssemblyGroup]:
        """Identify optimal groupings for assembly."""
        n_samples = len(sample_names)
        
        # Analyze overall similarity structure
        similarity_stats = self.analyze_similarity_structure(distance_matrix)
        
        groups = []
        
        # Strategy 1: Check if all samples are very similar (global assembly)
        if similarity_stats['similarity_pattern'] == 'very_high':
            group = AssemblyGroup(
                group_id=0,
                sample_names=sample_names,
                strategy='global',
                confidence=min(0.9, similarity_stats['mean_similarity']),
                similarity_stats=similarity_stats
            )
            groups.append(group)
            return groups
        
        # Strategy 2: Check if samples are too diverse (individual assembly)
        if similarity_stats['similarity_pattern'] == 'low':
            for i, sample_name in enumerate(sample_names):
                group = AssemblyGroup(
                    group_id=i,
                    sample_names=[sample_name],
                    strategy='individual',
                    confidence=0.7,  # Conservative confidence for individual assembly
                    similarity_stats={'mean_similarity': 0.0}  # No co-assembly
                )
                groups.append(group)
            return groups
        
        # Strategy 3: Hierarchical clustering for grouped assembly
        try:
            # Create hierarchical clustering
            condensed_distances = squareform(distance_matrix, checks=False)
            linkage_matrix = linkage(condensed_distances, method='ward')
            
            # Find optimal number of clusters
            best_n_clusters = self._find_optimal_clusters(distance_matrix, linkage_matrix)
            
            # Get cluster assignments
            cluster_labels = fcluster(linkage_matrix, best_n_clusters, criterion='maxclust') - 1
            
            # Create groups
            for cluster_id in range(best_n_clusters):
                cluster_mask = cluster_labels == cluster_id
                cluster_samples = [sample_names[i] for i in range(len(sample_names)) if cluster_mask[i]]
                
                if len(cluster_samples) == 0:
                    continue
                
                # Calculate within-cluster similarity statistics
                cluster_indices = np.where(cluster_mask)[0]
                if len(cluster_indices) > 1:
                    cluster_distances = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
                    cluster_similarities = 1 - cluster_distances
                    
                    # Remove diagonal
                    cluster_similarities_flat = cluster_similarities[np.triu(np.ones_like(cluster_similarities, dtype=bool), k=1)]
                    
                    cluster_stats = {
                        'mean_similarity': np.mean(cluster_similarities_flat),
                        'min_similarity': np.min(cluster_similarities_flat),
                        'std_similarity': np.std(cluster_similarities_flat)
                    }
                else:
                    cluster_stats = {'mean_similarity': 1.0, 'min_similarity': 1.0, 'std_similarity': 0.0}
                
                # Determine strategy for this group
                if len(cluster_samples) == 1:
                    strategy = 'individual'
                    confidence = 0.7
                elif len(cluster_samples) > self.max_group_size:
                    # Large group - check if should be global or split further
                    if cluster_stats['mean_similarity'] > self.similarity_threshold_high:
                        strategy = 'global'
                        confidence = cluster_stats['mean_similarity'] * 0.9
                    else:
                        strategy = 'grouped'  # Will be split in command generation
                        confidence = cluster_stats['mean_similarity'] * 0.8
                else:
                    strategy = 'grouped'
                    confidence = cluster_stats['mean_similarity'] * 0.9
                
                group = AssemblyGroup(
                    group_id=cluster_id,
                    sample_names=cluster_samples,
                    strategy=strategy,
                    confidence=confidence,
                    similarity_stats=cluster_stats
                )
                
                groups.append(group)
            
        except Exception as e:
            logging.warning(f"Clustering failed, falling back to individual assembly: {e}")
            # Fallback to individual assembly
            for i, sample_name in enumerate(sample_names):
                group = AssemblyGroup(
                    group_id=i,
                    sample_names=[sample_name],
                    strategy='individual',
                    confidence=0.6,  # Lower confidence due to fallback
                    similarity_stats={'mean_similarity': 0.0}
                )
                groups.append(group)
        
        return groups
    
    def _find_optimal_clusters(self, distance_matrix: np.ndarray, linkage_matrix: np.ndarray) -> int:
        """Find optimal number of clusters using silhouette analysis."""
        n_samples = distance_matrix.shape[0]
        max_clusters = min(8, n_samples - 1)
        
        if max_clusters <= 2:
            return min(2, n_samples)
        
        best_score = -1
        best_n_clusters = 2
        
        try:
            from sklearn.metrics import silhouette_score
            
            for n_clusters in range(2, max_clusters + 1):
                cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust') - 1
                
                # Calculate silhouette score
                score = silhouette_score(distance_matrix, cluster_labels, metric='precomputed')
                
                if score > best_score:
                    best_score = score
                    best_n_clusters = n_clusters
                    
        except ImportError:
            # Fallback: use elbow method with within-cluster sum of squares
            wcss_scores = []
            for n_clusters in range(2, max_clusters + 1):
                cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust') - 1
                wcss = 0
                
                for cluster_id in range(n_clusters):
                    cluster_indices = np.where(cluster_labels == cluster_id)[0]
                    if len(cluster_indices) > 1:
                        cluster_distances = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
                        wcss += np.sum(cluster_distances ** 2)
                
                wcss_scores.append(wcss)
            
            # Find elbow point (simple heuristic)
            if len(wcss_scores) > 2:
                diffs = np.diff(wcss_scores)
                second_diffs = np.diff(diffs)
                elbow_idx = np.argmax(second_diffs) + 2  # +2 because we start from 2 clusters
                best_n_clusters = min(elbow_idx + 2, max_clusters)
        
        return best_n_clusters
    
    def determine_overall_strategy(self, groups: List[AssemblyGroup]) -> str:
        """Determine overall assembly strategy."""
        if len(groups) == 1:
            return groups[0].strategy
        
        # Count strategies
        strategy_counts = {}
        for group in groups:
            strategy_counts[group.strategy] = strategy_counts.get(group.strategy, 0) + 1
        
        # If mostly individual assembly
        if strategy_counts.get('individual', 0) > len(groups) * 0.7:
            return 'individual'
        
        # If any global assembly
        if strategy_counts.get('global', 0) > 0:
            return 'mixed'  # Mix of global and other strategies
        
        # Otherwise grouped
        return 'grouped'
    
    def calculate_overall_confidence(self, groups: List[AssemblyGroup]) -> float:
        """Calculate overall confidence in the strategy."""
        if not groups:
            return 0.0
        
        # Weight by group size
        total_samples = sum(len(group.sample_names) for group in groups)
        weighted_confidence = sum(
            group.confidence * len(group.sample_names) / total_samples
            for group in groups
        )
        
        # Penalty for too many small groups (fragmentation)
        small_groups = sum(1 for group in groups if len(group.sample_names) == 1)
        fragmentation_penalty = min(0.2, small_groups / len(groups) * 0.3)
        
        return max(0.1, weighted_confidence - fragmentation_penalty)


class AssemblyCommandGenerator:
    """Generates assembly commands for different tools."""
    
    def __init__(self):
        # Default parameters for different assemblers
        self.megahit_params = {
            'individual': '--min-contig-len 1000 --k-min 21 --k-max 141 --k-step 12',
            'grouped': '--min-contig-len 1000 --k-min 21 --k-max 141 --k-step 12 --merge-level 20,0.95',
            'global': '--min-contig-len 1000 --k-min 21 --k-max 141 --k-step 12 --merge-level 20,0.98'
        }
        
        self.spades_params = {
            'individual': '--meta -k 21,33,55,77',
            'grouped': '--meta -k 21,33,55,77,99',
            'global': '--meta -k 21,33,55,77,99,121'
        }
        
        self.flye_params = {
            'individual': '--meta --genome-size 5m',
            'grouped': '--meta --genome-size 10m',
            'global': '--meta --genome-size 50m'
        }
    
    def generate_megahit_commands(self, groups: List[AssemblyGroup], 
                                 sample_files: Dict[str, str]) -> List[str]:
        """Generate MEGAHIT assembly commands."""
        commands = []
        max_group_size = 10  # Default max group size
        
        for group in groups:
            if group.strategy == 'individual':
                for sample in group.sample_names:
                    if sample in sample_files:
                        cmd = (f"megahit -1 {sample_files[sample]} "
                              f"{self.megahit_params['individual']} "
                              f"-o {sample}_megahit --out-prefix {sample}")
                        commands.append(cmd)
            else:
                # Grouped or global assembly
                input_files = [sample_files[sample] for sample in group.sample_names if sample in sample_files]
                if input_files:
                    if len(group.sample_names) > max_group_size:
                        # Split large groups
                        for i in range(0, len(input_files), max_group_size):
                            batch_files = input_files[i:i+max_group_size]
                            batch_samples = group.sample_names[i:i+max_group_size]
                            group_name = f"group_{group.group_id}_batch_{i//max_group_size}"
                            cmd = (f"megahit -1 {','.join(batch_files)} "
                                  f"{self.megahit_params[group.strategy]} "
                                  f"-o {group_name}_megahit --out-prefix {group_name}")
                            commands.append(cmd)
                    else:
                        group_name = f"group_{group.group_id}"
                        cmd = (f"megahit -1 {','.join(input_files)} "
                              f"{self.megahit_params[group.strategy]} "
                              f"-o {group_name}_megahit --out-prefix {group_name}")
                        commands.append(cmd)
        
        return commands
    
    def generate_spades_commands(self, groups: List[AssemblyGroup],
                               sample_files: Dict[str, str]) -> List[str]:
        """Generate SPAdes assembly commands."""
        commands = []
        
        for group in groups:
            if group.strategy == 'individual':
                for sample in group.sample_names:
                    if sample in sample_files:
                        cmd = (f"spades.py {self.spades_params['individual']} "
                              f"-1 {sample_files[sample]} "
                              f"-o {sample}_spades")
                        commands.append(cmd)
            else:
                # Grouped or global assembly - SPAdes doesn't support direct multi-sample input
                # Generate individual commands but note they're part of a group
                for sample in group.sample_names:
                    if sample in sample_files:
                        cmd = (f"spades.py {self.spades_params[group.strategy]} "
                              f"-1 {sample_files[sample]} "
                              f"-o {sample}_spades_group_{group.group_id}")
                        commands.append(cmd)
        
        return commands
    
    def generate_flye_commands(self, groups: List[AssemblyGroup],
                             sample_files: Dict[str, str]) -> List[str]:
        """Generate Flye assembly commands (for long reads)."""
        commands = []
        
        for group in groups:
            if group.strategy == 'individual':
                for sample in group.sample_names:
                    if sample in sample_files:
                        cmd = (f"flye --nano-raw {sample_files[sample]} "
                              f"{self.flye_params['individual']} "
                              f"--out-dir {sample}_flye")
                        commands.append(cmd)
            else:
                # Combine reads for grouped assembly
                input_files = [sample_files[sample] for sample in group.sample_names if sample in sample_files]
                if input_files:
                    group_name = f"group_{group.group_id}"
                    # Note: This assumes reads can be concatenated - may need adjustment
                    cmd = (f"cat {' '.join(input_files)} > {group_name}_combined.fastq && "
                          f"flye --nano-raw {group_name}_combined.fastq "
                          f"{self.flye_params[group.strategy]} "
                          f"--out-dir {group_name}_flye")
                    commands.append(cmd)
        
        return commands


class AssemblyRecommender:
    """Main class for generating assembly strategy recommendations."""
    
    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str]):
        """
        Initialize assembly recommender.
        
        Args:
            distance_matrix: Pairwise distance matrix between samples
            sample_names: List of sample names
        """
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names
        self.strategy_engine = AssemblyStrategyEngine()
        self.command_generator = AssemblyCommandGenerator()
    
    def generate_recommendations(self, metadata_results: Optional[pd.DataFrame] = None,
                               metadata: Optional[pd.DataFrame] = None) -> AssemblyRecommendation:
        """
        Generate complete assembly strategy recommendation.
        
        Args:
            metadata_results: Results from metadata analysis
            metadata: Original metadata DataFrame
            
        Returns:
            AssemblyRecommendation object
        """
        # Identify assembly groups
        groups = self.strategy_engine.identify_assembly_groups(
            self.distance_matrix, self.sample_names, metadata_results
        )
        
        # Enhance groups with metadata information
        if metadata_results is not None and metadata is not None:
            groups = self._enhance_groups_with_metadata(groups, metadata_results, metadata)
        
        # Determine overall strategy
        overall_strategy = self.strategy_engine.determine_overall_strategy(groups)
        
        # Calculate overall confidence
        overall_confidence = self.strategy_engine.calculate_overall_confidence(groups)
        
        # Generate rationale
        rationale = self._generate_rationale(groups, overall_strategy, metadata_results)
        
        # Generate summary
        summary = self._generate_summary(groups, overall_strategy, overall_confidence)
        
        # Generate assembly commands (placeholder file paths)
        sample_files = {sample: f"{sample}.fastq.gz" for sample in self.sample_names}
        
        assembly_commands = {
            'megahit': self.command_generator.generate_megahit_commands(groups, sample_files),
            'spades': self.command_generator.generate_spades_commands(groups, sample_files),
            'flye': self.command_generator.generate_flye_commands(groups, sample_files)
        }
        
        return AssemblyRecommendation(
            strategy=overall_strategy,
            groups=groups,
            overall_confidence=overall_confidence,
            summary=summary,
            assembly_commands=assembly_commands,
            rationale=rationale
        )
    
    def _enhance_groups_with_metadata(self, groups: List[AssemblyGroup],
                                    metadata_results: pd.DataFrame,
                                    metadata: pd.DataFrame) -> List[AssemblyGroup]:
        """Enhance groups with metadata information."""
        # Find most significant metadata variables
        significant_vars = metadata_results[metadata_results['p_value'] < 0.05]
        
        for group in groups:
            if len(group.sample_names) > 1:
                # Check metadata consistency within group
                group_metadata = metadata.loc[group.sample_names]
                
                metadata_support = {}
                for _, var_result in significant_vars.iterrows():
                    var_name = var_result['variable']
                    if var_name in group_metadata.columns:
                        var_values = group_metadata[var_name].dropna()
                        if len(var_values) > 0:
                            unique_values = var_values.nunique()
                            consistency = 1 - (unique_values - 1) / len(var_values)
                            metadata_support[var_name] = {
                                'consistency': consistency,
                                'unique_values': unique_values,
                                'total_samples': len(var_values)
                            }
                
                group.metadata_support = metadata_support
                
                # Adjust confidence based on metadata consistency
                if metadata_support:
                    avg_consistency = np.mean([info['consistency'] for info in metadata_support.values()])
                    group.confidence = min(0.95, group.confidence * (0.7 + 0.3 * avg_consistency))
        
        return groups
    
    def _generate_rationale(self, groups: List[AssemblyGroup], overall_strategy: str,
                           metadata_results: Optional[pd.DataFrame]) -> List[str]:
        """Generate rationale for the recommended strategy."""
        rationale = []
        
        # Overall strategy rationale
        if overall_strategy == 'individual':
            rationale.append("Individual assembly recommended due to high sample diversity")
        elif overall_strategy == 'global':
            rationale.append("Global assembly recommended due to high sample similarity")
        elif overall_strategy == 'grouped':
            rationale.append("Grouped assembly recommended based on similarity clustering")
        else:
            rationale.append("Mixed strategy recommended based on heterogeneous sample patterns")
        
        # Group-specific rationale
        grouped_count = sum(1 for g in groups if g.strategy == 'grouped' and len(g.sample_names) > 1)
        if grouped_count > 0:
            rationale.append(f"Identified {grouped_count} groups for co-assembly")
        
        individual_count = sum(1 for g in groups if len(g.sample_names) == 1)
        if individual_count > 0:
            rationale.append(f"{individual_count} samples recommended for individual assembly")
        
        # Metadata rationale
        if metadata_results is not None and not metadata_results.empty:
            significant_vars = metadata_results[metadata_results['p_value'] < 0.05]
            if len(significant_vars) > 0:
                top_var = significant_vars.iloc[0]
                rationale.append(f"Strategy supported by metadata variable '{top_var['variable']}' "
                               f"(R² = {top_var['r_squared']:.3f})")
        
        return rationale
    
    def _generate_summary(self, groups: List[AssemblyGroup], overall_strategy: str,
                         overall_confidence: float) -> str:
        """Generate human-readable summary."""
        n_groups = len([g for g in groups if len(g.sample_names) > 1])
        n_individual = len([g for g in groups if len(g.sample_names) == 1])
        
        summary_parts = [
            f"Recommended strategy: {overall_strategy.title()} Assembly",
            f"Confidence: {overall_confidence:.1%}",
        ]
        
        if n_groups > 0:
            summary_parts.append(f"Co-assembly groups: {n_groups}")
        
        if n_individual > 0:
            summary_parts.append(f"Individual assemblies: {n_individual}")
        
        return " | ".join(summary_parts)


def save_recommendations(recommendation: AssemblyRecommendation, output_dir: str):
    """Save assembly recommendations to files."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save JSON results
    json_data = {
        'strategy': recommendation.strategy,
        'overall_confidence': recommendation.overall_confidence,
        'summary': recommendation.summary,
        'rationale': recommendation.rationale,
        'groups': [
            {
                'group_id': group.group_id,
                'sample_names': group.sample_names,
                'strategy': group.strategy,
                'confidence': group.confidence,
                'similarity_stats': group.similarity_stats,
                'metadata_support': group.metadata_support
            }
            for group in recommendation.groups
        ]
    }
    
    with open(output_path / 'assembly_recommendations.json', 'w') as f:
        json.dump(json_data, f, indent=2)
    
    # Save human-readable strategy
    strategy_lines = [
        "# MetaGrouper Assembly Strategy Recommendation",
        "",
        f"**Strategy:** {recommendation.strategy.title()} Assembly",
        f"**Confidence:** {recommendation.overall_confidence:.1%}",
        "",
        "## Summary",
        recommendation.summary,
        "",
        "## Rationale",
    ]
    
    for reason in recommendation.rationale:
        strategy_lines.append(f"- {reason}")
    
    strategy_lines.extend([
        "",
        "## Groups",
        ""
    ])
    
    for group in recommendation.groups:
        strategy_lines.extend([
            f"### Group {group.group_id} ({group.strategy.title()})",
            f"- **Samples:** {len(group.sample_names)}",
            f"- **Confidence:** {group.confidence:.1%}",
            f"- **Mean similarity:** {group.similarity_stats.get('mean_similarity', 0):.3f}",
            f"- **Sample list:** {', '.join(group.sample_names[:10])}{'...' if len(group.sample_names) > 10 else ''}",
            ""
        ])
    
    with open(output_path / 'assembly_strategy.md', 'w') as f:
        f.write('\n'.join(strategy_lines))
    
    # Save assembly commands
    for tool, commands in recommendation.assembly_commands.items():
        if commands:
            script_path = output_path / f'run_{tool}_assemblies.sh'
            with open(script_path, 'w') as f:
                f.write(f"#!/bin/bash\n")
                f.write(f"# {tool.upper()} assembly commands generated by MetaGrouper\n")
                f.write(f"# Strategy: {recommendation.strategy}\n")
                f.write(f"# Confidence: {recommendation.overall_confidence:.1%}\n\n")
                
                for i, cmd in enumerate(commands):
                    f.write(f"# Command {i+1}\n")
                    f.write(f"{cmd}\n\n")
            
            # Make script executable
            script_path.chmod(0o755)
    
    logging.info(f"Assembly recommendations saved to {output_path}")


def visualize_assembly_strategy(recommendation: AssemblyRecommendation,
                              distance_matrix: np.ndarray,
                              sample_names: List[str],
                              output_path: str):
    """Create visualization of assembly strategy."""
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. Distance heatmap with group annotations
        group_colors = plt.cm.Set1(np.linspace(0, 1, len(recommendation.groups)))
        sample_to_group = {}
        sample_to_color = {}
        
        for i, group in enumerate(recommendation.groups):
            for sample in group.sample_names:
                sample_to_group[sample] = i
                sample_to_color[sample] = group_colors[i]
        
        # Reorder samples by group
        ordered_samples = []
        for group in recommendation.groups:
            ordered_samples.extend(group.sample_names)
        
        sample_indices = [sample_names.index(sample) for sample in ordered_samples]
        ordered_matrix = distance_matrix[np.ix_(sample_indices, sample_indices)]
        
        im1 = ax1.imshow(ordered_matrix, cmap='viridis', aspect='auto')
        ax1.set_title('Sample Distance Matrix\n(Ordered by Groups)')
        ax1.set_xticks(range(len(ordered_samples)))
        ax1.set_yticks(range(len(ordered_samples)))
        ax1.set_xticklabels(ordered_samples, rotation=90, fontsize=8)
        ax1.set_yticklabels(ordered_samples, fontsize=8)
        plt.colorbar(im1, ax=ax1, label='Distance')
        
        # 2. Group size distribution
        group_sizes = [len(group.sample_names) for group in recommendation.groups]
        ax2.hist(group_sizes, bins=max(1, len(set(group_sizes))), alpha=0.7, edgecolor='black')
        ax2.set_xlabel('Group Size')
        ax2.set_ylabel('Number of Groups')
        ax2.set_title('Assembly Group Size Distribution')
        
        # 3. Confidence scores by group
        confidences = [group.confidence for group in recommendation.groups]
        group_ids = [f"Group {group.group_id}" for group in recommendation.groups]
        
        bars = ax3.bar(group_ids, confidences, color=[sample_to_color[group.sample_names[0]] for group in recommendation.groups])
        ax3.set_ylabel('Confidence Score')
        ax3.set_title('Group Confidence Scores')
        ax3.set_ylim(0, 1)
        ax3.tick_params(axis='x', rotation=45)
        
        # Add horizontal line for overall confidence
        ax3.axhline(y=recommendation.overall_confidence, color='red', linestyle='--', 
                   label=f'Overall: {recommendation.overall_confidence:.2f}')
        ax3.legend()
        
        # 4. Strategy summary
        ax4.axis('off')
        summary_text = [
            f"Assembly Strategy: {recommendation.strategy.title()}",
            f"Overall Confidence: {recommendation.overall_confidence:.1%}",
            f"Number of Groups: {len(recommendation.groups)}",
            "",
            "Strategy Distribution:"
        ]
        
        strategy_counts = {}
        for group in recommendation.groups:
            strategy_counts[group.strategy] = strategy_counts.get(group.strategy, 0) + 1
        
        for strategy, count in strategy_counts.items():
            summary_text.append(f"  {strategy.title()}: {count} groups")
        
        summary_text.extend(["", "Top Rationale:"])
        for reason in recommendation.rationale[:3]:
            summary_text.append(f"• {reason}")
        
        ax4.text(0.05, 0.95, '\n'.join(summary_text), transform=ax4.transAxes,
                fontsize=10, verticalalignment='top', fontfamily='monospace')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"Assembly strategy visualization saved to {output_path}")
        
    except ImportError:
        logging.warning("Matplotlib not available - skipping visualization")
    except Exception as e:
        logging.error(f"Failed to create assembly strategy visualization: {e}")