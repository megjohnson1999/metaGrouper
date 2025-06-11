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
class AssemblyRecommendation:
    """Complete assembly strategy recommendation."""
    strategy: str  # 'individual', 'grouped', 'global'
    groups: List[AssemblyGroup]
    overall_confidence: float
    primary_criterion: Optional[str]
    decision_rationale: str
    assembly_commands: Dict[str, List[str]]
    performance_predictions: Dict[str, Any]

class AssemblyStrategyEngine:
    """Core engine for determining optimal assembly strategies."""
    
    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str]):
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names
        self.n_samples = len(sample_names)
        
        # Thresholds for assembly decisions (can be tuned)
        self.similarity_threshold_high = 0.15  # Very similar samples
        self.similarity_threshold_medium = 0.30  # Moderately similar samples
        self.significance_threshold = 0.05  # P-value threshold for metadata
        self.min_group_size = 2
        self.max_group_size = 10
        
    def _calculate_group_statistics(self, group_indices: List[int]) -> Tuple[float, float]:
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
        """Assess the quality of a proposed assembly group."""
        if len(group_indices) < 2:
            return 0.0
        
        avg_dist, max_dist = self._calculate_group_statistics(group_indices)
        
        # Calculate confidence based on distance metrics
        # Lower distances = higher confidence
        distance_score = max(0, 1 - (avg_dist / self.similarity_threshold_medium))
        
        # Size penalty/bonus
        size_score = 1.0
        if len(group_indices) < 3:
            size_score = 0.8  # Small groups less reliable
        elif len(group_indices) > 8:
            size_score = 0.7  # Very large groups may be problematic
        
        # Homogeneity score (max distance within group)
        homogeneity_score = max(0, 1 - (max_dist / self.similarity_threshold_medium))
        
        confidence = (distance_score * 0.4 + size_score * 0.2 + homogeneity_score * 0.4)
        return min(1.0, max(0.0, confidence))
    
    def recommend_by_similarity(self) -> List[AssemblyGroup]:
        """Recommend assembly groups based purely on k-mer similarity."""
        logging.info("Generating similarity-based assembly recommendations")
        
        # Use hierarchical clustering to identify natural groups
        condensed_distances = squareform(self.distance_matrix, checks=False)
        linkage_matrix = linkage(condensed_distances, method='average')
        
        groups = []
        
        # Try different numbers of clusters
        for n_clusters in range(2, min(self.n_samples, 8)):
            cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
            
            for cluster_id in range(1, n_clusters + 1):
                cluster_indices = [i for i, label in enumerate(cluster_labels) if label == cluster_id]
                
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
                            "Enhanced detection of strain variants"
                        ],
                        expected_challenges=[
                            "Potential strain mixing",
                            "Increased computational requirements",
                            f"Assembly size may be {len(cluster_indices)}x larger"
                        ]
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
    
    def recommend_by_metadata(self, metadata_results: pd.DataFrame, 
                             metadata: pd.DataFrame) -> List[AssemblyGroup]:
        """Recommend assembly groups based on metadata analysis."""
        logging.info("Generating metadata-based assembly recommendations")
        
        if metadata_results.empty or metadata is None:
            return []
        
        groups = []
        
        # Focus on significant and high-explaining variables
        significant_vars = metadata_results[
            (metadata_results['p_value'] < self.significance_threshold) |
            (metadata_results['r_squared'] > 0.20)  # High explanatory power
        ].sort_values('r_squared', ascending=False)
        
        for _, row in significant_vars.iterrows():
            variable = row['variable']
            
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
                sample_indices = [self.sample_names.index(name) for name in value_samples.index 
                                if name in self.sample_names]
                
                if len(sample_indices) < self.min_group_size:
                    continue
                
                avg_dist, max_dist = self._calculate_group_statistics(sample_indices)
                
                # Calculate confidence based on statistical significance and distance
                stat_confidence = min(1.0, row['r_squared'] * 2)  # R-squared based
                p_value_bonus = 1.2 if row['p_value'] < 0.01 else 1.0 if row['p_value'] < 0.05 else 0.8
                distance_penalty = max(0.5, 1 - (avg_dist / self.similarity_threshold_medium))
                
                confidence = stat_confidence * p_value_bonus * distance_penalty
                confidence = min(1.0, max(0.0, confidence))
                
                # Determine benefits and challenges
                benefits = [
                    f"Biologically meaningful grouping by {variable}",
                    "Reduced inter-sample contamination",
                    "Better representation of group-specific features"
                ]
                
                challenges = [
                    "May miss cross-group shared sequences",
                    f"Groups based on {variable} may have variable quality"
                ]
                
                if avg_dist > self.similarity_threshold_medium:
                    challenges.append(f"High within-group diversity (avg dist: {avg_dist:.3f})")
                
                group = AssemblyGroup(
                    group_id=f"{variable}_{value}",
                    sample_names=list(value_samples.index),
                    grouping_criterion=variable,
                    criterion_value=value,
                    avg_distance=avg_dist,
                    max_distance=max_dist,
                    confidence_score=confidence,
                    expected_benefits=benefits,
                    expected_challenges=challenges
                )
                groups.append(group)
        
        return groups
    
    def recommend_hybrid_strategy(self, similarity_groups: List[AssemblyGroup],
                                 metadata_groups: List[AssemblyGroup]) -> List[AssemblyGroup]:
        """Combine similarity and metadata-based recommendations."""
        logging.info("Generating hybrid assembly recommendations")
        
        # Score and combine approaches
        all_groups = similarity_groups + metadata_groups
        
        # Re-score considering both approaches
        for group in all_groups:
            if group.grouping_criterion == "k-mer_similarity":
                # Bonus for purely similarity-based groups with very low distances
                if group.avg_distance < self.similarity_threshold_high:
                    group.confidence_score *= 1.2
            else:
                # Metadata-based groups get bonus if they also have low distances
                if group.avg_distance < self.similarity_threshold_medium:
                    group.confidence_score *= 1.1
        
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
        
    def generate_megahit_commands(self, groups: List[AssemblyGroup]) -> Dict[str, List[str]]:
        """Generate MEGAHIT co-assembly commands."""
        commands = {}
        
        for group in groups:
            if len(group.sample_names) == 1:
                # Individual assembly
                sample = group.sample_names[0]
                cmd = f"megahit -r {sample}.fastq -o {sample}_assembly --min-contig-len 500"
            else:
                # Co-assembly
                input_files = ",".join([f"{sample}.fastq" for sample in group.sample_names])
                output_dir = f"{group.group_id}_coassembly"
                cmd = f"megahit -r {input_files} -o {output_dir} --min-contig-len 500 --k-list 21,29,39,59,79,99"
            
            commands[group.group_id] = [cmd]
        
        return commands
    
    def generate_spades_commands(self, groups: List[AssemblyGroup]) -> Dict[str, List[str]]:
        """Generate SPAdes metagenomic assembly commands."""
        commands = {}
        
        for group in groups:
            if len(group.sample_names) == 1:
                # Individual assembly
                sample = group.sample_names[0]
                cmd = f"spades.py --meta -s {sample}.fastq -o {sample}_spades_assembly"
            else:
                # Co-assembly (combine reads first)
                group_name = group.group_id
                combine_cmd = f"cat {' '.join([f'{sample}.fastq' for sample in group.sample_names])} > {group_name}_combined.fastq"
                assembly_cmd = f"spades.py --meta -s {group_name}_combined.fastq -o {group_name}_spades_assembly"
                commands[group.group_id] = [combine_cmd, assembly_cmd]
                continue
            
            commands[group.group_id] = [cmd]
        
        return commands
    
    def generate_flye_commands(self, groups: List[AssemblyGroup]) -> Dict[str, List[str]]:
        """Generate Flye assembly commands (for long reads)."""
        commands = {}
        
        for group in groups:
            if len(group.sample_names) == 1:
                sample = group.sample_names[0]
                cmd = f"flye --meta --nano-raw {sample}.fastq -o {sample}_flye_assembly"
            else:
                # Co-assembly
                group_name = group.group_id
                combine_cmd = f"cat {' '.join([f'{sample}.fastq' for sample in group.sample_names])} > {group_name}_combined.fastq"
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
        predictions['strategy_summary'] = {
            'total_assemblies': len(groups) + individual_samples,
            'co_assemblies': len(groups),
            'individual_assemblies': individual_samples,
            'samples_in_coassembly': grouped_samples,
            'coassembly_percentage': (grouped_samples / total_samples) * 100
        }
        
        # Per-group predictions
        group_predictions = {}
        for group in groups:
            n_samples = len(group.sample_names)
            avg_dist = group.avg_distance
            
            # Predict relative benefits/challenges
            expected_contiguity_improvement = max(1.0, 2.0 - avg_dist * 3)  # Lower distance = better contiguity
            expected_coverage_boost = min(n_samples * 0.8, 5.0)  # Diminishing returns
            computational_cost_multiplier = n_samples ** 1.5  # Non-linear scaling
            
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
                'n_samples': n_samples,
                'avg_distance': avg_dist,
                'predicted_quality': predicted_quality,
                'contamination_risk': contamination_risk,
                'expected_contiguity_improvement': f"{expected_contiguity_improvement:.1f}x",
                'expected_coverage_boost': f"{expected_coverage_boost:.1f}x",
                'computational_cost_multiplier': f"{computational_cost_multiplier:.1f}x",
                'confidence': group.confidence_score
            }
        
        predictions['group_predictions'] = group_predictions
        
        return predictions

class AssemblyRecommender:
    """Main class for generating comprehensive assembly recommendations."""
    
    def __init__(self, distance_matrix: np.ndarray, sample_names: List[str]):
        self.distance_matrix = distance_matrix
        self.sample_names = sample_names
        self.strategy_engine = AssemblyStrategyEngine(distance_matrix, sample_names)
        self.command_generator = AssemblyCommandGenerator()
        self.performance_predictor = PerformancePredictor(distance_matrix, sample_names)
    
    def generate_recommendations(self, metadata_results: Optional[pd.DataFrame] = None,
                               metadata: Optional[pd.DataFrame] = None) -> AssemblyRecommendation:
        """Generate comprehensive assembly recommendations."""
        logging.info("Generating comprehensive assembly recommendations")
        
        # Get recommendations from different approaches
        similarity_groups = self.strategy_engine.recommend_by_similarity()
        
        metadata_groups = []
        if metadata_results is not None and not metadata_results.empty:
            metadata_groups = self.strategy_engine.recommend_by_metadata(metadata_results, metadata)
        
        # Combine approaches
        if metadata_groups:
            final_groups = self.strategy_engine.recommend_hybrid_strategy(similarity_groups, metadata_groups)
        else:
            final_groups = similarity_groups
        
        # Determine overall strategy
        total_samples = len(self.sample_names)
        grouped_samples = sum(len(group.sample_names) for group in final_groups)
        
        if not final_groups or grouped_samples < 2:
            strategy = "individual"
            rationale = "No clear grouping patterns found. Individual assembly recommended."
            overall_confidence = 0.7
        elif grouped_samples == total_samples and len(final_groups) == 1:
            strategy = "global"
            rationale = "All samples show strong similarity. Global co-assembly recommended."
            overall_confidence = final_groups[0].confidence_score
        else:
            strategy = "grouped"
            rationale = f"Mixed strategy: {len(final_groups)} co-assembly groups covering {grouped_samples}/{total_samples} samples."
            overall_confidence = np.mean([group.confidence_score for group in final_groups])
        
        # Determine primary criterion
        primary_criterion = None
        if final_groups:
            # Find most common grouping criterion
            criteria = [group.grouping_criterion for group in final_groups]
            primary_criterion = max(set(criteria), key=criteria.count)
        
        # Generate assembly commands
        assembly_commands = {
            'megahit': self.command_generator.generate_megahit_commands(final_groups),
            'spades': self.command_generator.generate_spades_commands(final_groups),
            'flye': self.command_generator.generate_flye_commands(final_groups)
        }
        
        # Predict performance
        performance_predictions = self.performance_predictor.predict_assembly_metrics(final_groups)
        
        recommendation = AssemblyRecommendation(
            strategy=strategy,
            groups=final_groups,
            overall_confidence=overall_confidence,
            primary_criterion=primary_criterion,
            decision_rationale=rationale,
            assembly_commands=assembly_commands,
            performance_predictions=performance_predictions
        )
        
        return recommendation

def save_recommendations(recommendation: AssemblyRecommendation, output_path: str):
    """Save assembly recommendations to files."""
    output_dir = Path(output_path)
    output_dir.mkdir(exist_ok=True)
    
    # Save detailed recommendation as JSON
    rec_dict = asdict(recommendation)
    with open(output_dir / 'assembly_recommendations.json', 'w') as f:
        json.dump(rec_dict, f, indent=2, default=str)
    
    # Save human-readable summary
    with open(output_dir / 'assembly_strategy.md', 'w') as f:
        f.write("# MetaGrouper Assembly Strategy Recommendations\n\n")
        
        f.write(f"## Overall Strategy: {recommendation.strategy.title()}\n\n")
        f.write(f"**Confidence Score:** {recommendation.overall_confidence:.2f}\n\n")
        f.write(f"**Rationale:** {recommendation.decision_rationale}\n\n")
        
        if recommendation.primary_criterion:
            f.write(f"**Primary Grouping Criterion:** {recommendation.primary_criterion}\n\n")
        
        # Strategy summary
        perf = recommendation.performance_predictions['strategy_summary']
        f.write("## Strategy Summary\n\n")
        f.write(f"- **Total Assemblies:** {perf['total_assemblies']}\n")
        f.write(f"- **Co-assemblies:** {perf['co_assemblies']}\n")
        f.write(f"- **Individual Assemblies:** {perf['individual_assemblies']}\n")
        f.write(f"- **Samples in Co-assembly:** {perf['samples_in_coassembly']} ({perf['coassembly_percentage']:.1f}%)\n\n")
        
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
                
                f.write("**Expected Benefits:**\n")
                for benefit in group.expected_benefits:
                    f.write(f"- {benefit}\n")
                f.write("\n")
                
                f.write("**Expected Challenges:**\n")
                for challenge in group.expected_challenges:
                    f.write(f"- {challenge}\n")
                f.write("\n")
        
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
        script_content += f"# Generated assembly strategy: {recommendation.strategy}\n\n"
        
        for group_id, cmd_list in commands.items():
            script_content += f"# Group: {group_id}\n"
            for cmd in cmd_list:
                script_content += f"{cmd}\n"
            script_content += "\n"
        
        with open(output_dir / f'run_{tool}_assemblies.sh', 'w') as f:
            f.write(script_content)
    
    logging.info(f"Assembly recommendations saved to {output_path}")

def visualize_assembly_strategy(recommendation: AssemblyRecommendation, 
                              distance_matrix: np.ndarray, 
                              sample_names: List[str],
                              output_path: str):
    """Create visualization of the recommended assembly strategy."""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Distance matrix with group overlays
    im1 = ax1.imshow(distance_matrix, cmap='viridis', aspect='auto')
    ax1.set_title('Sample Distance Matrix with Assembly Groups')
    ax1.set_xticks(range(len(sample_names)))
    ax1.set_yticks(range(len(sample_names)))
    ax1.set_xticklabels(sample_names, rotation=45, ha='right')
    ax1.set_yticklabels(sample_names)
    
    # Overlay group boundaries
    colors = plt.cm.Set1(np.linspace(0, 1, len(recommendation.groups)))
    used_samples = set()
    
    for i, (group, color) in enumerate(zip(recommendation.groups, colors)):
        group_indices = [sample_names.index(name) for name in group.sample_names if name not in used_samples]
        if len(group_indices) > 1:
            for idx in group_indices:
                ax1.axhline(y=idx-0.5, color=color, linewidth=2, alpha=0.7)
                ax1.axvline(x=idx-0.5, color=color, linewidth=2, alpha=0.7)
        used_samples.update(group.sample_names)
    
    plt.colorbar(im1, ax=ax1, label='Distance')
    
    # 2. Group confidence scores
    if recommendation.groups:
        group_names = [f"{group.group_id}" for group in recommendation.groups]
        confidences = [group.confidence_score for group in recommendation.groups]
        
        bars = ax2.bar(range(len(group_names)), confidences, 
                      color=colors[:len(group_names)], alpha=0.7)
        ax2.set_title('Assembly Group Confidence Scores')
        ax2.set_xlabel('Assembly Groups')
        ax2.set_ylabel('Confidence Score')
        ax2.set_xticks(range(len(group_names)))
        ax2.set_xticklabels(group_names, rotation=45, ha='right')
        ax2.set_ylim(0, 1)
        ax2.grid(axis='y', alpha=0.3)
        
        # Add confidence values on bars
        for bar, conf in zip(bars, confidences):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                    f'{conf:.2f}', ha='center', va='bottom')
    else:
        ax2.text(0.5, 0.5, 'No groups recommended\nIndividual assembly suggested', 
                ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('Assembly Group Confidence Scores')
    
    # 3. Strategy overview pie chart
    perf = recommendation.performance_predictions['strategy_summary']
    
    labels = []
    sizes = []
    colors_pie = []
    
    if perf['co_assemblies'] > 0:
        labels.append(f"Co-assemblies ({perf['co_assemblies']})")
        sizes.append(perf['co_assemblies'])
        colors_pie.append('lightblue')
    
    if perf['individual_assemblies'] > 0:
        labels.append(f"Individual ({perf['individual_assemblies']})")
        sizes.append(perf['individual_assemblies'])
        colors_pie.append('lightcoral')
    
    if sizes:
        ax3.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%', startangle=90)
    ax3.set_title(f'Assembly Strategy Distribution\n({recommendation.strategy.title()})')
    
    # 4. Performance predictions
    if recommendation.groups:
        group_sizes = [len(group.sample_names) for group in recommendation.groups]
        avg_distances = [group.avg_distance for group in recommendation.groups]
        
        scatter = ax4.scatter(group_sizes, avg_distances, 
                            c=confidences, cmap='RdYlGn', 
                            s=100, alpha=0.7, edgecolors='black')
        
        ax4.set_xlabel('Group Size (number of samples)')
        ax4.set_ylabel('Average Intra-group Distance')
        ax4.set_title('Group Size vs. Distance (colored by confidence)')
        ax4.grid(True, alpha=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax4)
        cbar.set_label('Confidence Score')
        
        # Annotate points
        for i, group in enumerate(recommendation.groups):
            ax4.annotate(f'G{i+1}', (group_sizes[i], avg_distances[i]), 
                        xytext=(5, 5), textcoords='offset points')
    else:
        ax4.text(0.5, 0.5, 'No groups to analyze', 
                ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title('Group Size vs. Distance Analysis')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Assembly strategy visualization saved to {output_path}")