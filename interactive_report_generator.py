#!/usr/bin/env python3
"""
Interactive HTML Report Generator for MetaGrouper

Creates comprehensive, publication-ready interactive reports with:
- Dynamic visualizations
- Explained assembly strategies
- Interactive threshold exploration
- Professional layout with narrative explanations
"""

import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from datetime import datetime
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline as pyo
from jinja2 import Template

# Try to import existing components
try:
    from metagrouper_package.metagrouper.interactive_visualizer import InteractiveVisualizer
except ImportError:
    # Fallback for development
    import sys
    from pathlib import Path
    sys.path.append(str(Path(__file__).parent / "metagrouper_package" / "metagrouper"))
    from interactive_visualizer import InteractiveVisualizer


class InteractiveReportGenerator:
    """Generate comprehensive interactive HTML reports for MetaGrouper analysis."""
    
    def __init__(self, output_dir: str, title: str = "MetaGrouper Analysis Report"):
        """
        Initialize the report generator.
        
        Args:
            output_dir: Directory to save the report
            title: Title for the report
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.title = title
        self.report_data = {}
        self.interactive_plots = {}
        
        logging.info(f"Initialized InteractiveReportGenerator: {title}")
    
    def add_analysis_data(self, 
                         distance_matrix: np.ndarray,
                         sample_names: List[str],
                         metadata: Optional[pd.DataFrame] = None,
                         permanova_results: Optional[pd.DataFrame] = None,
                         assembly_recommendation: Optional[Any] = None,
                         kmer_data: Optional[Dict] = None):
        """
        Add analysis data to the report.
        
        Args:
            distance_matrix: Sample distance matrix
            sample_names: List of sample names
            metadata: Metadata DataFrame
            permanova_results: PERMANOVA analysis results
            assembly_recommendation: Assembly strategy recommendation
            kmer_data: K-mer analysis data
        """
        self.report_data.update({
            'distance_matrix': distance_matrix,
            'sample_names': sample_names,
            'metadata': metadata,
            'permanova_results': permanova_results,
            'assembly_recommendation': assembly_recommendation,
            'kmer_data': kmer_data,
            'timestamp': datetime.now().isoformat(),
            'n_samples': len(sample_names)
        })
        
        logging.info(f"Added analysis data for {len(sample_names)} samples")
    
    def create_interactive_threshold_explorer(self) -> str:
        """Create an interactive threshold exploration tool."""
        
        if 'assembly_recommendation' not in self.report_data:
            return "<p>No assembly recommendation data available</p>"
        
        recommendation = self.report_data['assembly_recommendation']
        distance_matrix = self.report_data['distance_matrix']
        sample_names = self.report_data['sample_names']
        
        # Create threshold exploration plot
        thresholds = np.arange(0.1, 0.8, 0.05)
        threshold_data = []
        
        for threshold in thresholds:
            # Simulate grouping at different thresholds
            n_groups, grouped_samples = self._simulate_grouping(distance_matrix, threshold)
            threshold_data.append({
                'threshold': threshold,
                'n_groups': n_groups,
                'grouped_samples': grouped_samples,
                'individual_samples': len(sample_names) - grouped_samples
            })
        
        df = pd.DataFrame(threshold_data)
        
        # Create the interactive plot
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=('Number of Groups vs Threshold', 'Grouping Strategy Distribution'),
            specs=[[{"secondary_y": False}], [{"type": "bar"}]]
        )
        
        # Line plot for groups
        fig.add_trace(
            go.Scatter(
                x=df['threshold'],
                y=df['n_groups'],
                mode='lines+markers',
                name='Number of Groups',
                line=dict(width=3, color='blue'),
                marker=dict(size=8)
            ),
            row=1, col=1
        )
        
        # Bar plot for strategy distribution
        fig.add_trace(
            go.Bar(
                x=df['threshold'],
                y=df['grouped_samples'],
                name='Grouped Samples',
                marker_color='lightblue'
            ),
            row=2, col=1
        )
        
        fig.add_trace(
            go.Bar(
                x=df['threshold'],
                y=df['individual_samples'],
                name='Individual Samples',
                marker_color='orange'
            ),
            row=2, col=1
        )
        
        # Add current threshold line
        if hasattr(recommendation, 'groups') and recommendation.groups:
            current_threshold = 0.45  # Current default
            fig.add_vline(
                x=current_threshold,
                line_dash="dash",
                line_color="red",
                annotation_text=f"Current: {current_threshold}",
                row=1, col=1
            )
            fig.add_vline(
                x=current_threshold,
                line_dash="dash", 
                line_color="red",
                row=2, col=1
            )
        
        fig.update_layout(
            title="Interactive Threshold Explorer",
            height=600,
            template='plotly_white',
            showlegend=True
        )
        
        fig.update_xaxes(title_text="Similarity Threshold", row=1, col=1)
        fig.update_yaxes(title_text="Number of Groups", row=1, col=1)
        fig.update_xaxes(title_text="Similarity Threshold", row=2, col=1)
        fig.update_yaxes(title_text="Number of Samples", row=2, col=1)
        
        # Convert to HTML
        return fig.to_html(include_plotlyjs='cdn', div_id="threshold-explorer")
    
    def _simulate_grouping(self, distance_matrix: np.ndarray, threshold: float) -> Tuple[int, int]:
        """Simulate grouping at a given threshold."""
        # Simple grouping simulation based on distance threshold
        # This is a simplified version - in practice would use the actual grouping algorithm
        
        similarities = 1 - distance_matrix
        np.fill_diagonal(similarities, 0)  # Remove self-similarities
        
        # Count pairs above threshold
        above_threshold = np.sum(similarities > threshold, axis=1)
        
        # Estimate number of groups (simplified)
        n_samples = len(distance_matrix)
        if np.max(above_threshold) == 0:
            return 0, 0  # No groups, all individual
        
        # Rough estimate of groups and grouped samples
        avg_connections = np.mean(above_threshold[above_threshold > 0])
        estimated_groups = max(1, int(n_samples / (avg_connections + 1)))
        grouped_samples = np.sum(above_threshold > 0)
        
        return estimated_groups, grouped_samples
    
    def create_assembly_strategy_explanation(self) -> str:
        """Create an interactive explanation of the assembly strategy."""
        
        if 'assembly_recommendation' not in self.report_data:
            return "<p>No assembly recommendation available</p>"
        
        recommendation = self.report_data['assembly_recommendation']
        
        # Create decision tree visualization
        decision_tree_html = self._create_decision_tree(recommendation)
        
        # Create strategy comparison table
        strategy_table = self._create_strategy_comparison_table(recommendation)
        
        # Create group details
        group_details = self._create_group_details(recommendation)
        
        explanation_html = f"""
        <div class="assembly-strategy-section">
            <h3>🎯 Assembly Strategy Explanation</h3>
            
            <div class="strategy-overview">
                <div class="strategy-card">
                    <h4>Recommended Strategy: {recommendation.strategy.title()}</h4>
                    <div class="confidence-bar">
                        <div class="confidence-fill" style="width: {recommendation.overall_confidence*100:.1f}%"></div>
                        <span class="confidence-text">{recommendation.overall_confidence:.1%} Confidence</span>
                    </div>
                    <p><strong>Rationale:</strong> {recommendation.decision_rationale}</p>
                </div>
            </div>
            
            <div class="decision-process">
                <h4>📊 Decision Process</h4>
                {decision_tree_html}
            </div>
            
            <div class="strategy-comparison">
                <h4>⚖️ Strategy Comparison</h4>
                {strategy_table}
            </div>
            
            <div class="group-details">
                <h4>👥 Group Details</h4>
                {group_details}
            </div>
        </div>
        """
        
        return explanation_html
    
    def _create_decision_tree(self, recommendation) -> str:
        """Create a visual decision tree showing the assembly strategy logic."""
        
        n_samples = self.report_data['n_samples']
        n_groups = len(recommendation.groups) if recommendation.groups else 0
        
        # Create a simple decision tree using HTML/CSS
        tree_html = f"""
        <div class="decision-tree">
            <div class="tree-node root">
                <div class="node-content">
                    <strong>{n_samples} Samples</strong><br>
                    Start Analysis
                </div>
            </div>
            
            <div class="tree-level">
                <div class="tree-node">
                    <div class="node-content">
                        K-mer Similarity<br>
                        Analysis
                    </div>
                </div>
                
                <div class="tree-node">
                    <div class="node-content">
                        Metadata<br>
                        Analysis
                    </div>
                </div>
            </div>
            
            <div class="tree-level">
                <div class="tree-node decision">
                    <div class="node-content">
                        <strong>Decision</strong><br>
                        {recommendation.strategy.title()}<br>
                        {n_groups} Groups
                    </div>
                </div>
            </div>
        </div>
        """
        
        return tree_html
    
    def _create_strategy_comparison_table(self, recommendation) -> str:
        """Create a comparison table of different assembly strategies."""
        
        n_samples = self.report_data['n_samples']
        
        strategies = [
            {
                'strategy': 'Individual',
                'assemblies': n_samples,
                'contamination_risk': 'Very Low',
                'coverage': 'Low',
                'computational_cost': 'High',
                'best_for': 'Highly diverse samples'
            },
            {
                'strategy': 'Grouped', 
                'assemblies': len(recommendation.groups) if recommendation.groups else 0,
                'contamination_risk': 'Medium',
                'coverage': 'High',
                'computational_cost': 'Medium',
                'best_for': 'Related samples with clear groupings'
            },
            {
                'strategy': 'Global',
                'assemblies': 1,
                'contamination_risk': 'High',
                'coverage': 'Very High',
                'computational_cost': 'Low',
                'best_for': 'Very similar samples'
            }
        ]
        
        table_html = """
        <table class="strategy-table">
            <thead>
                <tr>
                    <th>Strategy</th>
                    <th># Assemblies</th>
                    <th>Contamination Risk</th>
                    <th>Coverage</th>
                    <th>Computational Cost</th>
                    <th>Best For</th>
                </tr>
            </thead>
            <tbody>
        """
        
        for strategy in strategies:
            row_class = 'recommended' if strategy['strategy'].lower() == recommendation.strategy.lower() else ''
            table_html += f"""
                <tr class="{row_class}">
                    <td><strong>{strategy['strategy']}</strong></td>
                    <td>{strategy['assemblies']}</td>
                    <td>{strategy['contamination_risk']}</td>
                    <td>{strategy['coverage']}</td>
                    <td>{strategy['computational_cost']}</td>
                    <td>{strategy['best_for']}</td>
                </tr>
            """
        
        table_html += """
            </tbody>
        </table>
        """
        
        return table_html
    
    def _create_group_details(self, recommendation) -> str:
        """Create detailed information about each assembly group."""
        
        if not recommendation.groups:
            return "<p>No groups formed - individual assembly recommended</p>"
        
        details_html = """
        <div class="groups-container">
        """
        
        for i, group in enumerate(recommendation.groups):
            benefits_list = "".join([f"<li>{benefit}</li>" for benefit in group.expected_benefits])
            challenges_list = "".join([f"<li>{challenge}</li>" for challenge in group.expected_challenges])
            
            details_html += f"""
            <div class="group-card">
                <div class="group-header">
                    <h5>Group {i+1}: {group.group_id}</h5>
                    <div class="confidence-badge">{group.confidence_score:.1%}</div>
                </div>
                
                <div class="group-info">
                    <div class="group-stat">
                        <span class="stat-label">Samples:</span>
                        <span class="stat-value">{len(group.sample_names)}</span>
                    </div>
                    <div class="group-stat">
                        <span class="stat-label">Avg Distance:</span>
                        <span class="stat-value">{group.avg_distance:.3f}</span>
                    </div>
                    <div class="group-stat">
                        <span class="stat-label">Criterion:</span>
                        <span class="stat-value">{group.grouping_criterion}</span>
                    </div>
                </div>
                
                <div class="group-samples">
                    <strong>Samples:</strong> {', '.join(group.sample_names)}
                </div>
                
                <div class="group-assessment">
                    <div class="benefits">
                        <strong>Expected Benefits:</strong>
                        <ul>{benefits_list}</ul>
                    </div>
                    <div class="challenges">
                        <strong>Potential Challenges:</strong>
                        <ul>{challenges_list}</ul>
                    </div>
                </div>
            </div>
            """
        
        details_html += """
        </div>
        """
        
        return details_html
    
    def create_comprehensive_report(self, include_raw_data: bool = False) -> str:
        """
        Create the comprehensive interactive HTML report.
        
        Args:
            include_raw_data: Whether to include raw data downloads
            
        Returns:
            Path to generated HTML report
        """
        logging.info("Creating comprehensive interactive report")
        
        # Generate all interactive components
        visualizations = self._create_all_visualizations()
        threshold_explorer = self.create_interactive_threshold_explorer()
        strategy_explanation = self.create_assembly_strategy_explanation()
        summary_stats = self._create_summary_statistics()
        
        # Create the HTML template
        html_template = self._get_html_template()
        
        # Render the template
        template = Template(html_template)
        html_content = template.render(
            title=self.title,
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            summary_stats=summary_stats,
            visualizations=visualizations,
            threshold_explorer=threshold_explorer,
            strategy_explanation=strategy_explanation,
            include_raw_data=include_raw_data
        )
        
        # Save the report
        report_path = self.output_dir / "interactive_report.html"
        with open(report_path, 'w') as f:
            f.write(html_content)
        
        # Save supporting files
        self._save_supporting_files()
        
        logging.info(f"Comprehensive report saved to {report_path}")
        return str(report_path)
    
    def _create_all_visualizations(self) -> str:
        """Create all interactive visualizations."""
        
        if 'distance_matrix' not in self.report_data:
            return "<p>No visualization data available</p>"
        
        # Initialize interactive visualizer
        visualizer = InteractiveVisualizer(
            self.report_data['sample_names'],
            self.report_data.get('metadata')
        )
        
        # Create PCA if we have k-mer data
        pca_html = ""
        if 'kmer_data' in self.report_data and self.report_data['kmer_data']:
            try:
                from sklearn.decomposition import PCA
                
                # Extract k-mer profiles for PCA
                profiles = self.report_data['kmer_data'].get('profiles', {})
                if profiles:
                    # Convert to matrix
                    sample_names = list(profiles.keys())
                    all_kmers = set()
                    for profile in profiles.values():
                        all_kmers.update(profile.keys())
                    
                    kmer_matrix = np.zeros((len(sample_names), len(all_kmers)))
                    kmer_list = list(all_kmers)
                    
                    for i, sample in enumerate(sample_names):
                        for j, kmer in enumerate(kmer_list):
                            kmer_matrix[i, j] = profiles[sample].get(kmer, 0)
                    
                    # Perform PCA
                    pca = PCA(n_components=2)
                    pca_result = pca.fit_transform(kmer_matrix)
                    
                    # Create enhanced dimensionality plot
                    pca_html = visualizer.create_enhanced_dimensionality_plot(
                        kmer_matrix,
                        str(self.output_dir / "temp_pca.html"),
                        "Interactive K-mer Analysis"
                    )
                    
                    # Read the HTML content
                    with open(self.output_dir / "temp_pca.html", 'r') as f:
                        pca_html = f.read()
                    
                    # Extract just the plot div
                    start = pca_html.find('<div id="')
                    end = pca_html.find('</script>') + 9
                    if start != -1 and end != -1:
                        pca_html = pca_html[start:end]
                    
            except Exception as e:
                logging.warning(f"Could not create PCA plot: {e}")
                pca_html = "<p>PCA visualization not available</p>"
        
        # Create distance heatmap
        heatmap_html = visualizer.create_interactive_heatmap(
            self.report_data['distance_matrix'],
            str(self.output_dir / "temp_heatmap.html"),
            "Sample Distance Matrix"
        )
        
        # Read and extract heatmap HTML
        try:
            with open(self.output_dir / "temp_heatmap.html", 'r') as f:
                heatmap_content = f.read()
            
            # Extract just the plot div
            start = heatmap_content.find('<div id="')
            end = heatmap_content.find('</script>') + 9
            if start != -1 and end != -1:
                heatmap_html = heatmap_content[start:end]
                
        except Exception as e:
            logging.warning(f"Could not extract heatmap HTML: {e}")
            heatmap_html = "<p>Heatmap visualization not available</p>"
        
        # Combine visualizations
        viz_html = f"""
        <div class="visualizations-container">
            <div class="viz-section">
                <h4>📊 K-mer Principal Component Analysis</h4>
                <div class="viz-content">
                    {pca_html}
                </div>
            </div>
            
            <div class="viz-section">
                <h4>🔥 Sample Distance Heatmap</h4>
                <div class="viz-content">
                    {heatmap_html}
                </div>
            </div>
        </div>
        """
        
        return viz_html
    
    def _create_summary_statistics(self) -> Dict[str, Any]:
        """Create summary statistics for the report."""
        
        stats = {
            'n_samples': self.report_data.get('n_samples', 0),
            'analysis_date': datetime.now().strftime("%Y-%m-%d"),
            'has_metadata': self.report_data.get('metadata') is not None,
            'n_metadata_vars': len(self.report_data['metadata'].columns) if self.report_data.get('metadata') is not None else 0,
        }
        
        # Distance matrix stats
        if 'distance_matrix' in self.report_data:
            dm = self.report_data['distance_matrix']
            stats.update({
                'mean_distance': np.mean(dm),
                'min_distance': np.min(dm[dm > 0]),
                'max_distance': np.max(dm),
                'distance_std': np.std(dm)
            })
        
        # Assembly recommendation stats
        if 'assembly_recommendation' in self.report_data:
            rec = self.report_data['assembly_recommendation']
            stats.update({
                'strategy': rec.strategy,
                'n_groups': len(rec.groups) if rec.groups else 0,
                'confidence': rec.overall_confidence,
                'primary_criterion': rec.primary_criterion
            })
        
        return stats
    
    def _get_html_template(self) -> str:
        """Get the HTML template for the report."""
        
        return '''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{ title }}</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        /* Modern, clean styling for the report */
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background-color: #f8f9fa;
        }
        
        .container {
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }
        
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }
        
        .header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
            font-weight: 300;
        }
        
        .header .subtitle {
            font-size: 1.2em;
            opacity: 0.9;
        }
        
        .section {
            background: white;
            margin-bottom: 30px;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }
        
        .section h2 {
            color: #4a5568;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 1.8em;
        }
        
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }
        
        .stat-card {
            background: #f7fafc;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
            text-align: center;
        }
        
        .stat-value {
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
            display: block;
        }
        
        .stat-label {
            color: #718096;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        
        .confidence-bar {
            background: #e2e8f0;
            height: 20px;
            border-radius: 10px;
            overflow: hidden;
            position: relative;
            margin: 10px 0;
        }
        
        .confidence-fill {
            background: linear-gradient(90deg, #48bb78, #38a169);
            height: 100%;
            transition: width 0.3s ease;
        }
        
        .confidence-text {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            color: white;
            font-weight: bold;
            font-size: 0.9em;
        }
        
        .strategy-card {
            background: #edf2f7;
            padding: 20px;
            border-radius: 8px;
            border: 2px solid #cbd5e0;
            margin-bottom: 20px;
        }
        
        .strategy-table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }
        
        .strategy-table th,
        .strategy-table td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e2e8f0;
        }
        
        .strategy-table th {
            background: #f7fafc;
            font-weight: bold;
            color: #4a5568;
        }
        
        .strategy-table .recommended {
            background: #f0fff4;
            border-left: 4px solid #48bb78;
        }
        
        .groups-container {
            display: grid;
            gap: 20px;
        }
        
        .group-card {
            border: 1px solid #e2e8f0;
            border-radius: 8px;
            padding: 20px;
            background: white;
        }
        
        .group-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 15px;
        }
        
        .confidence-badge {
            background: #667eea;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            font-size: 0.9em;
            font-weight: bold;
        }
        
        .group-info {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px;
            margin-bottom: 15px;
        }
        
        .group-stat {
            display: flex;
            justify-content: space-between;
        }
        
        .stat-label {
            font-weight: bold;
            color: #718096;
        }
        
        .group-assessment {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin-top: 15px;
        }
        
        .benefits ul {
            color: #38a169;
        }
        
        .challenges ul {
            color: #e53e3e;
        }
        
        .viz-section {
            margin-bottom: 40px;
        }
        
        .viz-content {
            background: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }
        
        .decision-tree {
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
        }
        
        .tree-level {
            display: flex;
            gap: 40px;
        }
        
        .tree-node {
            background: #edf2f7;
            border: 2px solid #cbd5e0;
            border-radius: 8px;
            padding: 15px;
            text-align: center;
            min-width: 120px;
        }
        
        .tree-node.root {
            background: #667eea;
            color: white;
            border-color: #5a67d8;
        }
        
        .tree-node.decision {
            background: #48bb78;
            color: white;
            border-color: #38a169;
        }
        
        .footer {
            text-align: center;
            padding: 20px;
            color: #718096;
            border-top: 1px solid #e2e8f0;
            margin-top: 40px;
        }
        
        @media (max-width: 768px) {
            .container {
                padding: 10px;
            }
            
            .header h1 {
                font-size: 2em;
            }
            
            .stats-grid {
                grid-template-columns: 1fr;
            }
            
            .group-assessment {
                grid-template-columns: 1fr;
            }
        }
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <h1>🧬 {{ title }}</h1>
            <div class="subtitle">Comprehensive Metagenomic Analysis Report</div>
            <div class="subtitle">Generated on {{ timestamp }}</div>
        </div>
        
        <!-- Summary Statistics -->
        <div class="section">
            <h2>📊 Analysis Summary</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <span class="stat-value">{{ summary_stats.n_samples }}</span>
                    <span class="stat-label">Samples Analyzed</span>
                </div>
                {% if summary_stats.strategy %}
                <div class="stat-card">
                    <span class="stat-value">{{ summary_stats.strategy.title() }}</span>
                    <span class="stat-label">Assembly Strategy</span>
                </div>
                {% endif %}
                {% if summary_stats.confidence %}
                <div class="stat-card">
                    <span class="stat-value">{{ "%.1f"|format(summary_stats.confidence * 100) }}%</span>
                    <span class="stat-label">Confidence</span>
                </div>
                {% endif %}
                {% if summary_stats.n_groups %}
                <div class="stat-card">
                    <span class="stat-value">{{ summary_stats.n_groups }}</span>
                    <span class="stat-label">Assembly Groups</span>
                </div>
                {% endif %}
                {% if summary_stats.mean_distance %}
                <div class="stat-card">
                    <span class="stat-value">{{ "%.3f"|format(summary_stats.mean_distance) }}</span>
                    <span class="stat-label">Mean Distance</span>
                </div>
                {% endif %}
                {% if summary_stats.has_metadata %}
                <div class="stat-card">
                    <span class="stat-value">{{ summary_stats.n_metadata_vars }}</span>
                    <span class="stat-label">Metadata Variables</span>
                </div>
                {% endif %}
            </div>
        </div>
        
        <!-- Interactive Visualizations -->
        <div class="section">
            <h2>📈 Interactive Visualizations</h2>
            {{ visualizations|safe }}
        </div>
        
        <!-- Threshold Explorer -->
        <div class="section">
            <h2>🎚️ Threshold Explorer</h2>
            <p>Explore how different similarity thresholds affect sample grouping decisions:</p>
            {{ threshold_explorer|safe }}
        </div>
        
        <!-- Assembly Strategy Explanation -->
        <div class="section">
            <h2>🎯 Assembly Strategy</h2>
            {{ strategy_explanation|safe }}
        </div>
        
        <!-- Footer -->
        <div class="footer">
            <p>Generated by MetaGrouper | 🤖 Enhanced with research-based thresholds</p>
            <p>Report created on {{ timestamp }}</p>
        </div>
    </div>
</body>
</html>
        '''
    
    def _save_supporting_files(self):
        """Save supporting files and data."""
        
        # Save configuration used
        config_data = {
            'analysis_timestamp': self.report_data.get('timestamp'),
            'n_samples': self.report_data.get('n_samples'),
            'similarity_thresholds': {
                'high': 0.25,
                'medium': 0.45,
                'default': 0.45
            },
            'version': 'MetaGrouper v2.0 (Enhanced)',
            'report_generator_version': '1.0'
        }
        
        with open(self.output_dir / 'analysis_config.json', 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Save sample names
        if 'sample_names' in self.report_data:
            with open(self.output_dir / 'sample_names.txt', 'w') as f:
                for sample in self.report_data['sample_names']:
                    f.write(f"{sample}\n")
        
        logging.info("Supporting files saved")


def create_interactive_report(distance_matrix: np.ndarray,
                            sample_names: List[str],
                            output_dir: str,
                            metadata: Optional[pd.DataFrame] = None,
                            permanova_results: Optional[pd.DataFrame] = None,
                            assembly_recommendation: Optional[Any] = None,
                            kmer_data: Optional[Dict] = None,
                            title: str = "MetaGrouper Analysis Report") -> str:
    """
    Convenience function to create a comprehensive interactive report.
    
    Args:
        distance_matrix: Sample distance matrix
        sample_names: List of sample names
        output_dir: Directory to save the report
        metadata: Optional metadata DataFrame
        permanova_results: Optional PERMANOVA results
        assembly_recommendation: Optional assembly recommendation
        kmer_data: Optional k-mer analysis data
        title: Report title
        
    Returns:
        Path to generated HTML report
    """
    
    generator = InteractiveReportGenerator(output_dir, title)
    
    generator.add_analysis_data(
        distance_matrix=distance_matrix,
        sample_names=sample_names,
        metadata=metadata,
        permanova_results=permanova_results,
        assembly_recommendation=assembly_recommendation,
        kmer_data=kmer_data
    )
    
    return generator.create_comprehensive_report()


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)
    
    # Create example data
    n_samples = 6
    sample_names = [f"sample_{i:02d}" for i in range(n_samples)]
    distance_matrix = np.random.rand(n_samples, n_samples)
    np.fill_diagonal(distance_matrix, 0)
    distance_matrix = (distance_matrix + distance_matrix.T) / 2  # Make symmetric
    
    # Create example metadata
    metadata = pd.DataFrame({
        'sample_id': sample_names,
        'treatment': ['control', 'control', 'treated', 'treated', 'control', 'treated'],
        'timepoint': ['baseline', 'week4', 'baseline', 'week4', 'baseline', 'week4'],
        'batch': ['A', 'A', 'B', 'B', 'A', 'B']
    }).set_index('sample_id')
    
    # Create a mock assembly recommendation for testing
    from dataclasses import dataclass
    
    @dataclass 
    class MockAssemblyGroup:
        group_id: str = "test_group"
        sample_names: list = None
        grouping_criterion: str = "similarity"
        criterion_value: str = "0.3"
        avg_distance: float = 0.3
        max_distance: float = 0.4
        confidence_score: float = 0.8
        expected_benefits: list = None
        expected_challenges: list = None
        
        def __post_init__(self):
            if self.sample_names is None:
                self.sample_names = sample_names[:3]
            if self.expected_benefits is None:
                self.expected_benefits = ["Better coverage", "Improved assembly"]
            if self.expected_challenges is None:
                self.expected_challenges = ["Potential contamination"]
    
    @dataclass
    class MockAssemblyRecommendation:
        strategy: str = "grouped"
        groups: list = None
        overall_confidence: float = 0.75
        primary_criterion: str = "similarity"
        decision_rationale: str = "Samples show moderate similarity"
        assembly_commands: dict = None
        performance_predictions: dict = None
        
        def __post_init__(self):
            if self.groups is None:
                self.groups = [MockAssemblyGroup()]
            if self.assembly_commands is None:
                self.assembly_commands = {"megahit": {}, "spades": {}}
            if self.performance_predictions is None:
                self.performance_predictions = {}
    
    mock_recommendation = MockAssemblyRecommendation()
    
    # Generate report
    report_path = create_interactive_report(
        distance_matrix=distance_matrix,
        sample_names=sample_names,
        output_dir="example_interactive_report",
        metadata=metadata,
        assembly_recommendation=mock_recommendation,
        title="Example MetaGrouper Report"
    )
    
    print(f"Example report generated: {report_path}")