#!/usr/bin/env python3
"""
Interactive HTML Report Generator for MetaGrouper

Creates comprehensive, publication-ready interactive reports with:
- Dynamic visualizations
- Explained assembly strategies
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
                         kmer_data: Optional[Dict] = None,
                         grouping_recommendations: Optional[List[Dict]] = None,
                         sample_id_column: str = "sample_id",
                         analyzed_variables: Optional[List[str]] = None):
        """
        Add analysis data to the report.
        
        Args:
            distance_matrix: Sample distance matrix
            sample_names: List of sample names
            metadata: Metadata DataFrame
            permanova_results: PERMANOVA analysis results
            assembly_recommendation: Assembly strategy recommendation
            kmer_data: K-mer analysis data
            grouping_recommendations: Intelligent metadata grouping recommendations
            sample_id_column: Column name for sample IDs in metadata
            analyzed_variables: Optional list of specific variables to include in plot dropdowns
        """
        self.report_data.update({
            'distance_matrix': distance_matrix,
            'sample_names': sample_names,
            'metadata': metadata,
            'permanova_results': permanova_results,
            'assembly_recommendation': assembly_recommendation,
            'kmer_data': kmer_data,
            'grouping_recommendations': grouping_recommendations,
            'sample_id_column': sample_id_column,
            'analyzed_variables': analyzed_variables,
            'timestamp': datetime.now().isoformat(),
            'n_samples': len(sample_names)
        })
        
        logging.info(f"Added analysis data for {len(sample_names)} samples")
    
    
    def create_assembly_strategy_explanation(self) -> str:
        """Create an interactive explanation of the assembly strategy."""
        
        if 'assembly_recommendation' not in self.report_data or self.report_data['assembly_recommendation'] is None:
            return "<p>No assembly recommendation available</p>"
        
        recommendation = self.report_data['assembly_recommendation']
        
        # Create decision tree visualization
        decision_tree_html = self._create_decision_tree(recommendation)
        
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
    
    
    def _create_group_details(self, recommendation) -> str:
        """Create detailed information about each assembly group."""
        
        if not recommendation.groups:
            return "<p>No groups formed - individual assembly recommended</p>"
        
        details_html = """
        <div class="groups-container">
        """
        
        for i, group in enumerate(recommendation.groups):
            
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
                
            </div>
            """
        
        details_html += """
        </div>
        """
        
        return details_html
    
    def create_permanova_section(self) -> str:
        """Create PERMANOVA analysis section for the report."""
        
        if 'permanova_results' not in self.report_data or self.report_data['permanova_results'] is None:
            return """
            <div class="permanova-explanation">
                <p><strong>⚠️ No PERMANOVA results available</strong></p>
                <p>PERMANOVA (Permutational Multivariate Analysis of Variance) analysis was not performed or no results were generated. This analysis helps identify which metadata variables significantly explain differences in sample composition.</p>
            </div>
            """
        
        permanova_df = self.report_data['permanova_results']
        
        if permanova_df.empty:
            return """
            <div class="permanova-explanation">
                <p><strong>⚠️ No significant PERMANOVA results found</strong></p>
                <p>None of the metadata variables showed significant associations with sample composition differences.</p>
            </div>
            """
        
        # Create PERMANOVA results visualization
        permanova_html = f"""
        <div class="permanova-section">
            <div class="permanova-explanation">
                <h4>📊 What is PERMANOVA?</h4>
                <p><strong>PERMANOVA</strong> (Permutational Multivariate Analysis of Variance) tests which metadata variables significantly explain differences in sample composition. Higher R² values indicate variables that better explain sample groupings.</p>
                
                <div class="permanova-guide">
                    <div class="guide-item">
                        <strong>R² (Effect Size):</strong> Proportion of variation explained (higher = more important)
                    </div>
                    <div class="guide-item">
                        <strong>p-value:</strong> Statistical significance (< 0.05 = significant)
                    </div>
                    <div class="guide-item">
                        <strong>Significance:</strong> *** p<0.001, ** p<0.01, * p<0.05, . p<0.1
                    </div>
                </div>
            </div>
            
            <div class="permanova-results">
                <h4>🎯 Variable Importance Rankings</h4>
                <div class="variables-container">
        """
        
        # Sort by R-squared (most important first)
        sorted_results = permanova_df.dropna(subset=['r_squared']).sort_values('r_squared', ascending=False)
        
        for _, row in sorted_results.iterrows():
            # Determine significance level
            p_val = row['p_value']
            if p_val < 0.001:
                significance = "***"
                sig_class = "highly-significant"
            elif p_val < 0.01:
                significance = "**" 
                sig_class = "very-significant"
            elif p_val < 0.05:
                significance = "*"
                sig_class = "significant"
            elif p_val < 0.1:
                significance = "."
                sig_class = "marginally-significant"
            else:
                significance = ""
                sig_class = "not-significant"
            
            # Create progress bar for R-squared
            r_squared_percent = row['r_squared'] * 100
            
            permanova_html += f"""
                <div class="variable-card {sig_class}">
                    <div class="variable-header">
                        <h5>{row['variable']}</h5>
                        <div class="significance-badge {sig_class}">{significance if significance else 'ns'}</div>
                    </div>
                    
                    <div class="variable-stats">
                        <div class="stat-item">
                            <span class="stat-label">R² (Effect Size):</span>
                            <div class="r-squared-bar">
                                <div class="r-squared-fill {sig_class}" style="width: {r_squared_percent:.1f}%"></div>
                                <span class="r-squared-value">{row['r_squared']:.3f} ({r_squared_percent:.1f}%)</span>
                            </div>
                        </div>
                        
                        <div class="stat-grid">
                            <div class="stat-item">
                                <span class="stat-label">p-value:</span>
                                <span class="stat-value">{row['p_value']:.4f}</span>
                            </div>
                            <div class="stat-item">
                                <span class="stat-label">Type:</span>
                                <span class="stat-value">{row['variable_type']}</span>
                            </div>
                            <div class="stat-item">
                                <span class="stat-label">Valid samples:</span>
                                <span class="stat-value">{row['n_samples']}</span>
                            </div>
                            <div class="stat-item">
                                <span class="stat-label">Groups:</span>
                                <span class="stat-value">{row['n_groups']}</span>
                            </div>
                        </div>
                    </div>
                    
                    <div class="variable-interpretation">
                        {self._get_permanova_interpretation(row)}
                    </div>
                </div>
            """
        
        permanova_html += """
                </div>
            </div>
            
            <div class="permanova-recommendations">
                <h4>💡 Recommendations</h4>
        """
        
        # Add recommendations based on top variables
        if not sorted_results.empty:
            top_var = sorted_results.iloc[0]
            if top_var['p_value'] < 0.05:
                permanova_html += f"""
                    <div class="recommendation-item significant">
                        <strong>🎯 Primary grouping variable:</strong> <em>{top_var['variable']}</em> explains {top_var['r_squared']:.1%} of sample composition differences (p = {top_var['p_value']:.4f}).
                        <br><strong>Recommendation:</strong> Consider grouping samples by this variable for co-assembly.
                    </div>
                """
            else:
                permanova_html += """
                    <div class="recommendation-item not-significant">
                        <strong>⚠️ No strong metadata associations found.</strong>
                        <br><strong>Recommendation:</strong> Consider individual sample assembly or use similarity-based grouping instead.
                    </div>
                """
        
        permanova_html += """
            </div>
        </div>
        """
        
        return permanova_html
    
    def _get_permanova_interpretation(self, row) -> str:
        """Generate interpretation text for a PERMANOVA result."""
        
        r_squared = row['r_squared']
        p_value = row['p_value']
        variable = row['variable']
        
        # Effect size interpretation
        if r_squared >= 0.3:
            effect_size = "large effect"
        elif r_squared >= 0.1:
            effect_size = "medium effect" 
        elif r_squared >= 0.05:
            effect_size = "small effect"
        else:
            effect_size = "very small effect"
            
        # Significance interpretation
        if p_value < 0.05:
            significance_text = "statistically significant"
            action = f"Strong evidence that {variable} influences sample composition."
        else:
            significance_text = "not statistically significant"
            action = f"Insufficient evidence that {variable} influences sample composition."
            
        return f"""
        <p><strong>Interpretation:</strong> This variable has a <em>{effect_size}</em> on sample composition and is <em>{significance_text}</em>. {action}</p>
        """
    
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
        strategy_explanation = self.create_assembly_strategy_explanation()
        permanova_section = self.create_permanova_section()
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
            strategy_explanation=strategy_explanation,
            permanova_section=permanova_section,
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
            self.report_data.get('metadata'),
            self.report_data.get('analyzed_variables')
        )
        
        # Create PCA if we have k-mer data
        pca_html = ""
        if 'kmer_data' in self.report_data and self.report_data['kmer_data']:
            try:
                from sklearn.decomposition import PCA
                print(f"🔬 Found k-mer data, generating PCA plot...")
                logging.info(f"Found k-mer data, generating PCA plot...")
                
                # Extract k-mer profiles for PCA
                profiles = self.report_data['kmer_data'].get('profiles', {})
                print(f"📊 K-mer profiles: {len(profiles)} samples found")
                print(f"🔍 Profile keys: {list(profiles.keys())[:3] if profiles else 'None'}...")
                logging.info(f"K-mer profiles: {len(profiles)} samples found")
                if profiles:
                    # Check profile structure
                    first_sample = list(profiles.keys())[0]
                    first_profile = profiles[first_sample]
                    print(f"📋 First sample '{first_sample}' has {len(first_profile)} k-mers")
                    print(f"🧬 Sample k-mers: {list(first_profile.keys())[:3]}...")
                    print(f"📊 Sample values: {list(first_profile.values())[:3]}...")
                    # Convert to matrix
                    sample_names = list(profiles.keys())
                    all_kmers = set()
                    for profile in profiles.values():
                        all_kmers.update(profile.keys())
                    
                    print(f"🧮 Creating matrix: {len(sample_names)} samples × {len(all_kmers)} k-mers")
                    kmer_matrix = np.zeros((len(sample_names), len(all_kmers)))
                    kmer_list = list(all_kmers)
                    
                    for i, sample in enumerate(sample_names):
                        for j, kmer in enumerate(kmer_list):
                            kmer_matrix[i, j] = profiles[sample].get(kmer, 0)
                    
                    print(f"✅ Matrix created: shape {kmer_matrix.shape}")
                    print(f"📊 Matrix stats: min={kmer_matrix.min():.3f}, max={kmer_matrix.max():.3f}, mean={kmer_matrix.mean():.3f}")
                    
                    # Perform PCA
                    print(f"🔄 Computing PCA...")
                    pca = PCA(n_components=2)
                    pca_result = pca.fit_transform(kmer_matrix)
                    print(f"✅ PCA completed: {pca_result.shape}")
                    print(f"📊 Explained variance: {pca.explained_variance_ratio_}")
                    
                    print(f"🎨 Creating enhanced plot...")
                    
                    # Create enhanced plot with multiple dimensionality reduction methods
                    fig = self._create_enhanced_kmer_plot(kmer_matrix, sample_names, pca, pca_result)
                    
                    # Convert to HTML div without full page structure
                    pca_html = fig.to_html(include_plotlyjs='cdn', div_id="enhanced-plot", config={'displayModeBar': True})
                    print(f"✅ Enhanced plot HTML generated successfully")
                else:
                    print(f"❌ No k-mer profiles found in data")
                    pca_html = "<p>No k-mer profiles available for PCA</p>"
                    
            except Exception as e:
                print(f"❌ PCA plot creation failed: {e}")
                print(f"📋 Error type: {type(e).__name__}")
                import traceback
                print(f"🔍 Full traceback:")
                traceback.print_exc()
                logging.warning(f"Could not create PCA plot: {e}")
                pca_html = f"<p>PCA visualization not available - Error: {e}</p>"
        else:
            print(f"⚠️  No k-mer data available - using fallback MDS visualization")
            # Fallback: Use distance matrix for PCA-like visualization
            try:
                from sklearn.decomposition import PCA
                from sklearn.manifold import MDS
                logging.info("No k-mer data available, using distance matrix for MDS visualization...")
                
                # Use MDS (multidimensional scaling) on distance matrix
                mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
                mds_result = mds.fit_transform(self.report_data['distance_matrix'])
                
                import plotly.express as px
                
                pca_df = pd.DataFrame({
                    'MDS1': mds_result[:, 0],
                    'MDS2': mds_result[:, 1],
                    'sample_id': self.report_data['sample_names']
                })
                
                # Add metadata if available
                if self.report_data.get('metadata') is not None:
                    # Use same metadata processing logic as enhanced plot path
                    sample_id_column = self.report_data.get('sample_id_column', 'sample_id')
                    metadata_orig = self.report_data['metadata']
                    
                    # Handle metadata preparation consistently
                    if sample_id_column == metadata_orig.index.name:
                        metadata_for_merge = metadata_orig.reset_index()
                        if sample_id_column != 'sample_id':
                            metadata_for_merge['sample_id'] = metadata_for_merge[sample_id_column].astype(str)
                    elif sample_id_column in metadata_orig.columns:
                        metadata_for_merge = metadata_orig.reset_index(drop=True)
                        if sample_id_column != 'sample_id':
                            metadata_for_merge['sample_id'] = metadata_for_merge[sample_id_column].astype(str)
                        else:
                            metadata_for_merge['sample_id'] = metadata_for_merge['sample_id'].astype(str)
                    else:
                        metadata_for_merge = metadata_orig.reset_index()
                        metadata_for_merge['sample_id'] = metadata_for_merge.index.astype(str)
                    
                    # Fix data type mismatch by converting pca_df to strings (metadata already converted)
                    pca_df['sample_id'] = pca_df['sample_id'].astype(str)
                    # metadata_for_merge['sample_id'] already converted to string above
                    
                    # Check for and handle duplicate sample IDs in metadata
                    if metadata_for_merge['sample_id'].duplicated().any():
                        duplicates = metadata_for_merge['sample_id'][metadata_for_merge['sample_id'].duplicated(keep=False)]
                        logging.warning(f"Found {len(duplicates)} duplicate sample IDs in metadata for merge: {list(duplicates.unique())}")
                        
                        # Remove duplicates, keeping the first occurrence
                        metadata_for_merge = metadata_for_merge[~metadata_for_merge['sample_id'].duplicated(keep='first')]
                        logging.info(f"Removed duplicates for merge, kept first occurrence for each sample ID")
                    
                    pca_df = pca_df.merge(metadata_for_merge, on='sample_id', how='left')
                    
                    # Check matching for fallback case too
                    metadata_cols_for_checking = [col for col in pca_df.columns if col not in ['sample_id', 'MDS1', 'MDS2']]
                    if metadata_cols_for_checking:
                        first_metadata_col = metadata_cols_for_checking[0]
                        successful_matches = pca_df[first_metadata_col].notna().sum()
                        print(f"🔍 Fallback MDS merge results: {successful_matches}/{len(pca_df)} samples matched")
                
                # Get metadata columns for coloring options
                metadata_cols = []
                if self.report_data.get('metadata') is not None:
                    analyzed_variables = self.report_data.get('analyzed_variables')
                    
                    if analyzed_variables:
                        # Use user-specified variables from --variables flag
                        for col in analyzed_variables:
                            if col in pca_df.columns and col not in ['sample_id', 'MDS1', 'MDS2']:
                                n_unique = pca_df[col].nunique()
                                non_null = pca_df[col].notna().sum()
                                if n_unique > 1 and non_null > 0:  # Relaxed criteria for user-specified variables
                                    metadata_cols.append(col)
                    else:
                        # Fallback to auto-detection when no variables specified
                        for col in pca_df.columns:
                            if col not in ['sample_id', 'MDS1', 'MDS2']:
                                n_unique = pca_df[col].nunique()
                                non_null = pca_df[col].notna().sum()
                                if n_unique > 1 and n_unique <= 20 and non_null > 0:
                                    metadata_cols.append(col)
                
                # Create base plot with first metadata variable as default color
                color_col = metadata_cols[0] if metadata_cols else None
                hover_cols = ['sample_id'] + metadata_cols if metadata_cols else ['sample_id']
                
                fig = px.scatter(
                    pca_df, 
                    x='MDS1', 
                    y='MDS2',
                    color=color_col,
                    hover_data=hover_cols,
                    title="Multi-Dimensional Scaling (Distance-based)",
                    labels={
                        'MDS1': 'MDS Dimension 1',
                        'MDS2': 'MDS Dimension 2'
                    }
                )
                
                # Add color selection dropdown if metadata available
                if metadata_cols:
                    buttons = []
                    for col in metadata_cols:
                        buttons.append(
                            dict(
                                label=col.replace('_', ' ').title(),
                                method="restyle",
                                args=[{"marker.color": pca_df[col]}]
                            )
                        )
                    
                    # Add no coloring option
                    buttons.append(
                        dict(
                            label="No Coloring",
                            method="restyle", 
                            args=[{"marker.color": "blue"}]
                        )
                    )
                    
                    fig.update_layout(
                        updatemenus=[
                            dict(
                                buttons=buttons,
                                direction="down",
                                showactive=True,
                                x=0.1,
                                y=1.15,
                                xanchor="left",
                                yanchor="top"
                            )
                        ],
                        annotations=[
                            dict(
                                text="Color by:",
                                x=0.05, y=1.18,
                                xref="paper", yref="paper",
                                align="left",
                                showarrow=False
                            )
                        ]
                    )
                
                fig.update_traces(marker=dict(size=10, opacity=0.8))
                fig.update_layout(height=500, template='plotly_white')
                
                # Convert to HTML div without full page structure
                pca_html = fig.to_html(include_plotlyjs='cdn', div_id="mds-plot", config={'displayModeBar': True})
                
            except Exception as e:
                logging.warning(f"Could not create MDS plot: {e}")
                pca_html = "<p>No visualization data available</p>"
        
        # Combine visualizations (removed heatmap)
        viz_html = f"""
        <div class="visualizations-container">
            <div class="viz-section">
                <h4>📊 Interactive Dimensionality Analysis</h4>
                <div class="viz-content">
                    {pca_html}
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
        if 'assembly_recommendation' in self.report_data and self.report_data['assembly_recommendation'] is not None:
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
            
            .interpretation-grid {
                grid-template-columns: 1fr;
            }
            
            .threshold-explanation {
                padding: 15px;
            }
            
            .guide-section {
                padding: 10px;
            }
        }
        
        /* PERMANOVA Section Styles */
        .permanova-section {
            margin: 20px 0;
        }
        
        .permanova-explanation {
            background: #f8f9fa;
            border-radius: 8px;
            padding: 20px;
            margin-bottom: 20px;
            border-left: 4px solid #667eea;
        }
        
        .permanova-guide {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }
        
        .guide-item {
            background: white;
            padding: 10px;
            border-radius: 6px;
            border: 1px solid #e9ecef;
        }
        
        .variables-container {
            display: grid;
            gap: 20px;
        }
        
        .variable-card {
            background: white;
            border-radius: 10px;
            padding: 20px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
            border-left: 4px solid #cccccc;
        }
        
        .variable-card.highly-significant {
            border-left-color: #dc3545;
        }
        
        .variable-card.very-significant {
            border-left-color: #fd7e14;
        }
        
        .variable-card.significant {
            border-left-color: #ffc107;
        }
        
        .variable-card.marginally-significant {
            border-left-color: #20c997;
        }
        
        .variable-card.not-significant {
            border-left-color: #6c757d;
        }
        
        .variable-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 15px;
        }
        
        .variable-header h5 {
            margin: 0;
            color: #495057;
        }
        
        .significance-badge {
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 0.9em;
            color: white;
        }
        
        .significance-badge.highly-significant {
            background-color: #dc3545;
        }
        
        .significance-badge.very-significant {
            background-color: #fd7e14;
        }
        
        .significance-badge.significant {
            background-color: #ffc107;
            color: #212529;
        }
        
        .significance-badge.marginally-significant {
            background-color: #20c997;
        }
        
        .significance-badge.not-significant {
            background-color: #6c757d;
        }
        
        .variable-stats {
            margin-bottom: 15px;
        }
        
        .r-squared-bar {
            position: relative;
            background: #e9ecef;
            height: 25px;
            border-radius: 12px;
            overflow: hidden;
            margin: 8px 0;
        }
        
        .r-squared-fill {
            height: 100%;
            border-radius: 12px;
            transition: width 0.3s ease;
        }
        
        .r-squared-fill.highly-significant {
            background: linear-gradient(90deg, #dc3545, #c82333);
        }
        
        .r-squared-fill.very-significant {
            background: linear-gradient(90deg, #fd7e14, #e8640f);
        }
        
        .r-squared-fill.significant {
            background: linear-gradient(90deg, #ffc107, #e0a800);
        }
        
        .r-squared-fill.marginally-significant {
            background: linear-gradient(90deg, #20c997, #1aa179);
        }
        
        .r-squared-fill.not-significant {
            background: linear-gradient(90deg, #6c757d, #5a6268);
        }
        
        .r-squared-value {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            color: white;
            font-weight: bold;
            font-size: 0.9em;
            text-shadow: 1px 1px 2px rgba(0, 0, 0, 0.3);
        }
        
        .stat-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
            gap: 10px;
            margin-top: 10px;
        }
        
        .stat-item {
            display: flex;
            flex-direction: column;
        }
        
        .stat-label {
            font-weight: bold;
            color: #6c757d;
            font-size: 0.9em;
            margin-bottom: 2px;
        }
        
        .stat-value {
            color: #495057;
        }
        
        .variable-interpretation {
            background: #f8f9fa;
            padding: 12px;
            border-radius: 6px;
            margin-top: 15px;
            border: 1px solid #e9ecef;
        }
        
        .variable-interpretation p {
            margin: 0;
            font-size: 0.95em;
            line-height: 1.4;
        }
        
        .permanova-recommendations {
            background: #fff3cd;
            border: 1px solid #ffeaa7;
            border-radius: 8px;
            padding: 20px;
            margin-top: 20px;
        }
        
        .permanova-recommendations h4 {
            color: #856404;
            margin-bottom: 15px;
        }
        
        .recommendation-item {
            padding: 15px;
            border-radius: 6px;
            margin-bottom: 10px;
        }
        
        .recommendation-item.significant {
            background: #d4edda;
            border: 1px solid #c3e6cb;
            color: #155724;
        }
        
        .recommendation-item.not-significant {
            background: #f8d7da;
            border: 1px solid #f5c6cb;
            color: #721c24;
        }
        
        @media (max-width: 768px) {
            .permanova-guide {
                grid-template-columns: 1fr;
            }
            
            .stat-grid {
                grid-template-columns: 1fr 1fr;
            }
            
            .variable-header {
                flex-direction: column;
                align-items: flex-start;
            }
            
            .significance-badge {
                margin-top: 5px;
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
        
        <!-- PERMANOVA Analysis -->
        <div class="section">
            <h2>🧬 Metadata Analysis (PERMANOVA)</h2>
            <p>Statistical analysis of which metadata variables significantly explain differences in sample composition.</p>
            {{ permanova_section|safe }}
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
    
    def _try_match_sample_names(self, fastq_names, metadata_ids):
        """Try to create a mapping between FASTQ names and metadata IDs using common transformations."""
        from pathlib import Path
        print(f"🔧 Attempting smart name matching...")
        
        # Try different transformations
        transformations = {
            'exact': lambda x: x,
            'remove_last_underscore': lambda x: '_'.join(x.split('_')[:-1]) if '_' in x else x,  # sample_001_hr -> sample_001, sample_clean -> sample
            'remove_suffix': lambda x: x.split('_')[0].split('.')[0],  # SRR123_1.fastq -> SRR123
            'extract_srr': lambda x: x if x.startswith('SRR') else ('SRR' + x if x.isdigit() else x),  # 123 -> SRR123
            'remove_srr_prefix': lambda x: x[3:] if x.startswith('SRR') else x,  # SRR123 -> 123
            'basename_only': lambda x: Path(x).stem.split('_')[0]  # /path/to/SRR123_1.fastq.gz -> SRR123
        }
        
        best_matches = 0
        best_transform = None
        best_mapping = {}
        
        fastq_set = set(fastq_names)
        metadata_set = set(metadata_ids)
        
        # Smart per-sample matching: preserve existing matches, transform others
        final_mapping = {}
        
        for transform_name, transform_func in transformations.items():
            try:
                current_mapping = {}
                current_matches = 0
                
                for orig_name in fastq_names:
                    # If this is exact matching or name doesn't already have an exact match
                    if transform_name == 'exact' or orig_name not in metadata_set:
                        transformed = transform_func(orig_name)
                        current_mapping[orig_name] = transformed
                        if transformed in metadata_set:
                            current_matches += 1
                    else:
                        # Preserve existing exact matches
                        current_mapping[orig_name] = orig_name
                        current_matches += 1
                
                print(f"   {transform_name}: {current_matches}/{len(fastq_names)} matches")
                
                if current_matches > best_matches:
                    best_matches = current_matches
                    best_transform = transform_name
                    best_mapping = current_mapping
                    
            except Exception as e:
                print(f"   {transform_name}: failed ({e})")
        
        if best_matches > 0:
            print(f"✅ Best transformation: '{best_transform}' with {best_matches} matches")
            return best_mapping, best_transform
        else:
            print(f"❌ No transformations yielded matches")
            return {}, None

    def _create_enhanced_kmer_plot(self, kmer_matrix, sample_names, pca, pca_result):
        """Create enhanced plot with method switching and metadata coloring."""
        import plotly.graph_objects as go
        from sklearn.manifold import TSNE
        
        # Compute multiple dimensionality reduction methods
        projections = {'pca': (pca_result, pca)}
        
        # Compute t-SNE
        try:
            perplexity = min(30, max(5, len(sample_names) // 4))
            tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42, max_iter=1000)
            tsne_result = tsne.fit_transform(kmer_matrix)
            projections['tsne'] = (tsne_result, tsne)
        except Exception as e:
            logging.warning(f"t-SNE computation failed: {e}")
        
        # Compute UMAP if available
        try:
            import umap
            n_neighbors = min(15, max(2, len(sample_names) // 3))
            umap_reducer = umap.UMAP(n_components=2, n_neighbors=n_neighbors, random_state=42)
            umap_result = umap_reducer.fit_transform(kmer_matrix)
            projections['umap'] = (umap_result, umap_reducer)
        except ImportError:
            logging.info("UMAP not available, skipping")
        except Exception as e:
            logging.warning(f"UMAP computation failed: {e}")
        
        # Create base DataFrame
        plot_data = {'sample_id': sample_names}
        method_models = {}
        
        for method, (projection, model) in projections.items():
            plot_data[f'{method}_x'] = projection[:, 0]
            plot_data[f'{method}_y'] = projection[:, 1]
            method_models[method] = model
        
        plot_df = pd.DataFrame(plot_data)
        
        # Add metadata if available
        metadata_cols = []
        if self.report_data.get('metadata') is not None:
            print(f"🔍 Processing metadata for interactive plot...")
            print(f"📊 Found metadata with {len(self.report_data['metadata'])} rows and columns: {list(self.report_data['metadata'].columns)}")
            logging.info(f"Found metadata with {len(self.report_data['metadata'])} rows and columns: {list(self.report_data['metadata'].columns)}")
            
            # Get sample ID column name
            sample_id_column = self.report_data.get('sample_id_column', 'sample_id')
            print(f"🔍 Using '{sample_id_column}' column for sample matching")
            
            # Handle metadata preparation more carefully
            metadata_orig = self.report_data['metadata']
            
            # CRITICAL DEBUG: Check what's actually in the original metadata
            print(f"🔍 CRITICAL DEBUG - Original metadata shape: {metadata_orig.shape}")
            print(f"🔍 CRITICAL DEBUG - Index name: {metadata_orig.index.name}")
            print(f"🔍 CRITICAL DEBUG - First few index values: {metadata_orig.index.tolist()[:3]}")
            debug_cols = ['FC_categories', 'current_medications', 'patient_ID']
            for col in debug_cols:
                if col in metadata_orig.columns:
                    non_null = metadata_orig[col].notna().sum()
                    unique_vals = metadata_orig[col].nunique()
                    print(f"🔍 CRITICAL DEBUG - {col}: {non_null} non-null, {unique_vals} unique")
                    if non_null > 0:
                        sample_vals = metadata_orig[col].dropna().unique()[:3]
                        print(f"   Sample values: {sample_vals.tolist()}")
                else:
                    print(f"🔍 CRITICAL DEBUG - {col}: NOT FOUND in original metadata")
            
            # Check if the sample_id_column is the index
            if sample_id_column == metadata_orig.index.name:
                print(f"📋 '{sample_id_column}' is the DataFrame index")
                metadata_for_merge = metadata_orig.reset_index()
                # After reset_index, the index becomes a column
                if sample_id_column != 'sample_id':
                    metadata_for_merge['sample_id'] = metadata_for_merge[sample_id_column].astype(str)  # Convert to string immediately!
                    print(f"✅ Using index '{sample_id_column}' as sample_id (converted to string)")
            elif sample_id_column in metadata_orig.columns:
                print(f"📋 '{sample_id_column}' is a regular column")
                metadata_for_merge = metadata_orig.reset_index(drop=True)  # Don't add index as column
                if sample_id_column != 'sample_id':
                    metadata_for_merge['sample_id'] = metadata_for_merge[sample_id_column].astype(str)  # Convert to string immediately!
                    print(f"✅ Using column '{sample_id_column}' as sample_id (converted to string)")
                else:
                    # Even if it's already called sample_id, make sure it's a string
                    metadata_for_merge['sample_id'] = metadata_for_merge['sample_id'].astype(str)
                    print(f"✅ Using existing 'sample_id' column (converted to string)")
            else:
                print(f"❌ Column '{sample_id_column}' not found in metadata!")
                print(f"📋 Available columns: {list(metadata_orig.columns)}")
                print(f"📋 Index name: {metadata_orig.index.name}")
                # Fallback to index
                metadata_for_merge = metadata_orig.reset_index()
                metadata_for_merge['sample_id'] = metadata_for_merge.index
                print(f"🔧 Falling back to row index as 'sample_id'")
            
            print(f"🔗 Sample names for merging: {sample_names[:3]}...")
            print(f"🔗 Metadata sample IDs: {metadata_for_merge['sample_id'].tolist()[:3]}...")
            print(f"📊 Sample names type: {type(sample_names[0]) if sample_names else 'None'}")
            print(f"📊 Metadata sample_id type: {metadata_for_merge['sample_id'].dtype}")
            
            # Fix data type mismatch by converting plot_df to strings (metadata already converted)
            plot_df['sample_id'] = plot_df['sample_id'].astype(str)
            # metadata_for_merge['sample_id'] already converted to string above
            
            # Check for and handle duplicate sample IDs in metadata
            if metadata_for_merge['sample_id'].duplicated().any():
                duplicates = metadata_for_merge['sample_id'][metadata_for_merge['sample_id'].duplicated(keep=False)]
                print(f"⚠️  Found {len(duplicates)} duplicate sample IDs in metadata: {list(duplicates.unique())}")
                logging.warning(f"Found {len(duplicates)} duplicate sample IDs in metadata for enhanced plot: {list(duplicates.unique())}")
                
                # Remove duplicates, keeping the first occurrence
                metadata_for_merge = metadata_for_merge[~metadata_for_merge['sample_id'].duplicated(keep='first')]
                print(f"✅ Removed duplicates, kept first occurrence for each sample ID")
                logging.info(f"Removed duplicates for enhanced plot, kept first occurrence for each sample ID")
            
            print(f"✅ Converted both to strings for merging")
            logging.info(f"Sample names for merging: {sample_names[:3]}...")
            logging.info(f"Metadata sample IDs: {metadata_for_merge['sample_id'].tolist()[:3]}...")
            
            # DEBUG: Check a few specific variables before merge
            debug_vars = ['FC_categories', 'current_medications', 'patient_ID']
            for var in debug_vars:
                if var in metadata_for_merge.columns:
                    non_null_count = metadata_for_merge[var].notna().sum()
                    unique_count = metadata_for_merge[var].nunique()
                    print(f"🔍 PRE-MERGE: {var} has {non_null_count} non-null, {unique_count} unique values")
                    if non_null_count > 0:
                        print(f"   Sample values: {metadata_for_merge[var].dropna().unique()[:3].tolist()}")
                else:
                    print(f"🔍 PRE-MERGE: {var} NOT FOUND in metadata_for_merge")
            
            plot_df = plot_df.merge(metadata_for_merge, on='sample_id', how='left')
            print(f"✅ After merge, plot_df columns: {list(plot_df.columns)}")
            
            # Check how many samples actually matched
            metadata_cols_for_checking = [col for col in plot_df.columns if col not in ['sample_id'] and not col.endswith(('_x', '_y'))]
            if metadata_cols_for_checking:
                first_metadata_col = metadata_cols_for_checking[0]
                successful_matches = plot_df[first_metadata_col].notna().sum()
                print(f"🔍 Merge results: {successful_matches}/{len(plot_df)} samples successfully matched")
                
                if successful_matches == 0:
                    print(f"❌ NO SAMPLES MATCHED! Checking name format differences...")
                    print(f"📋 FASTQ sample names: {plot_df['sample_id'].tolist()[:5]}")
                    print(f"📋 Metadata sample IDs: {metadata_for_merge['sample_id'].tolist()[:5]}")
                    
                    # Try smart name matching
                    name_mapping, best_transform = self._try_match_sample_names(
                        plot_df['sample_id'].tolist(),
                        metadata_for_merge['sample_id'].tolist()
                    )
                    
                    if name_mapping and best_transform:
                        print(f"🔧 Applying '{best_transform}' transformation and re-merging...")
                        
                        # Clean plot_df to just the essential columns before re-merge
                        essential_cols = ['sample_id', 'pca_x', 'pca_y', 'tsne_x', 'tsne_y', 'umap_x', 'umap_y']
                        plot_df_clean = plot_df[essential_cols].copy()
                        
                        # Apply the best transformation
                        plot_df_clean['sample_id_transformed'] = plot_df_clean['sample_id'].map(name_mapping)
                        plot_df_clean = plot_df_clean.drop('sample_id', axis=1).rename(columns={'sample_id_transformed': 'sample_id'})
                        
                        # Re-merge with transformed names (no column conflicts now)
                        plot_df = plot_df_clean.merge(metadata_for_merge, on='sample_id', how='left')
                        
                        # Check success - recalculate metadata columns after transformation
                        updated_metadata_cols = [col for col in plot_df.columns if col not in ['sample_id'] and not col.endswith(('_x', '_y'))]
                        if updated_metadata_cols:
                            updated_first_col = updated_metadata_cols[0]
                            new_matches = plot_df[updated_first_col].notna().sum()
                            print(f"✅ After transformation: {new_matches}/{len(plot_df)} samples matched!")
                        else:
                            print(f"⚠️ No metadata columns available after transformation")
                    else:
                        print(f"❌ Smart matching failed - no suitable transformations found")
                        print(f"🔍 Name format analysis:")
                        fastq_names = set(plot_df['sample_id'].tolist())
                        metadata_ids = set(metadata_for_merge['sample_id'].tolist())
                        print(f"   FASTQ names example: {list(fastq_names)[:3]}")
                        print(f"   Metadata IDs example: {list(metadata_ids)[:3]}")
                else:
                    print(f"✅ {successful_matches} samples matched successfully")
                    print(f"📋 Sample with metadata: {plot_df[plot_df[first_metadata_col].notna()]['sample_id'].tolist()[:3]}")
            
            logging.info(f"After merge, plot_df columns: {list(plot_df.columns)}")
            
            # Get metadata columns for coloring options
            analyzed_variables = self.report_data.get('analyzed_variables')
            
            if analyzed_variables:
                # Use user-specified variables from --variables flag
                print(f"🎯 Using user-specified variables: {analyzed_variables}")
                print(f"🔍 Available columns in plot_df after merge: {list(plot_df.columns)}")
                for col in analyzed_variables:
                    print(f"🧪 Checking variable '{col}'...")
                    if col in plot_df.columns:
                        if col not in ['sample_id'] and not col.endswith(('_x', '_y')):
                            n_unique = plot_df[col].nunique()
                            non_null = plot_df[col].notna().sum()
                            print(f"🎨 Column {col}: {n_unique} unique values, {non_null} non-null values")
                            if non_null > 0:
                                print(f"   📋 Sample values: {plot_df[col].dropna().unique()[:5].tolist()}")
                            else:
                                print(f"   ❌ All values are NaN! This suggests merge failed for this column")
                            logging.info(f"Column {col}: {n_unique} unique values, {non_null} non-null values")
                            if n_unique > 1 and non_null > 0:  # Relaxed criteria for user-specified variables
                                metadata_cols.append(col)
                                print(f"✅ Added {col} to coloring options")
                            else:
                                print(f"❌ Skipped {col}: n_unique={n_unique}, non_null={non_null}")
                        else:
                            print(f"   ⏭️  Skipped {col}: excluded column type")
                    else:
                        print(f"   ❌ Column '{col}' not found in plot_df!")
            else:
                # Fallback to auto-detection when no variables specified
                print(f"🔍 Auto-detecting metadata columns for coloring")
                for col in plot_df.columns:
                    if col not in ['sample_id'] and not col.endswith(('_x', '_y')):
                        n_unique = plot_df[col].nunique()
                        non_null = plot_df[col].notna().sum()
                        print(f"🎨 Column {col}: {n_unique} unique values, {non_null} non-null values")
                        logging.info(f"Column {col}: {n_unique} unique values, {non_null} non-null values")
                        if n_unique > 1 and n_unique <= 20 and non_null > 0:
                            metadata_cols.append(col)
                            print(f"✅ Added {col} to coloring options")
            
            print(f"🎨 Selected metadata columns for coloring: {metadata_cols}")
            logging.info(f"Selected metadata columns for coloring: {metadata_cols}")
        else:
            print(f"❌ No metadata available for interactive plot")
        
        # Create figure with traces for each method
        fig = go.Figure()
        
        default_method = list(projections.keys())[0]
        default_color = metadata_cols[0] if metadata_cols else None
        
        # Add traces for each method
        for i, (method_name, (projection, model)) in enumerate(projections.items()):
            visible = method_name == default_method
            
            # Determine axis labels
            if method_name == 'pca':
                x_label = f'PC1 ({model.explained_variance_ratio_[0]:.1%} variance)'
                y_label = f'PC2 ({model.explained_variance_ratio_[1]:.1%} variance)'
            else:
                x_label = f'{method_name.upper()} 1'
                y_label = f'{method_name.upper()} 2'
            
            # Prepare hover data
            hover_text = []
            for idx in range(len(sample_names)):
                hover_info = f"<b>{sample_names[idx]}</b><br>"
                hover_info += f"{x_label}: {projection[idx, 0]:.3f}<br>"
                hover_info += f"{y_label}: {projection[idx, 1]:.3f}"
                
                # Add metadata to hover
                if len(metadata_cols) > 0:
                    hover_info += "<br>--- Metadata ---"
                    for col in metadata_cols:
                        if col in plot_df.columns:
                            val = plot_df.iloc[idx][col]
                            if pd.notna(val):  # Only show non-null values
                                hover_info += f"<br>{col.replace('_', ' ').title()}: {val}"
                
                hover_text.append(hover_info)
            
            # Determine color array
            if default_color and default_color in plot_df.columns:
                # Check if categorical or continuous
                is_categorical = (plot_df[default_color].dtype == 'object' or 
                                pd.api.types.is_categorical_dtype(plot_df[default_color]) or
                                plot_df[default_color].nunique() <= 10)  # Treat <= 10 unique values as categorical
                
                if is_categorical:
                    # For categorical data, use discrete color mapping
                    import plotly.colors as pc
                    categories = plot_df[default_color].dropna().unique()  # Remove NaN values
                    color_discrete_map = {cat: pc.qualitative.Set3[i % len(pc.qualitative.Set3)] 
                                        for i, cat in enumerate(categories)}
                    # Add a color for missing values
                    color_discrete_map[pd.NA] = '#cccccc'  # Light gray for missing
                    color_discrete_map[None] = '#cccccc'
                    # Handle NaN values properly
                    color_data = []
                    for val in plot_df[default_color]:
                        if pd.isna(val):
                            color_data.append('#cccccc')  # Gray for missing values
                        else:
                            color_data.append(color_discrete_map[val])
                    colorscale = None
                    colorbar = None
                    showscale = False  # Don't show colorbar for categorical
                else:
                    # For numeric data, use continuous colorscale
                    color_data = plot_df[default_color].tolist()
                    colorscale = 'viridis'
                    colorbar = dict(title=default_color.replace('_', ' ').title())
                    showscale = True
            else:
                color_data = 'blue'
                colorscale = None
                colorbar = None
                showscale = False
            
            fig.add_trace(go.Scatter(
                x=projection[:, 0],
                y=projection[:, 1],
                mode='markers',
                marker=dict(
                    size=10,
                    opacity=0.8,
                    line=dict(width=1, color='DarkSlateGrey'),
                    color=color_data,
                    colorscale=colorscale,
                    colorbar=colorbar,
                    showscale=showscale
                ),
                name=method_name.upper(),
                visible=visible,
                hovertemplate="%{text}<extra></extra>",
                text=hover_text
            ))
        
        # Create method switching buttons
        method_buttons = []
        for i, (method_name, (projection, model)) in enumerate(projections.items()):
            if method_name == 'pca':
                x_title = f'PC1 ({model.explained_variance_ratio_[0]:.1%} variance)'
                y_title = f'PC2 ({model.explained_variance_ratio_[1]:.1%} variance)'
            else:
                x_title = f'{method_name.upper()} 1'
                y_title = f'{method_name.upper()} 2'
            
            visibility = [False] * len(projections)
            visibility[i] = True
            
            method_buttons.append(
                dict(
                    label=method_name.upper(),
                    method="update",
                    args=[
                        {"visible": visibility},
                        {"xaxis.title": x_title, "yaxis.title": y_title}
                    ]
                )
            )
        
        # Create color buttons if metadata available
        color_buttons = []
        if metadata_cols:
            for col in metadata_cols:
                # Check if categorical or continuous
                is_categorical = (plot_df[col].dtype == 'object' or 
                                pd.api.types.is_categorical_dtype(plot_df[col]) or
                                plot_df[col].nunique() <= 10)  # Treat <= 10 unique values as categorical
                
                if is_categorical:
                    # For categorical data, use discrete color mapping
                    import plotly.colors as pc
                    categories = plot_df[col].dropna().unique()  # Remove NaN values
                    color_discrete_map = {cat: pc.qualitative.Set3[i % len(pc.qualitative.Set3)] 
                                        for i, cat in enumerate(categories)}
                    # Add a color for missing values
                    color_discrete_map[pd.NA] = '#cccccc'  # Light gray for missing
                    color_discrete_map[None] = '#cccccc'
                    # Handle NaN values properly
                    color_data = []
                    for val in plot_df[col]:
                        if pd.isna(val):
                            color_data.append('#cccccc')  # Gray for missing values
                        else:
                            color_data.append(color_discrete_map[val])
                    colorscale = None
                    showscale_setting = False
                else:
                    # For numeric data, use continuous colorscale
                    color_data = plot_df[col].tolist()
                    colorscale = "viridis"
                    showscale_setting = True
                
                restyle_args = {
                    "marker.color": [color_data] * len(projections),
                    "marker.colorbar.title.text": col.replace('_', ' ').title(),
                    "marker.colorscale": colorscale,
                    "marker.showscale": showscale_setting
                }
                
                color_buttons.append(
                    dict(
                        label=col.replace('_', ' ').title(),
                        method="restyle",
                        args=[restyle_args]
                    )
                )
            
            # Add no coloring option
            color_buttons.append(
                dict(
                    label="No Coloring",
                    method="restyle",
                    args=[{
                        "marker.color": ["blue"] * len(projections),
                        "marker.colorbar.title.text": "",
                        "marker.colorscale": None
                    }]
                )
            )
        
        # Create layout with dropdowns
        updatemenus = []
        annotations = []
        
        # Method selector
        updatemenus.append(
            dict(
                buttons=method_buttons,
                direction="down",
                showactive=True,
                x=0.1,
                y=1.15,
                xanchor="left",
                yanchor="top"
            )
        )
        annotations.append(
            dict(
                text="Method:",
                x=0.05, y=1.18,
                xref="paper", yref="paper",
                align="left",
                showarrow=False
            )
        )
        
        # Color selector (if metadata available)
        if color_buttons:
            logging.info(f"Creating color dropdown with {len(color_buttons)} options")
            updatemenus.append(
                dict(
                    buttons=color_buttons,
                    direction="down", 
                    showactive=True,
                    x=0.4,
                    y=1.15,
                    xanchor="left",
                    yanchor="top"
                )
            )
            annotations.append(
                dict(
                    text="Color by:",
                    x=0.35, y=1.18,
                    xref="paper", yref="paper", 
                    align="left",
                    showarrow=False
                )
            )
        else:
            logging.warning("No color buttons created - no suitable metadata columns found")
        
        # Set initial axis labels
        initial_model = method_models[default_method]
        if default_method == 'pca':
            x_title = f'PC1 ({initial_model.explained_variance_ratio_[0]:.1%} variance)'
            y_title = f'PC2 ({initial_model.explained_variance_ratio_[1]:.1%} variance)'
        else:
            x_title = f'{default_method.upper()} 1'
            y_title = f'{default_method.upper()} 2'
        
        fig.update_layout(
            title="Interactive K-mer Analysis",
            xaxis_title=x_title,
            yaxis_title=y_title,
            height=600,
            width=1000,
            hovermode='closest',
            template='plotly_white',
            font=dict(size=12),
            updatemenus=updatemenus,
            annotations=annotations,
            showlegend=True,
            legend=dict(
                x=1.02,
                y=1,
                xanchor='left',
                yanchor='top',
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='rgba(0,0,0,0.2)',
                borderwidth=1
            )
        )
        
        return fig

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
                            grouping_recommendations: Optional[List[Dict]] = None,
                            title: str = "MetaGrouper Analysis Report",
                            sample_id_column: str = "sample_id",
                            analyzed_variables: Optional[List[str]] = None) -> str:
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
        grouping_recommendations: Optional metadata grouping recommendations
        title: Report title
        sample_id_column: Column name for sample IDs in metadata
        analyzed_variables: Optional list of specific variables to include in plot dropdowns
        
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
        kmer_data=kmer_data,
        grouping_recommendations=grouping_recommendations,
        sample_id_column=sample_id_column,
        analyzed_variables=analyzed_variables
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