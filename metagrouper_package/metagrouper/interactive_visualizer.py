"""
Interactive visualization functionality for metagenomic analysis results.

This module contains the InteractiveVisualizer class for generating interactive HTML plots
and visualizations of sample relationships, distance matrices, and dimensionality reduction results.
"""

import logging
import numpy as np
import pandas as pd
from typing import List, Dict, Optional, Any, Tuple
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline as pyo

# Optional advanced dimensionality reduction
try:
    import umap
    UMAP_AVAILABLE = True
except ImportError:
    UMAP_AVAILABLE = False
    logging.warning("UMAP not available - install umap-learn for UMAP projections")


class InteractiveVisualizer:
    """Generate interactive visualizations for sample relationships."""

    def __init__(self, sample_names: List[str], metadata: Optional[pd.DataFrame] = None):
        """
        Initialize the InteractiveVisualizer.
        
        Args:
            sample_names: List of sample names
            metadata: Optional metadata DataFrame with sample information
        """
        self.sample_names = sample_names
        self.metadata = metadata
        self.figures = {}
        self.projections_cache = {}  # Cache for expensive computations
        
        # Validate metadata if provided
        if self.metadata is not None:
            self._validate_metadata()
    
    def _validate_metadata(self):
        """Validate that metadata contains samples and prepare for visualization."""
        if self.metadata is None:
            return
        
        # Check if sample_id column exists or if index contains sample names
        if 'sample_id' in self.metadata.columns:
            sample_col = 'sample_id'
        elif self.metadata.index.name == 'sample_id' or any(name in self.metadata.index for name in self.sample_names):
            sample_col = self.metadata.index.name or 'index'
        else:
            # Try to find samples in any column
            for col in self.metadata.columns:
                if any(name in self.metadata[col].values for name in self.sample_names):
                    sample_col = col
                    break
            else:
                logging.warning("Could not find sample names in metadata")
                return
        
        # Filter metadata to only include samples we have
        if sample_col == 'index' or sample_col == self.metadata.index.name:
            self.metadata = self.metadata[self.metadata.index.isin(self.sample_names)]
        else:
            self.metadata = self.metadata[self.metadata[sample_col].isin(self.sample_names)]
        
        logging.info(f"Loaded metadata for {len(self.metadata)} samples with {len(self.metadata.columns)} variables")

    def compute_projections(self, data: np.ndarray, methods: List[str] = None) -> Dict[str, Tuple[np.ndarray, Any]]:
        """
        Compute multiple dimensionality reduction projections.
        
        Args:
            data: Input data matrix (samples x features)
            methods: List of methods to compute ['pca', 'tsne', 'umap']
            
        Returns:
            Dictionary mapping method names to (projection, fitted_model) tuples
        """
        if methods is None:
            methods = ['pca', 'tsne']
            if UMAP_AVAILABLE:
                methods.append('umap')
        
        projections = {}
        
        for method in methods:
            cache_key = f"{method}_{data.shape}_{hash(data.tobytes())}"
            if cache_key in self.projections_cache:
                projections[method] = self.projections_cache[cache_key]
                continue
                
            logging.info(f"Computing {method.upper()} projection...")
            
            try:
                if method == 'pca':
                    model = PCA(n_components=2)
                    projection = model.fit_transform(data)
                    
                elif method == 'tsne':
                    # Use appropriate perplexity based on sample size
                    perplexity = min(30, max(5, len(self.sample_names) // 4))
                    model = TSNE(n_components=2, perplexity=perplexity, random_state=42, 
                               max_iter=1000, learning_rate='auto', init='pca')
                    projection = model.fit_transform(data)
                    
                elif method == 'umap' and UMAP_AVAILABLE:
                    # Use appropriate n_neighbors based on sample size
                    n_neighbors = min(15, max(2, len(self.sample_names) // 3))
                    model = umap.UMAP(n_components=2, n_neighbors=n_neighbors, 
                                    random_state=42, min_dist=0.1)
                    projection = model.fit_transform(data)
                    
                else:
                    logging.warning(f"Method {method} not available, skipping")
                    continue
                
                projections[method] = (projection, model)
                self.projections_cache[cache_key] = (projection, model)
                logging.info(f"✅ {method.upper()} projection computed")
                
            except Exception as e:
                logging.warning(f"Failed to compute {method} projection: {e}")
                continue
        
        return projections

    def create_enhanced_dimensionality_plot(self, data: np.ndarray, output_path: str, 
                                           title: str = "Enhanced Dimensionality Analysis",
                                           methods: List[str] = None) -> str:
        """
        Create an enhanced interactive plot with multiple dimensionality reduction methods.
        
        Args:
            data: Input data matrix (samples x features)
            output_path: Path to save HTML file
            title: Plot title
            methods: List of methods to compute ['pca', 'tsne', 'umap']
            
        Returns:
            Path to generated HTML file
        """
        logging.info("Creating enhanced multi-method dimensionality plot")
        
        # Compute all projections
        projections = self.compute_projections(data, methods)
        
        if not projections:
            logging.error("No valid projections computed")
            return None
        
        # Create base DataFrame
        plot_data = {}
        method_models = {}
        
        for method, (projection, model) in projections.items():
            plot_data[f'{method}_x'] = projection[:, 0]
            plot_data[f'{method}_y'] = projection[:, 1]
            method_models[method] = model
        
        plot_df = pd.DataFrame({
            'sample_id': self.sample_names,
            **plot_data
        })
        
        # Merge with metadata if available
        if self.metadata is not None:
            metadata_for_merge = self.metadata.reset_index()
            if 'sample_id' not in metadata_for_merge.columns:
                metadata_for_merge['sample_id'] = metadata_for_merge.index
            plot_df = plot_df.merge(metadata_for_merge, on='sample_id', how='left')
        
        # Create the enhanced plot
        fig = self._create_enhanced_plot(plot_df, projections, method_models, title)
        
        # Save as HTML
        html_path = Path(output_path)
        pyo.plot(fig, filename=str(html_path), auto_open=False)
        
        logging.info(f"Enhanced dimensionality plot saved to {html_path}")
        return str(html_path)

    def _create_enhanced_plot(self, plot_df: pd.DataFrame, projections: Dict, 
                             method_models: Dict, title: str) -> go.Figure:
        """Create enhanced plot with method switching and advanced interactions."""
        
        # Get metadata columns for coloring options and convert categorical to numerical
        metadata_cols = []
        color_mappings = {}
        if self.metadata is not None:
            for col in plot_df.columns:
                if col not in ['sample_id'] and not col.endswith(('_x', '_y')):
                    n_unique = plot_df[col].nunique()
                    if n_unique > 1 and n_unique <= 20:
                        metadata_cols.append(col)
                        # Convert categorical to numerical for Plotly
                        if plot_df[col].dtype == 'object':
                            unique_vals = plot_df[col].unique()
                            mapping = {val: i for i, val in enumerate(unique_vals)}
                            plot_df[f'{col}_numeric'] = plot_df[col].map(mapping)
                            color_mappings[col] = mapping
                    elif n_unique > 20 and plot_df[col].dtype in ['int64', 'float64']:
                        metadata_cols.append(col)
        
        # Start with the first available method
        default_method = list(projections.keys())[0]
        default_color = metadata_cols[0] if metadata_cols else None
        
        # Create initial scatter plot
        fig = go.Figure()
        
        # Add traces for each method (initially hidden except default)
        for method_name, (projection, model) in projections.items():
            visible = method_name == default_method
            
            # Determine axis labels
            if method_name == 'pca':
                x_label = f'PC1 ({model.explained_variance_ratio_[0]:.1%} variance)'
                y_label = f'PC2 ({model.explained_variance_ratio_[1]:.1%} variance)'
            else:
                x_label = f'{method_name.upper()} 1'
                y_label = f'{method_name.upper()} 2'
            
            # Prepare hover data
            hover_data = []
            if metadata_cols:
                for col in metadata_cols:
                    hover_data.append(plot_df[col].tolist())
            
            # Base trace
            scatter_data = dict(
                x=plot_df[f'{method_name}_x'],
                y=plot_df[f'{method_name}_y'],
                text=plot_df['sample_id'],
                customdata=list(zip(*hover_data)) if hover_data else None,
                mode='markers',
                marker=dict(
                    size=10,
                    opacity=0.8,
                    line=dict(width=1, color='DarkSlateGrey'),
                    color=self._get_color_array(plot_df, default_color, color_mappings),
                    colorscale='viridis' if default_color else None,
                    colorbar=dict(
                        title=dict(text=default_color.replace('_', ' ').title()) if default_color else None
                    ) if default_color else None
                ),
                name=method_name.upper(),
                visible=visible,
                hovertemplate=self._create_hover_template(plot_df, method_name, metadata_cols)
            )
            
            fig.add_trace(go.Scatter(**scatter_data))
        
        # Create buttons for method switching
        method_buttons = []
        for i, (method_name, (projection, model)) in enumerate(projections.items()):
            # Determine axis labels
            if method_name == 'pca':
                x_label = f'PC1 ({model.explained_variance_ratio_[0]:.1%} variance)'
                y_label = f'PC2 ({model.explained_variance_ratio_[1]:.1%} variance)'
            else:
                x_label = f'{method_name.upper()} 1'
                y_label = f'{method_name.upper()} 2'
            
            visibility = [False] * len(projections)
            visibility[i] = True
            
            method_buttons.append(
                dict(
                    label=method_name.upper(),
                    method="update",
                    args=[
                        {"visible": visibility},
                        {"xaxis.title": x_label, "yaxis.title": y_label}
                    ]
                )
            )
        
        # Create color buttons if metadata available
        color_buttons = []
        if metadata_cols:
            for col in metadata_cols:
                # Create color arrays for each trace/method
                color_arrays = []
                color_data = self._get_color_array(plot_df, col, color_mappings)
                for _ in projections:
                    color_arrays.append(color_data)
                
                # Create colorbar update for this variable
                colorbar_update = {
                    "marker.color": color_arrays,
                    "marker.colorbar.title.text": col.replace('_', ' ').title()
                }
                
                color_buttons.append(
                    dict(
                        label=col.replace('_', ' ').title(),
                        method="restyle",
                        args=[colorbar_update]
                    )
                )
            
            # Add no coloring option
            color_buttons.append(
                dict(
                    label="No Coloring",
                    method="restyle",
                    args=[{
                        "marker.color": ['blue'] * len(projections),
                        "marker.colorbar.title.text": ""
                    }]
                )
            )
        
        # Update layout with dropdowns
        updatemenus = []
        
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
        
        # Color selector (if metadata available)
        if color_buttons:
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
        
        annotations = [
            dict(
                text="Method:",
                x=0.05, y=1.18,
                xref="paper", yref="paper",
                align="left",
                showarrow=False
            )
        ]
        
        if color_buttons:
            annotations.append(
                dict(
                    text="Color by:",
                    x=0.35, y=1.18,
                    xref="paper", yref="paper", 
                    align="left",
                    showarrow=False
                )
            )
        
        # Set initial axis labels
        initial_model = method_models[default_method]
        if default_method == 'pca':
            x_title = f'PC1 ({initial_model.explained_variance_ratio_[0]:.1%} variance)'
            y_title = f'PC2 ({initial_model.explained_variance_ratio_[1]:.1%} variance)'
        else:
            x_title = f'{default_method.upper()} 1'
            y_title = f'{default_method.upper()} 2'
        
        fig.update_layout(
            title=title,
            xaxis_title=x_title,
            yaxis_title=y_title,
            height=700,
            width=1000,
            hovermode='closest',
            template='plotly_white',
            font=dict(size=12),
            updatemenus=updatemenus,
            annotations=annotations
        )
        
        return fig

    def _get_color_array(self, plot_df: pd.DataFrame, col: str, color_mappings: Dict) -> List:
        """Get appropriate color array for a column, handling categorical vs numerical."""
        if col is None:
            return 'blue'
        
        if col in color_mappings:
            # Use numeric mapping for categorical data
            return plot_df[f'{col}_numeric'].tolist()
        else:
            # Use original values for numerical data
            return plot_df[col].tolist()

    def _create_categorical_legend_traces(self, plot_df: pd.DataFrame, col: str, color_mappings: Dict) -> List:
        """Create invisible traces for categorical legend."""
        if col not in color_mappings:
            return []
        
        legend_traces = []
        unique_vals = list(color_mappings[col].keys())
        
        # Use discrete colors for categories
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
        
        for i, val in enumerate(unique_vals):
            legend_traces.append(
                go.Scatter(
                    x=[None], y=[None],  # Invisible points
                    mode='markers',
                    marker=dict(size=10, color=colors[i % len(colors)]),
                    name=str(val),
                    showlegend=True,
                    visible='legendonly'  # Only show in legend
                )
            )
        
        return legend_traces

    def _create_hover_template(self, plot_df: pd.DataFrame, method_name: str, metadata_cols: List[str]) -> str:
        """Create hover template for enhanced plots."""
        template = f"<b>%{{text}}</b><br>"
        template += f"{method_name.upper()} 1: %{{x:.3f}}<br>"
        template += f"{method_name.upper()} 2: %{{y:.3f}}<br>"
        
        # Add metadata info to hover
        for i, col in enumerate(metadata_cols):
            template += f"{col.replace('_', ' ').title()}: %{{customdata[{i}]}}<br>"
        
        template += "<extra></extra>"
        return template

    def create_interactive_pca(self, pca_result: np.ndarray, pca: PCA, 
                              output_path: str, title: str = "Interactive PCA Analysis") -> str:
        """
        Create an interactive PCA plot with metadata coloring options.
        
        Args:
            pca_result: PCA transformed data
            pca: Fitted PCA object
            output_path: Path to save HTML file
            title: Plot title
            
        Returns:
            Path to generated HTML file
        """
        logging.info("Creating interactive PCA plot")
        
        # Create base DataFrame for plotting
        plot_df = pd.DataFrame({
            'sample_id': self.sample_names,
            'PC1': pca_result[:, 0],
            'PC2': pca_result[:, 1],
        })
        
        # Add PC3 if available
        if pca_result.shape[1] > 2:
            plot_df['PC3'] = pca_result[:, 2]
        
        # Merge with metadata if available
        if self.metadata is not None:
            # Reset index to ensure proper merging
            metadata_for_merge = self.metadata.reset_index()
            if 'sample_id' not in metadata_for_merge.columns:
                metadata_for_merge['sample_id'] = metadata_for_merge.index
            
            plot_df = plot_df.merge(metadata_for_merge, on='sample_id', how='left')
        
        # Create the interactive plot
        fig = self._create_pca_plot(plot_df, pca, title)
        
        # Save as HTML
        html_path = Path(output_path)
        pyo.plot(fig, filename=str(html_path), auto_open=False)
        
        logging.info(f"Interactive PCA plot saved to {html_path}")
        return str(html_path)
    
    def _create_pca_plot(self, plot_df: pd.DataFrame, pca: PCA, title: str) -> go.Figure:
        """Create the PCA plot figure with interactive elements."""
        
        # Get metadata columns for coloring options
        metadata_cols = []
        if self.metadata is not None:
            # Get categorical and numerical columns
            for col in plot_df.columns:
                if col not in ['sample_id', 'PC1', 'PC2', 'PC3']:
                    # Check if column has reasonable number of unique values for coloring
                    n_unique = plot_df[col].nunique()
                    if n_unique > 1 and n_unique <= 20:  # Good for categorical coloring
                        metadata_cols.append(col)
                    elif n_unique > 20 and plot_df[col].dtype in ['int64', 'float64']:  # Continuous
                        metadata_cols.append(col)
        
        # Create base scatter plot
        if len(metadata_cols) > 0:
            # Create plot with default coloring by first metadata variable
            default_color = metadata_cols[0]
            fig = px.scatter(
                plot_df, 
                x='PC1', 
                y='PC2',
                color=default_color,
                hover_data=['sample_id'] + metadata_cols,
                title=title,
                labels={
                    'PC1': f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)',
                    'PC2': f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)'
                }
            )
            
            # Add dropdown for color selection
            buttons = []
            for col in metadata_cols:
                buttons.append(
                    dict(
                        label=col.replace('_', ' ').title(),
                        method="restyle",
                        args=[{"marker.color": plot_df[col]}]
                    )
                )
            
            # Add option for no coloring
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
        else:
            # Simple plot without metadata coloring
            fig = px.scatter(
                plot_df,
                x='PC1',
                y='PC2', 
                hover_data=['sample_id'],
                title=title,
                labels={
                    'PC1': f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)',
                    'PC2': f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)'
                }
            )
        
        # Enhance the plot appearance
        fig.update_traces(
            marker=dict(size=10, opacity=0.8, line=dict(width=1, color='DarkSlateGrey')),
            selector=dict(mode='markers')
        )
        
        fig.update_layout(
            height=600,
            width=900,
            hovermode='closest',
            template='plotly_white',
            font=dict(size=12)
        )
        
        return fig

    def create_interactive_heatmap(self, distance_matrix: np.ndarray, 
                                 output_path: str, title: str = "Interactive Distance Heatmap") -> str:
        """
        Create an interactive heatmap of the distance matrix.
        
        Args:
            distance_matrix: Distance matrix between samples
            output_path: Path to save HTML file
            title: Plot title
            
        Returns:
            Path to generated HTML file
        """
        logging.info("Creating interactive distance heatmap")
        
        # Create the heatmap
        fig = go.Figure(data=go.Heatmap(
            z=distance_matrix,
            x=self.sample_names,
            y=self.sample_names,
            colorscale='Viridis',
            hoverongaps=False,
            hovertemplate='Sample 1: %{y}<br>Sample 2: %{x}<br>Distance: %{z:.3f}<extra></extra>'
        ))
        
        fig.update_layout(
            title=title,
            xaxis_title="Samples",
            yaxis_title="Samples", 
            height=600,
            width=700,
            template='plotly_white'
        )
        
        # Save as HTML
        html_path = Path(output_path)
        pyo.plot(fig, filename=str(html_path), auto_open=False)
        
        logging.info(f"Interactive heatmap saved to {html_path}")
        return str(html_path)

    def create_dashboard(self, pca_result: np.ndarray, pca: PCA, 
                        distance_matrix: np.ndarray, output_path: str,
                        title: str = "MetaGrouper Interactive Dashboard") -> str:
        """
        Create a unified dashboard with multiple interactive visualizations.
        
        Args:
            pca_result: PCA transformed data
            pca: Fitted PCA object
            distance_matrix: Distance matrix between samples
            output_path: Path to save HTML file
            title: Dashboard title
            
        Returns:
            Path to generated HTML file
        """
        logging.info("Creating interactive dashboard")
        
        # Create base DataFrame for plotting
        plot_df = pd.DataFrame({
            'sample_id': self.sample_names,
            'PC1': pca_result[:, 0],
            'PC2': pca_result[:, 1],
        })
        
        # Merge with metadata if available
        if self.metadata is not None:
            metadata_for_merge = self.metadata.reset_index()
            if 'sample_id' not in metadata_for_merge.columns:
                metadata_for_merge['sample_id'] = metadata_for_merge.index
            plot_df = plot_df.merge(metadata_for_merge, on='sample_id', how='left')
        
        # Create subplots
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=('PCA Analysis', 'Distance Heatmap'),
            specs=[[{"type": "scatter"}, {"type": "heatmap"}]]
        )
        
        # Add PCA plot
        pca_fig = self._create_pca_plot(plot_df, pca, "")
        for trace in pca_fig.data:
            fig.add_trace(trace, row=1, col=1)
        
        # Add heatmap
        fig.add_trace(
            go.Heatmap(
                z=distance_matrix,
                x=self.sample_names,
                y=self.sample_names,
                colorscale='Viridis',
                hoverongaps=False,
                hovertemplate='Sample 1: %{y}<br>Sample 2: %{x}<br>Distance: %{z:.3f}<extra></extra>',
                showscale=True
            ),
            row=1, col=2
        )
        
        # Update layout
        fig.update_layout(
            title=title,
            height=600,
            width=1400,
            template='plotly_white',
            showlegend=True
        )
        
        # Update axis labels
        fig.update_xaxes(title_text=f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)', row=1, col=1)
        fig.update_yaxes(title_text=f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)', row=1, col=1)
        fig.update_xaxes(title_text="Samples", row=1, col=2)
        fig.update_yaxes(title_text="Samples", row=1, col=2)
        
        # Save as HTML
        html_path = Path(output_path)
        pyo.plot(fig, filename=str(html_path), auto_open=False)
        
        logging.info(f"Interactive dashboard saved to {html_path}")
        return str(html_path)