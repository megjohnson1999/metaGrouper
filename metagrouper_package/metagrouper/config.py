"""
Configuration management for MetaGrouper.

This module provides configuration classes and functions for managing
MetaGrouper settings, defaults, and parameter validation.
"""

import os
import json
import logging
import multiprocessing
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Dict, Any, List


@dataclass
class ProfilingConfig:
    """Configuration for k-mer profiling."""
    k_size: int = 21
    max_reads: Optional[int] = None
    min_kmer_freq: int = 1
    memory_efficient: bool = True
    
    def __post_init__(self):
        if self.k_size < 1 or self.k_size > 32:
            raise ValueError(f"k_size must be between 1 and 32, got {self.k_size}")
        if self.max_reads is not None and self.max_reads < 1:
            raise ValueError(f"max_reads must be positive, got {self.max_reads}")
        if self.min_kmer_freq < 1:
            raise ValueError(f"min_kmer_freq must be positive, got {self.min_kmer_freq}")


@dataclass
class ProcessingConfig:
    """Configuration for processing parameters."""
    n_processes: Optional[int] = None
    sequential: bool = False
    show_progress: bool = True
    
    def __post_init__(self):
        if self.n_processes is None:
            self.n_processes = multiprocessing.cpu_count()
        elif self.n_processes < 1:
            raise ValueError(f"n_processes must be positive, got {self.n_processes}")


@dataclass
class AnalysisConfig:
    """Configuration for similarity analysis."""
    distance_metric: str = "braycurtis"
    memory_efficient: bool = True
    sparse_threshold: int = 50  # Use sparse computation for >50 samples
    
    valid_metrics = ["braycurtis", "jaccard", "cosine", "euclidean", "hamming", "manhattan"]
    
    def __post_init__(self):
        if self.distance_metric not in self.valid_metrics:
            raise ValueError(f"distance_metric must be one of {self.valid_metrics}, got {self.distance_metric}")
        if self.sparse_threshold < 1:
            raise ValueError(f"sparse_threshold must be positive, got {self.sparse_threshold}")


@dataclass
class VisualizationConfig:
    """Configuration for visualization."""
    figure_size: tuple = (10, 8)
    dpi: int = 300
    file_format: str = "png"
    color_map: str = "viridis"
    show_sample_labels: bool = True
    
    valid_formats = ["png", "pdf", "svg", "jpg"]
    
    def __post_init__(self):
        if self.file_format not in self.valid_formats:
            raise ValueError(f"file_format must be one of {self.valid_formats}, got {self.file_format}")
        if self.dpi < 50 or self.dpi > 1200:
            raise ValueError(f"dpi must be between 50 and 1200, got {self.dpi}")


@dataclass
class OutputConfig:
    """Configuration for output settings."""
    output_dir: str = "metagrouper_output"
    save_profiles: bool = True
    save_distance_matrix: bool = True
    save_plots: bool = True
    overwrite: bool = False
    log_file: str = "metagrouper.log"
    verbose: bool = False


class MetaGrouperConfig:
    """Main configuration class for MetaGrouper."""
    
    def __init__(self, config_file: Optional[str] = None):
        # Initialize with defaults
        self.profiling = ProfilingConfig()
        self.processing = ProcessingConfig()
        self.analysis = AnalysisConfig()
        self.visualization = VisualizationConfig()
        self.output = OutputConfig()
        
        # Load from file if provided
        if config_file:
            self.load_from_file(config_file)
    
    def load_from_file(self, config_file: str):
        """Load configuration from JSON file."""
        config_path = Path(config_file)
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_file}")
        
        try:
            with open(config_path, 'r') as f:
                config_data = json.load(f)
            
            # Update configurations
            if 'profiling' in config_data:
                self.profiling = ProfilingConfig(**config_data['profiling'])
            if 'processing' in config_data:
                self.processing = ProcessingConfig(**config_data['processing'])
            if 'analysis' in config_data:
                self.analysis = AnalysisConfig(**config_data['analysis'])
            if 'visualization' in config_data:
                self.visualization = VisualizationConfig(**config_data['visualization'])
            if 'output' in config_data:
                self.output = OutputConfig(**config_data['output'])
                
            logging.info(f"Configuration loaded from {config_file}")
            
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in configuration file: {e}")
        except TypeError as e:
            raise ValueError(f"Invalid configuration parameters: {e}")
    
    def save_to_file(self, config_file: str):
        """Save current configuration to JSON file."""
        config_data = {
            'profiling': asdict(self.profiling),
            'processing': asdict(self.processing),
            'analysis': asdict(self.analysis),
            'visualization': asdict(self.visualization),
            'output': asdict(self.output)
        }
        
        config_path = Path(config_file)
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(config_path, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        logging.info(f"Configuration saved to {config_file}")
    
    def update_from_args(self, args):
        """Update configuration from command-line arguments."""
        # Profiling parameters
        if hasattr(args, 'kmer_size') and args.kmer_size:
            self.profiling.k_size = args.kmer_size
        if hasattr(args, 'max_reads') and args.max_reads:
            self.profiling.max_reads = args.max_reads
        if hasattr(args, 'min_kmer_freq') and args.min_kmer_freq:
            self.profiling.min_kmer_freq = args.min_kmer_freq
        if hasattr(args, 'memory_efficient'):
            self.profiling.memory_efficient = args.memory_efficient and not getattr(args, 'no_memory_efficient', False)
        
        # Processing parameters
        if hasattr(args, 'processes') and args.processes:
            self.processing.n_processes = args.processes
        if hasattr(args, 'sequential') and args.sequential:
            self.processing.sequential = args.sequential
        
        # Analysis parameters
        if hasattr(args, 'distance_metric') and args.distance_metric:
            self.analysis.distance_metric = args.distance_metric
        
        # Output parameters
        if hasattr(args, 'output') and args.output:
            self.output.output_dir = args.output
        if hasattr(args, 'verbose') and args.verbose:
            self.output.verbose = args.verbose
    
    def validate(self):
        """Validate all configuration parameters."""
        # Validate each configuration section
        try:
            # This will trigger __post_init__ validation
            ProfilingConfig(**asdict(self.profiling))
            ProcessingConfig(**asdict(self.processing))
            AnalysisConfig(**asdict(self.analysis))
            VisualizationConfig(**asdict(self.visualization))
        except ValueError as e:
            raise ValueError(f"Configuration validation failed: {e}")
        
        # Additional cross-validation
        if self.processing.sequential and self.processing.n_processes > 1:
            logging.warning("Sequential mode enabled but n_processes > 1. Will use sequential processing.")
    
    def get_summary(self) -> str:
        """Get a human-readable summary of the configuration."""
        summary = []
        summary.append("MetaGrouper Configuration Summary:")
        summary.append("=" * 40)
        
        summary.append(f"Profiling:")
        summary.append(f"  K-mer size: {self.profiling.k_size}")
        summary.append(f"  Max reads per sample: {self.profiling.max_reads or 'unlimited'}")
        summary.append(f"  Min k-mer frequency: {self.profiling.min_kmer_freq}")
        summary.append(f"  Memory efficient: {self.profiling.memory_efficient}")
        
        summary.append(f"Processing:")
        summary.append(f"  Processes: {self.processing.n_processes}")
        summary.append(f"  Sequential mode: {self.processing.sequential}")
        summary.append(f"  Show progress: {self.processing.show_progress}")
        
        summary.append(f"Analysis:")
        summary.append(f"  Distance metric: {self.analysis.distance_metric}")
        summary.append(f"  Memory efficient: {self.analysis.memory_efficient}")
        summary.append(f"  Sparse threshold: {self.analysis.sparse_threshold} samples")
        
        summary.append(f"Output:")
        summary.append(f"  Output directory: {self.output.output_dir}")
        summary.append(f"  Verbose logging: {self.output.verbose}")
        
        return "\n".join(summary)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'profiling': asdict(self.profiling),
            'processing': asdict(self.processing),
            'analysis': asdict(self.analysis),
            'visualization': asdict(self.visualization),
            'output': asdict(self.output)
        }


def create_default_config(output_file: str = "metagrouper_config.json"):
    """Create a default configuration file."""
    config = MetaGrouperConfig()
    config.save_to_file(output_file)
    return config


def load_config_from_env() -> MetaGrouperConfig:
    """Load configuration from environment variables."""
    config = MetaGrouperConfig()
    
    # Profiling configuration from environment
    if 'METAGROUPER_K_SIZE' in os.environ:
        config.profiling.k_size = int(os.environ['METAGROUPER_K_SIZE'])
    if 'METAGROUPER_MAX_READS' in os.environ:
        config.profiling.max_reads = int(os.environ['METAGROUPER_MAX_READS'])
    if 'METAGROUPER_MIN_KMER_FREQ' in os.environ:
        config.profiling.min_kmer_freq = int(os.environ['METAGROUPER_MIN_KMER_FREQ'])
    if 'METAGROUPER_MEMORY_EFFICIENT' in os.environ:
        config.profiling.memory_efficient = os.environ['METAGROUPER_MEMORY_EFFICIENT'].lower() == 'true'
    
    # Processing configuration from environment
    if 'METAGROUPER_PROCESSES' in os.environ:
        config.processing.n_processes = int(os.environ['METAGROUPER_PROCESSES'])
    if 'METAGROUPER_SEQUENTIAL' in os.environ:
        config.processing.sequential = os.environ['METAGROUPER_SEQUENTIAL'].lower() == 'true'
    
    # Analysis configuration from environment
    if 'METAGROUPER_DISTANCE_METRIC' in os.environ:
        config.analysis.distance_metric = os.environ['METAGROUPER_DISTANCE_METRIC']
    
    # Output configuration from environment
    if 'METAGROUPER_OUTPUT_DIR' in os.environ:
        config.output.output_dir = os.environ['METAGROUPER_OUTPUT_DIR']
    if 'METAGROUPER_VERBOSE' in os.environ:
        config.output.verbose = os.environ['METAGROUPER_VERBOSE'].lower() == 'true'
    
    return config


def get_recommended_config(num_samples: int, total_file_size_gb: float) -> MetaGrouperConfig:
    """Get recommended configuration based on dataset characteristics."""
    config = MetaGrouperConfig()
    
    # Adjust k-mer size based on dataset complexity
    if num_samples > 100:
        config.profiling.k_size = 19  # Smaller k for large datasets
        config.profiling.min_kmer_freq = 3  # More aggressive filtering
    elif num_samples > 50:
        config.profiling.k_size = 21  # Standard k
        config.profiling.min_kmer_freq = 2
    else:
        config.profiling.k_size = 23  # Larger k for small datasets
        config.profiling.min_kmer_freq = 1
    
    # Adjust processing based on dataset size
    if total_file_size_gb > 10:
        config.profiling.memory_efficient = True
        config.analysis.memory_efficient = True
        config.processing.n_processes = min(multiprocessing.cpu_count(), 8)
    else:
        config.processing.n_processes = min(multiprocessing.cpu_count(), 4)
    
    # Sequential for very small datasets
    if num_samples <= 2:
        config.processing.sequential = True
    
    return config