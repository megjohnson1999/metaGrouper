"""
MetaGrouper: K-mer-based Analysis for Optimal Metagenomic Assembly Grouping

A modular bioinformatics tool for metagenomic sample analysis and assembly grouping.
"""

from .profiler import KmerProfiler, process_sample_worker
from .analyzer import SimilarityAnalyzer
from .visualizer import Visualizer
from .utils import find_fastq_files, setup_logging, save_results
from .config import MetaGrouperConfig

__version__ = "2.0.0"
__author__ = "MetaGrouper Development Team"

__all__ = [
    "KmerProfiler",
    "process_sample_worker", 
    "SimilarityAnalyzer",
    "Visualizer",
    "find_fastq_files",
    "setup_logging",
    "save_results",
    "MetaGrouperConfig"
]