"""
Utility functions for MetaGrouper.

This module contains shared utility functions for file handling, logging setup,
result saving, and other common operations used across the package.
"""

import os
import sys
import logging
import json
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Union


def setup_logging(verbose: bool = False, log_file: str = "metagrouper.log"):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file),
        ],
    )


def find_fastq_files(input_dir: str) -> List[Tuple[Union[str, List[str]], str]]:
    """
    Find FASTQ files in input directory with paired-end support.
    
    Returns:
        List of tuples containing (filepath_or_list, sample_name)
        - For single-end: (str_path, sample_name)
        - For paired-end: ([R1_path, R2_path], sample_name)
    """
    fastq_extensions = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    files = []

    input_path = Path(input_dir)
    for ext in fastq_extensions:
        files.extend(list(input_path.glob(f"*{ext}")))
        files.extend(list(input_path.glob(f"**/*{ext}")))

    # Group files by potential sample names
    sample_files = {}
    
    for file_path in files:
        # Extract base name without extensions
        name_parts = file_path.name
        for ext in fastq_extensions:
            if name_parts.endswith(ext):
                name_parts = name_parts[:-len(ext)]
                break
        
        # Check for paired-end patterns
        paired_patterns = [
            ("_R1", "_R2"), ("_1", "_2"), (".R1", ".R2"), 
            (".1", ".2"), ("_F", "_R"), ("_forward", "_reverse")
        ]
        
        sample_name = None
        read_type = None
        
        for r1_pattern, r2_pattern in paired_patterns:
            if name_parts.endswith(r1_pattern):
                sample_name = name_parts[:-len(r1_pattern)]
                read_type = "R1"
                break
            elif name_parts.endswith(r2_pattern):
                sample_name = name_parts[:-len(r2_pattern)]
                read_type = "R2"
                break
        
        # If no paired pattern found, treat as single-end
        if sample_name is None:
            sample_name = name_parts
            read_type = "single"
        
        # Store file information
        if sample_name not in sample_files:
            sample_files[sample_name] = {}
        sample_files[sample_name][read_type] = str(file_path)

    # Create final file pairs list
    file_pairs = []
    for sample_name, file_dict in sample_files.items():
        if "R1" in file_dict and "R2" in file_dict:
            # Paired-end reads
            file_pairs.append(([file_dict["R1"], file_dict["R2"]], sample_name))
            logging.info(f"Found paired-end sample: {sample_name}")
        elif "R1" in file_dict or "R2" in file_dict:
            # Only one mate found - warn but process as single-end
            read_file = file_dict.get("R1") or file_dict.get("R2")
            file_pairs.append((read_file, sample_name))
            logging.warning(f"Only one mate found for {sample_name}, processing as single-end")
        elif "single" in file_dict:
            # Single-end reads
            file_pairs.append((file_dict["single"], sample_name))
        else:
            logging.warning(f"No valid reads found for sample: {sample_name}")

    logging.info(f"Found {len(file_pairs)} samples total")
    return sorted(file_pairs, key=lambda x: x[1])


def save_results(profiles: Dict, distance_matrix: np.ndarray, sample_names: List[str], output_dir: str):
    """Save analysis results to files."""
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    # Save profiles
    with open(output_path / "kmer_profiles.pkl", "wb") as f:
        pickle.dump(profiles, f)

    # Save distance matrix
    np.save(output_path / "distance_matrix.npy", distance_matrix)

    # Save sample names
    with open(output_path / "sample_names.json", "w") as f:
        json.dump(sample_names, f)

    # Save distance matrix as CSV
    df = pd.DataFrame(distance_matrix, index=sample_names, columns=sample_names)
    df.to_csv(output_path / "distance_matrix.csv")

    logging.info(f"Results saved to {output_dir}")


def load_results(results_dir: str) -> Tuple[Dict, np.ndarray, List[str]]:
    """Load previously saved analysis results."""
    results_path = Path(results_dir)
    
    # Load profiles
    with open(results_path / "kmer_profiles.pkl", "rb") as f:
        profiles = pickle.load(f)
    
    # Load distance matrix
    distance_matrix = np.load(results_path / "distance_matrix.npy")
    
    # Load sample names
    with open(results_path / "sample_names.json", "r") as f:
        sample_names = json.load(f)
    
    logging.info(f"Results loaded from {results_dir}")
    return profiles, distance_matrix, sample_names


def validate_input_directory(input_dir: str) -> bool:
    """Validate that input directory exists and contains FASTQ files."""
    if not os.path.exists(input_dir):
        logging.error(f"Input directory does not exist: {input_dir}")
        return False
    
    if not os.path.isdir(input_dir):
        logging.error(f"Input path is not a directory: {input_dir}")
        return False
    
    fastq_files = find_fastq_files(input_dir)
    if not fastq_files:
        logging.error(f"No FASTQ files found in {input_dir}")
        return False
    
    return True


def create_output_directory(output_dir: str, overwrite: bool = False) -> bool:
    """Create output directory with optional overwrite protection."""
    output_path = Path(output_dir)
    
    if output_path.exists() and not overwrite:
        if any(output_path.iterdir()):  # Directory exists and is not empty
            logging.warning(f"Output directory {output_dir} already exists and is not empty")
            response = input("Overwrite existing results? (y/n): ").lower().strip()
            if response != 'y':
                logging.info("Analysis cancelled by user")
                return False
    
    try:
        output_path.mkdir(parents=True, exist_ok=True)
        logging.info(f"Output directory created: {output_dir}")
        return True
    except Exception as e:
        logging.error(f"Failed to create output directory {output_dir}: {e}")
        return False


def format_file_size(size_bytes: int) -> str:
    """Format file size in human-readable format."""
    if size_bytes == 0:
        return "0 B"
    
    size_names = ["B", "KB", "MB", "GB", "TB"]
    i = 0
    while size_bytes >= 1024 and i < len(size_names) - 1:
        size_bytes /= 1024.0
        i += 1
    
    return f"{size_bytes:.1f} {size_names[i]}"


def estimate_memory_usage(num_samples: int, avg_reads_per_sample: int, k_size: int) -> dict:
    """Estimate memory usage for k-mer profiling."""
    # Rough estimates based on typical values
    avg_read_length = 150  # bp
    kmers_per_read = max(1, avg_read_length - k_size + 1)
    total_kmers = num_samples * avg_reads_per_sample * kmers_per_read
    
    # Assuming 4^k possible k-mers, but actual usage much lower
    estimated_unique_kmers = min(total_kmers * 0.1, 4**min(k_size, 15))
    
    # Memory estimates (very rough)
    kmer_memory = estimated_unique_kmers * (k_size + 8)  # k-mer string + count
    profile_memory = num_samples * estimated_unique_kmers * 8  # float64 per k-mer
    distance_matrix_memory = num_samples * num_samples * 8  # float64 matrix
    
    return {
        "estimated_kmers": int(estimated_unique_kmers),
        "kmer_storage_mb": kmer_memory / (1024 * 1024),
        "profile_storage_mb": profile_memory / (1024 * 1024),
        "distance_matrix_mb": distance_matrix_memory / (1024 * 1024),
        "total_estimated_mb": (kmer_memory + profile_memory + distance_matrix_memory) / (1024 * 1024)
    }


def summarize_fastq_files(file_list: List[Tuple[Union[str, List[str]], str]]) -> dict:
    """Summarize information about FASTQ files."""
    total_files = 0
    paired_end_samples = 0
    single_end_samples = 0
    total_size = 0
    
    for filepath, sample_name in file_list:
        if isinstance(filepath, list):
            paired_end_samples += 1
            total_files += len(filepath)
            for fp in filepath:
                if os.path.exists(fp):
                    total_size += os.path.getsize(fp)
        else:
            single_end_samples += 1
            total_files += 1
            if os.path.exists(filepath):
                total_size += os.path.getsize(filepath)
    
    return {
        "total_samples": len(file_list),
        "paired_end_samples": paired_end_samples,
        "single_end_samples": single_end_samples,
        "total_files": total_files,
        "total_size_bytes": total_size,
        "total_size_formatted": format_file_size(total_size)
    }


def check_dependencies() -> dict:
    """Check if required dependencies are available."""
    dependencies = {
        "numpy": False,
        "pandas": False,
        "matplotlib": False,
        "seaborn": False,
        "sklearn": False,
        "scipy": False
    }
    
    for dep in dependencies:
        try:
            __import__(dep)
            dependencies[dep] = True
        except ImportError:
            dependencies[dep] = False
    
    return dependencies


def get_system_info() -> dict:
    """Get system information for debugging and performance optimization."""
    import platform
    import multiprocessing
    
    return {
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "cpu_count": multiprocessing.cpu_count(),
        "architecture": platform.architecture()[0]
    }