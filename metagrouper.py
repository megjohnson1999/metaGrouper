#!/usr/bin/env python3
"""
MetaGrouper: K-mer-based Analysis for Optimal Metagenomic Assembly Grouping

This script computes k-mer profiles for metagenomic samples, analyzes their
similarities, and uses metadata to recommend optimal grouping strategies for assembly.

Phases:
1. K-mer profiling and sample similarity analysis
2. Metadata variable analysis with PERMANOVA
"""

import os
import sys
import argparse
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import gzip
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.metrics.pairwise import pairwise_distances
import pickle
import json

# Import Phase 2 functionality
try:
    from metadata_analyzer import (
        MetadataAnalyzer,
        MetadataVisualizer,
        generate_summary_report,
    )

    PHASE2_AVAILABLE = True
except ImportError:
    PHASE2_AVAILABLE = False
    logging.warning(
        "Phase 2 metadata analysis not available. Install required dependencies."
    )

# Import Phase 3 functionality
try:
    from assembly_recommender import (
        AssemblyRecommender,
        save_recommendations,
        visualize_assembly_strategy,
    )

    PHASE3_AVAILABLE = True
except ImportError:
    PHASE3_AVAILABLE = False
    logging.warning(
        "Phase 3 assembly recommendations not available. Install required dependencies."
    )


class KmerProfiler:
    """Handle k-mer profiling of metagenomic samples."""

    def __init__(self, k: int = 21, max_reads: Optional[int] = None):
        self.k = k
        self.max_reads = max_reads
        self.profiles = {}
        self.sample_names = []

    def _open_file(self, filepath: str):
        """Open file, handling gzip compression."""
        if filepath.endswith(".gz"):
            return gzip.open(filepath, "rt")
        return open(filepath, "r")

    def _extract_kmers(self, sequence: str) -> Dict[str, int]:
        """Extract k-mers from a sequence."""
        kmers = defaultdict(int)
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i : i + self.k]
            if "N" not in kmer:  # Skip k-mers with ambiguous bases
                # Use canonical k-mer (lexicographically smaller of forward/reverse)
                rev_comp = self._reverse_complement(kmer)
                canonical = min(kmer, rev_comp)
                kmers[canonical] += 1
        return kmers

    def _reverse_complement(self, seq: str) -> str:
        """Return reverse complement of a DNA sequence."""
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return "".join(complement.get(base, base) for base in reversed(seq))

    def _validate_fastq_record(self, header: str, sequence: str, plus: str, quality: str, line_num: int, filepath: str) -> bool:
        """Validate a single FASTQ record."""
        if not header.startswith('@'):
            logging.error(f"Invalid FASTQ header at line {line_num} in {filepath}: {header[:50]}...")
            return False
        
        if not plus.startswith('+'):
            logging.error(f"Invalid FASTQ plus line at line {line_num + 2} in {filepath}: {plus[:50]}...")
            return False
        
        if len(sequence) != len(quality):
            logging.error(f"Sequence and quality length mismatch at line {line_num} in {filepath}")
            return False
        
        # Check for valid DNA bases
        valid_bases = set('ATCGN')
        if not all(base in valid_bases for base in sequence.upper()):
            logging.warning(f"Non-standard bases found in sequence at line {line_num} in {filepath}")
        
        return True

    def _parse_fastq(self, filepath: Union[str, List[str]]) -> List[str]:
        """
        Parse FASTQ file(s) and return sequences with validation.
        
        Args:
            filepath: Single file path (str) or list of paired-end file paths
            
        Returns:
            List of sequences from all input files
        """
        sequences = []
        read_count = 0
        
        # Handle both single files and paired-end file lists
        file_paths = [filepath] if isinstance(filepath, str) else filepath
        
        for file_path in file_paths:
            line_num = 0
            try:
                with self._open_file(file_path) as f:
                    while True:
                        if self.max_reads and read_count >= self.max_reads:
                            break

                        line_num += 1
                        header = f.readline().strip()
                        if not header:
                            break

                        line_num += 1
                        sequence = f.readline().strip()
                        line_num += 1
                        plus = f.readline().strip()
                        line_num += 1
                        quality = f.readline().strip()

                        if header and sequence and plus and quality:
                            if self._validate_fastq_record(header, sequence, plus, quality, line_num - 3, file_path):
                                sequences.append(sequence.upper())
                                read_count += 1
                        elif header:  # Incomplete record at end of file
                            logging.warning(f"Incomplete FASTQ record at end of {file_path}")
                            break
                            
            except Exception as e:
                raise ValueError(f"Error parsing FASTQ file {file_path} at line {line_num}: {str(e)}")

        return sequences

    def profile_sample(self, filepath: Union[str, List[str]], sample_name: str) -> Dict[str, int]:
        """Generate k-mer profile for a single sample (single-end or paired-end)."""
        if isinstance(filepath, list):
            logging.info(f"Processing paired-end sample: {sample_name} ({len(filepath)} files)")
        else:
            logging.info(f"Processing single-end sample: {sample_name}")

        sequences = self._parse_fastq(filepath)
        logging.info(f"Loaded {len(sequences)} sequences from {sample_name}")

        total_kmers = defaultdict(int)
        for seq in sequences:
            seq_kmers = self._extract_kmers(seq)
            for kmer, count in seq_kmers.items():
                total_kmers[kmer] += count

        # Normalize by total k-mer count
        total_count = sum(total_kmers.values())
        normalized_profile = {
            kmer: count / total_count for kmer, count in total_kmers.items()
        }

        self.profiles[sample_name] = normalized_profile
        if sample_name not in self.sample_names:
            self.sample_names.append(sample_name)

        logging.info(f"Generated profile with {len(normalized_profile)} unique k-mers")
        return normalized_profile


class SimilarityAnalyzer:
    """Analyze similarities between k-mer profiles."""

    def __init__(self, profiles: Dict[str, Dict[str, float]]):
        self.profiles = profiles
        self.sample_names = list(profiles.keys())
        self.distance_matrix = None
        self.similarity_matrix = None

    def compute_distance_matrix(self, metric: str = "braycurtis") -> np.ndarray:
        """Compute pairwise distance matrix between samples."""
        logging.info(f"Computing distance matrix using {metric} metric")

        # Get all unique k-mers across all samples
        all_kmers = set()
        for profile in self.profiles.values():
            all_kmers.update(profile.keys())
        all_kmers = sorted(list(all_kmers))

        # Create feature matrix
        feature_matrix = np.zeros((len(self.sample_names), len(all_kmers)))
        for i, sample in enumerate(self.sample_names):
            for j, kmer in enumerate(all_kmers):
                feature_matrix[i, j] = self.profiles[sample].get(kmer, 0)

        # Compute pairwise distances
        self.distance_matrix = pairwise_distances(feature_matrix, metric=metric)
        self.similarity_matrix = 1 - self.distance_matrix

        return self.distance_matrix

    def perform_pca(self, n_components: int = 2) -> Tuple[np.ndarray, PCA]:
        """Perform PCA on k-mer profiles."""
        logging.info("Performing PCA analysis")

        # Get all unique k-mers
        all_kmers = set()
        for profile in self.profiles.values():
            all_kmers.update(profile.keys())
        all_kmers = sorted(list(all_kmers))

        # Create feature matrix
        feature_matrix = np.zeros((len(self.sample_names), len(all_kmers)))
        for i, sample in enumerate(self.sample_names):
            for j, kmer in enumerate(all_kmers):
                feature_matrix[i, j] = self.profiles[sample].get(kmer, 0)

        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(feature_matrix)

        logging.info(f"PCA explained variance ratio: {pca.explained_variance_ratio_}")
        return pca_result, pca

    def perform_mds(self, n_components: int = 2) -> np.ndarray:
        """Perform MDS on distance matrix."""
        logging.info("Performing MDS analysis")

        if self.distance_matrix is None:
            self.compute_distance_matrix()

        mds = MDS(
            n_components=n_components, dissimilarity="precomputed", random_state=42
        )
        mds_result = mds.fit_transform(self.distance_matrix)

        return mds_result


class Visualizer:
    """Generate visualizations for sample relationships."""

    def __init__(self, sample_names: List[str]):
        self.sample_names = sample_names

    def plot_distance_heatmap(self, distance_matrix: np.ndarray, output_path: str):
        """Plot distance matrix as heatmap."""
        plt.figure(figsize=(10, 8))

        # Create DataFrame for better labeling
        df = pd.DataFrame(
            distance_matrix, index=self.sample_names, columns=self.sample_names
        )

        sns.heatmap(df, annot=True, cmap="viridis", fmt=".3f")
        plt.title("Sample Distance Matrix")
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"Distance heatmap saved to {output_path}")

    def plot_pca(self, pca_result: np.ndarray, pca: PCA, output_path: str):
        """Plot PCA results."""
        plt.figure(figsize=(10, 8))

        plt.scatter(pca_result[:, 0], pca_result[:, 1], s=100, alpha=0.7)

        # Add sample labels
        for i, sample in enumerate(self.sample_names):
            plt.annotate(
                sample,
                (pca_result[i, 0], pca_result[i, 1]),
                xytext=(5, 5),
                textcoords="offset points",
            )

        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
        plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
        plt.title("PCA of K-mer Profiles")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"PCA plot saved to {output_path}")

    def plot_mds(self, mds_result: np.ndarray, output_path: str):
        """Plot MDS results."""
        plt.figure(figsize=(10, 8))

        plt.scatter(mds_result[:, 0], mds_result[:, 1], s=100, alpha=0.7)

        # Add sample labels
        for i, sample in enumerate(self.sample_names):
            plt.annotate(
                sample,
                (mds_result[i, 0], mds_result[i, 1]),
                xytext=(5, 5),
                textcoords="offset points",
            )

        plt.xlabel("MDS Dimension 1")
        plt.ylabel("MDS Dimension 2")
        plt.title("MDS of Sample Distances")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"MDS plot saved to {output_path}")


def setup_logging(verbose: bool = False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler("metagrouper.log"),
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


def save_results(
    profiles: Dict,
    distance_matrix: np.ndarray,
    sample_names: List[str],
    output_dir: str,
):
    """Save analysis results."""
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


def main():
    parser = argparse.ArgumentParser(
        description="MetaGrouper: K-mer profiling and metadata-driven assembly grouping"
    )
    parser.add_argument("input_dir", help="Directory containing FASTQ files")
    parser.add_argument(
        "-o",
        "--output",
        default="metagrouper_output",
        help="Output directory (default: metagrouper_output)",
    )
    parser.add_argument(
        "-k", "--kmer-size", type=int, default=21, help="K-mer size (default: 21)"
    )
    parser.add_argument(
        "--max-reads", type=int, help="Maximum reads per sample (for testing)"
    )
    parser.add_argument(
        "--distance-metric",
        default="braycurtis",
        choices=["braycurtis", "jaccard", "cosine", "euclidean"],
        help="Distance metric (default: braycurtis)",
    )

    # Phase 2 arguments
    parser.add_argument(
        "-m", "--metadata", help="Metadata file (CSV/TSV) for Phase 2 analysis"
    )
    parser.add_argument(
        "--sample-id-column",
        default="sample_id",
        help="Column name for sample IDs in metadata (default: sample_id)",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        help="Specific metadata variables to analyze (default: all)",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=999,
        help="Number of permutations for PERMANOVA (default: 999)",
    )
    parser.add_argument(
        "--cluster-range",
        nargs=2,
        type=int,
        default=[2, 8],
        help="Range for number of clusters to test (default: 2 8)",
    )

    # Phase 3 arguments
    parser.add_argument(
        "--assembly-tools",
        nargs="+",
        choices=["megahit", "spades", "flye", "all"],
        default=["megahit", "spades"],
        help="Assembly tools to generate commands for (default: megahit spades)",
    )
    parser.add_argument(
        "--similarity-threshold",
        type=float,
        default=0.30,
        help="Distance threshold for grouping samples (default: 0.30)",
    )
    parser.add_argument(
        "--min-group-size",
        type=int,
        default=2,
        help="Minimum samples per assembly group (default: 2)",
    )
    parser.add_argument(
        "--max-group-size",
        type=int,
        default=10,
        help="Maximum samples per assembly group (default: 10)",
    )

    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    setup_logging(args.verbose)

    # Find FASTQ files
    logging.info(f"Searching for FASTQ files in {args.input_dir}")
    fastq_files = find_fastq_files(args.input_dir)

    if not fastq_files:
        logging.error(f"No FASTQ files found in {args.input_dir}")
        sys.exit(1)

    logging.info(f"Found {len(fastq_files)} FASTQ files")

    # Initialize profiler
    profiler = KmerProfiler(k=args.kmer_size, max_reads=args.max_reads)

    # Process each sample
    failed_samples = []
    for filepath, sample_name in fastq_files:
        try:
            profiler.profile_sample(filepath, sample_name)
        except ValueError as e:
            logging.error(f"Input validation error for {sample_name}: {e}")
            failed_samples.append(sample_name)
            continue
        except FileNotFoundError as e:
            logging.error(f"File not found for {sample_name}: {e}")
            failed_samples.append(sample_name)
            continue
        except PermissionError as e:
            logging.error(f"Permission denied for {sample_name}: {e}")
            failed_samples.append(sample_name)
            continue
        except Exception as e:
            logging.error(f"Unexpected error processing {sample_name}: {e}")
            logging.debug(f"Full error details: {type(e).__name__}: {e}")
            failed_samples.append(sample_name)
            continue
    
    if failed_samples:
        logging.warning(f"Failed to process {len(failed_samples)} samples: {', '.join(failed_samples)}")
    
    if not profiler.profiles:
        logging.error("No samples were successfully processed. Check input files and error messages above.")
        sys.exit(1)

    logging.info(f"Successfully processed {len(profiler.profiles)} samples")

    # Analyze similarities
    analyzer = SimilarityAnalyzer(profiler.profiles)
    distance_matrix = analyzer.compute_distance_matrix(metric=args.distance_metric)

    # Perform dimensionality reduction
    pca_result, pca = analyzer.perform_pca()
    mds_result = analyzer.perform_mds()

    # Generate visualizations
    visualizer = Visualizer(profiler.sample_names)
    output_path = Path(args.output)
    output_path.mkdir(exist_ok=True)

    visualizer.plot_distance_heatmap(
        distance_matrix, output_path / "distance_heatmap.png"
    )
    visualizer.plot_pca(pca_result, pca, output_path / "pca_plot.png")
    visualizer.plot_mds(mds_result, output_path / "mds_plot.png")

    # Save results
    save_results(profiler.profiles, distance_matrix, profiler.sample_names, args.output)

    # Phase 2: Metadata Analysis
    metadata_results_df = None
    cluster_results = {}

    if args.metadata and PHASE2_AVAILABLE:
        print(f"\nStarting Phase 2: Metadata Analysis...")

        try:
            # Initialize metadata analyzer
            meta_analyzer = MetadataAnalyzer(distance_matrix, profiler.sample_names)
            meta_analyzer.load_metadata(args.metadata, args.sample_id_column)

            # Analyze variables
            metadata_results_df = meta_analyzer.analyze_variables(
                variables=args.variables, n_permutations=args.permutations
            )

            # Identify clusters
            cluster_results = meta_analyzer.identify_clusters(
                n_clusters_range=tuple(args.cluster_range)
            )

            # Generate Phase 2 visualizations
            meta_visualizer = MetadataVisualizer(
                profiler.sample_names, meta_analyzer.metadata
            )

            # Variable importance plot
            meta_visualizer.plot_variable_importance(
                metadata_results_df, output_path / "variable_importance.png"
            )

            # Plot PCA colored by top variables
            if not metadata_results_df.empty:
                valid_results = metadata_results_df.dropna(subset=["r_squared"])
                top_variables = valid_results.head(3)["variable"].tolist()

                for var in top_variables:
                    safe_var_name = var.replace(" ", "_").replace("/", "_")
                    meta_visualizer.plot_samples_by_variable(
                        pca_result,
                        var,
                        output_path / f"pca_by_{safe_var_name}.png",
                        pca,
                    )

            # Plot clustering results
            for method, method_results in cluster_results.items():
                if "optimal" in method_results:
                    optimal = method_results["optimal"]
                    meta_visualizer.plot_clustering_results(
                        pca_result,
                        optimal["labels"],
                        method,
                        optimal["n_clusters"],
                        output_path / f"clustering_{method}.png",
                    )

            # Generate summary report
            generate_summary_report(
                metadata_results_df, cluster_results, output_path / "analysis_report.md"
            )

            # Save metadata analysis results
            if not metadata_results_df.empty:
                metadata_results_df.to_csv(
                    output_path / "permanova_results.csv", index=False
                )

            print(f"Phase 2 Analysis Complete!")

        except Exception as e:
            logging.error(f"Phase 2 analysis failed: {e}")
            print(f"Phase 2 analysis failed: {e}")

    elif args.metadata and not PHASE2_AVAILABLE:
        print(f"\nPhase 2 metadata analysis requested but not available.")
        print(f"Please install required dependencies.")

    # Phase 3: Assembly Strategy Recommendations
    assembly_recommendation = None

    if PHASE3_AVAILABLE:
        print(f"\nStarting Phase 3: Assembly Strategy Recommendations...")

        try:
            # Initialize assembly recommender
            recommender = AssemblyRecommender(distance_matrix, profiler.sample_names)

            # Configure thresholds based on arguments
            recommender.strategy_engine.similarity_threshold_medium = (
                args.similarity_threshold
            )
            recommender.strategy_engine.min_group_size = args.min_group_size
            recommender.strategy_engine.max_group_size = args.max_group_size

            # Generate recommendations
            if "all" in args.assembly_tools:
                tools = ["megahit", "spades", "flye"]
            else:
                tools = args.assembly_tools

            # Generate recommendations using both Phase 1 and Phase 2 results
            assembly_recommendation = recommender.generate_recommendations(
                metadata_results=metadata_results_df,
                metadata=(
                    meta_analyzer.metadata if "meta_analyzer" in locals() else None
                ),
            )

            # Filter assembly commands to requested tools
            filtered_commands = {
                tool: commands
                for tool, commands in assembly_recommendation.assembly_commands.items()
                if tool in tools
            }
            assembly_recommendation.assembly_commands = filtered_commands

            # Save recommendations
            save_recommendations(
                assembly_recommendation, output_path / "assembly_recommendations"
            )

            # Create visualization
            visualize_assembly_strategy(
                assembly_recommendation,
                distance_matrix,
                profiler.sample_names,
                output_path / "assembly_strategy_overview.png",
            )

            print(f"Phase 3 Analysis Complete!")
            print(f"Assembly Strategy: {assembly_recommendation.strategy.title()}")
            print(f"Confidence Score: {assembly_recommendation.overall_confidence:.2f}")
            print(f"Recommended Groups: {len(assembly_recommendation.groups)}")

        except Exception as e:
            logging.error(f"Phase 3 analysis failed: {e}")
            print(f"Phase 3 analysis failed: {e}")

    else:
        print(f"\nPhase 3 assembly recommendations not available.")
        print(f"Please install required dependencies.")

    # Print summary
    print(f"\nMetaGrouper Analysis Complete!")
    print(f"Processed {len(profiler.profiles)} samples")
    print(f"K-mer size: {args.kmer_size}")
    print(f"Distance metric: {args.distance_metric}")
    print(f"Results saved to: {args.output}")

    print(f"\nPhase 1 files:")
    print(f"  - distance_heatmap.png: Sample distance matrix visualization")
    print(f"  - pca_plot.png: PCA analysis of k-mer profiles")
    print(f"  - mds_plot.png: MDS analysis of sample distances")
    print(f"  - distance_matrix.csv: Pairwise distance matrix")
    print(f"  - kmer_profiles.pkl: K-mer profiles for all samples")

    if args.metadata and PHASE2_AVAILABLE and metadata_results_df is not None:
        print(f"\nPhase 2 files:")
        print(f"  - analysis_report.md: Comprehensive analysis report")
        print(f"  - variable_importance.png: Variable importance ranking")
        print(f"  - permanova_results.csv: PERMANOVA statistical results")
        print(f"  - pca_by_*.png: PCA plots colored by metadata variables")
        print(f"  - clustering_*.png: Clustering analysis results")

    if PHASE3_AVAILABLE and assembly_recommendation is not None:
        print(f"\nPhase 3 files:")
        print(f"  - assembly_recommendations/: Directory with detailed recommendations")
        print(f"  - assembly_strategy.md: Human-readable assembly strategy")
        print(f"  - assembly_recommendations.json: Detailed recommendations (JSON)")
        print(f"  - run_*_assemblies.sh: Executable assembly scripts")
        print(f"  - assembly_strategy_overview.png: Strategy visualization")


if __name__ == "__main__":
    main()
