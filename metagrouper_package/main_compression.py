#!/usr/bin/env python3
"""
MetaGrouper with Compression-based Similarity Analysis

Main entry point for MetaGrouper using compression-based similarity measures
instead of k-mer profiling. This should be significantly faster for large datasets.
"""

import sys
import argparse
import logging
import time
from pathlib import Path
from typing import List, Tuple

# Add current directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from metagrouper.compression_profiler import CompressionProfiler
from metagrouper.compression_analyzer import CompressionAnalyzer
from metagrouper.visualizer import Visualizer
from metagrouper.utils import find_fastq_files, setup_logging, save_results
from metagrouper.config import MetaGrouperConfig
from metadata_analyzer import MetadataAnalyzer
from assembly_recommender import AssemblyRecommender


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='MetaGrouper: Compression-based Analysis for Optimal Metagenomic Assembly Grouping',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic compression-based analysis
  python main_compression.py samples/ -o results/
  
  # With metadata integration
  python main_compression.py samples/ --metadata metadata.csv -o results/
  
  # Fast analysis for large datasets
  python main_compression.py samples/ -o results/ --compressor zstd --max-reads 5000
  
  # Detailed analysis with assembly recommendations
  python main_compression.py samples/ --metadata metadata.csv -o results/ \\
      --similarity-threshold 0.3 --generate-commands
        """
    )
    
    # Input/Output arguments
    parser.add_argument('input_dir', help='Directory containing FASTQ files')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('--metadata', help='Metadata CSV file (optional)')
    
    # Compression arguments
    parser.add_argument('--compressor', choices=['gzip', 'zstd', 'lzma', 'bz2'], 
                       default='zstd', help='Compression algorithm (default: zstd)')
    parser.add_argument('--compression-level', type=int, default=3,
                       help='Compression level (1-9, default: 3)')
    parser.add_argument('--max-reads', type=int, 
                       help='Maximum reads per sample (for speed)')
    
    # Processing arguments
    parser.add_argument('--processes', type=int, 
                       help='Number of parallel processes (default: CPU count)')
    parser.add_argument('--sequential', action='store_true',
                       help='Use sequential processing instead of parallel')
    
    # Analysis arguments
    parser.add_argument('--similarity-threshold', type=float, default=0.25,
                       help='Similarity threshold for grouping (default: 0.25)')
    parser.add_argument('--permutations', type=int, default=999,
                       help='Number of PERMANOVA permutations (default: 999)')
    
    # Output arguments
    parser.add_argument('--generate-commands', action='store_true',
                       help='Generate assembly commands for recommended groups')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')
    parser.add_argument('--no-progress', action='store_true',
                       help='Disable progress updates during processing')
    
    return parser.parse_args()


def main():
    """Main analysis pipeline using compression-based similarity."""
    args = parse_arguments()
    
    # Setup logging
    setup_logging(verbose=args.verbose, log_file=Path(args.output) / "compression_analysis.log")
    
    logging.info("🗜️  Starting MetaGrouper compression-based analysis")
    logging.info(f"Using {args.compressor} compression with level {args.compression_level}")
    
    start_time = time.time()
    
    try:
        # Create output directory
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Find FASTQ files
        logging.info(f"Scanning for FASTQ files in {args.input_dir}")
        fastq_files = find_fastq_files(args.input_dir)
        
        if not fastq_files:
            logging.error(f"No FASTQ files found in {args.input_dir}")
            return 1
        
        logging.info(f"Found {len(fastq_files)} samples")
        
        # Initialize compression profiler
        profiler = CompressionProfiler(
            compressor=args.compressor,
            max_reads=args.max_reads,
            compression_level=args.compression_level
        )
        
        # Process samples
        if args.sequential or len(fastq_files) <= 2:
            logging.info("Using sequential processing")
            failed_samples = []
            for filepath, sample_name in fastq_files:
                try:
                    profiler.profile_sample(filepath, sample_name)
                except Exception as e:
                    logging.error(f"Failed to process {sample_name}: {e}")
                    failed_samples.append(sample_name)
        else:
            logging.info("Using parallel processing")
            try:
                compression_sizes, failed_samples = profiler.process_samples_parallel(
                    fastq_files,
                    n_processes=args.processes,
                    show_progress=not args.no_progress
                )
            except Exception as e:
                logging.error(f"Parallel processing failed: {e}, falling back to sequential")
                failed_samples = []
                for filepath, sample_name in fastq_files:
                    try:
                        profiler.profile_sample(filepath, sample_name)
                    except Exception as e:
                        logging.error(f"Failed to process {sample_name}: {e}")
                        failed_samples.append(sample_name)
        
        if not profiler.sample_names:
            logging.error("No samples were successfully processed")
            return 1
        
        if failed_samples:
            logging.warning(f"Failed to process {len(failed_samples)} samples: {failed_samples}")
        
        # Compute compression-based distance matrix
        logging.info("Computing compression-based distance matrix")
        distance_matrix = profiler.compute_distance_matrix()
        
        # Initialize analyzer
        analyzer = CompressionAnalyzer(distance_matrix, profiler.sample_names)
        
        # Compute summary statistics
        stats = analyzer.compute_summary_statistics()
        
        # Perform dimensionality reduction
        logging.info("Performing dimensionality reduction")
        pca_result, pca_model = analyzer.perform_pca()
        mds_result = analyzer.perform_mds()
        
        # Perform clustering
        logging.info("Performing clustering analysis")
        cluster_results = analyzer.perform_clustering()
        
        # Identify outliers
        outliers = analyzer.identify_outliers()
        if outliers:
            logging.info(f"Identified outlier samples: {outliers}")
        
        # Save compression-based results
        save_compression_results(
            analyzer, profiler, stats, cluster_results, outliers, output_dir
        )
        
        # Metadata analysis (if provided)
        metadata_results = None
        if args.metadata:
            logging.info(f"Performing metadata analysis using {args.metadata}")
            
            meta_analyzer = MetadataAnalyzer(distance_matrix, profiler.sample_names)
            meta_analyzer.load_metadata(args.metadata)
            
            metadata_results = meta_analyzer.analyze_variables(
                n_permutations=args.permutations
            )
            
            # Save metadata results
            metadata_results.to_csv(output_dir / "metadata_analysis.csv", index=False)
            logging.info("Metadata analysis results saved")
        
        # Assembly recommendations (if requested)
        if args.generate_commands:
            logging.info("Generating assembly recommendations")
            
            recommender = AssemblyRecommender(distance_matrix, profiler.sample_names)
            recommendation = recommender.generate_recommendations(
                metadata_results=metadata_results,
                metadata=meta_analyzer.metadata if args.metadata else None,
                similarity_threshold=args.similarity_threshold
            )
            
            # Save assembly recommendations
            save_assembly_recommendations(recommendation, output_dir)
        
        # Generate visualizations
        logging.info("Generating visualizations")
        create_compression_visualizations(analyzer, profiler, output_dir)
        
        # Final summary
        elapsed_time = time.time() - start_time
        logging.info(f"🎉 Compression-based analysis completed in {elapsed_time:.1f}s")
        logging.info(f"Results saved to {output_dir}")
        logging.info(f"Processed {len(profiler.sample_names)} samples using {args.compressor} compression")
        
        return 0
        
    except KeyboardInterrupt:
        logging.warning("Analysis interrupted by user")
        return 1
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


def save_compression_results(analyzer, profiler, stats, cluster_results, outliers, output_dir):
    """Save compression-based analysis results."""
    
    # Save distance matrix
    analyzer.export_distance_matrix(output_dir / "compression_distance_matrix.csv")
    
    # Save compression statistics
    import json
    with open(output_dir / "compression_stats.json", 'w') as f:
        json.dump(stats, f, indent=2)
    
    # Save compression sizes
    import pandas as pd
    compression_df = pd.DataFrame([
        {'sample_name': name, 'compressed_size': size}
        for name, size in profiler.compression_sizes.items()
    ])
    compression_df.to_csv(output_dir / "compression_sizes.csv", index=False)
    
    # Save clustering results
    for method, results in cluster_results.items():
        if 'optimal' in results:
            optimal = results['optimal']
            cluster_df = pd.DataFrame({
                'sample_name': analyzer.sample_names,
                'cluster': optimal['labels']
            })
            cluster_df.to_csv(output_dir / f"clusters_{method}.csv", index=False)
    
    # Save outliers
    if outliers:
        outlier_df = pd.DataFrame({'outlier_samples': outliers})
        outlier_df.to_csv(output_dir / "outliers.csv", index=False)
    
    # Save PCA and MDS coordinates
    if analyzer.pca_result is not None:
        pca_df = pd.DataFrame(
            analyzer.pca_result,
            columns=[f'PC{i+1}' for i in range(analyzer.pca_result.shape[1])],
            index=analyzer.sample_names
        )
        pca_df.to_csv(output_dir / "pca_coordinates.csv")
    
    if analyzer.mds_result is not None:
        mds_df = pd.DataFrame(
            analyzer.mds_result,
            columns=[f'MDS{i+1}' for i in range(analyzer.mds_result.shape[1])],
            index=analyzer.sample_names
        )
        mds_df.to_csv(output_dir / "mds_coordinates.csv")


def save_assembly_recommendations(recommendation, output_dir):
    """Save assembly recommendations to files."""
    
    # Save assembly commands
    import json
    with open(output_dir / "assembly_commands.json", 'w') as f:
        json.dump(recommendation.assembly_commands, f, indent=2)
    
    # Save recommendation summary
    summary = {
        'strategy': recommendation.strategy,
        'overall_confidence': recommendation.overall_confidence,
        'groups': [
            {
                'group_id': i,
                'sample_names': group.sample_names,
                'avg_distance': group.avg_distance,
                'confidence': group.confidence
            }
            for i, group in enumerate(recommendation.groups)
        ]
    }
    
    with open(output_dir / "assembly_recommendation.json", 'w') as f:
        json.dump(summary, f, indent=2)


def create_compression_visualizations(analyzer, profiler, output_dir):
    """Create visualizations for compression-based analysis."""
    
    # Initialize visualizer
    visualizer = Visualizer(analyzer.sample_names)
    
    # Distance matrix heatmap
    heatmap_path = output_dir / "compression_distance_heatmap.png"
    visualizer.plot_distance_heatmap(analyzer.distance_matrix, str(heatmap_path))
    
    # PCA plot
    if analyzer.pca_result is not None:
        pca_path = output_dir / "compression_pca.png"
        visualizer.plot_pca(analyzer.pca_result, str(pca_path))
    
    # MDS plot
    if analyzer.mds_result is not None:
        mds_path = output_dir / "compression_mds.png"
        visualizer.plot_mds(analyzer.mds_result, str(mds_path))
    
    logging.info("Visualizations saved")


if __name__ == "__main__":
    sys.exit(main())