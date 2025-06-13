#!/usr/bin/env python3
"""
MetaGrouper: K-mer-based Analysis for Optimal Metagenomic Assembly Grouping

This is the main entry point for MetaGrouper, using the modular architecture.
It handles command-line argument parsing and orchestrates the analysis workflow.
"""

import sys
import argparse
import logging
from pathlib import Path

# Add the package to the path for development
sys.path.insert(0, str(Path(__file__).parent))

from metagrouper import (
    KmerProfiler, 
    SimilarityAnalyzer, 
    Visualizer,
    find_fastq_files,
    setup_logging,
    save_results,
    MetaGrouperConfig
)

# Import Phase 2 and 3 functionality if available
try:
    from metadata_analyzer import (
        MetadataAnalyzer,
        MetadataVisualizer,
        generate_summary_report,
    )
    PHASE2_AVAILABLE = True
except ImportError:
    PHASE2_AVAILABLE = False
    logging.warning("Phase 2 metadata analysis not available. Install required dependencies.")

try:
    from assembly_recommender import (
        AssemblyRecommender,
        save_recommendations,
        visualize_assembly_strategy,
    )
    PHASE3_AVAILABLE = True
except ImportError:
    PHASE3_AVAILABLE = False
    logging.warning("Phase 3 assembly recommendations not available. Install required dependencies.")


def parse_arguments():
    """Parse command-line arguments."""
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
    
    # Configuration file
    parser.add_argument(
        "--config",
        help="Configuration file (JSON format)"
    )
    
    # Profiling arguments
    parser.add_argument(
        "-k", "--kmer-size", type=int, default=21, help="K-mer size (default: 21)"
    )
    parser.add_argument(
        "--max-reads", type=int, help="Maximum reads per sample (for testing)"
    )
    parser.add_argument(
        "--min-kmer-freq",
        type=int,
        default=1,
        help="Minimum k-mer frequency to keep (filters rare k-mers, default: 1)"
    )
    parser.add_argument(
        "--memory-efficient",
        action="store_true",
        default=True,
        help="Use memory-efficient processing (default: enabled)"
    )
    parser.add_argument(
        "--no-memory-efficient",
        action="store_true",
        help="Disable memory-efficient processing (load all data into memory)"
    )
    
    # Processing arguments
    parser.add_argument(
        "--processes",
        type=int,
        help="Number of parallel processes (default: CPU count)"
    )
    parser.add_argument(
        "--sequential",
        action="store_true",
        help="Disable parallel processing (use sequential mode)"
    )
    
    # Analysis arguments
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

    return parser.parse_args()


def run_phase1_analysis(config: MetaGrouperConfig, fastq_files: list) -> tuple:
    """Run Phase 1: K-mer profiling and similarity analysis."""
    print("\n🧬 Phase 1: K-mer Profiling and Similarity Analysis")
    print("=" * 60)
    
    # Initialize profiler
    profiler = KmerProfiler(
        k=config.profiling.k_size,
        max_reads=config.profiling.max_reads,
        min_kmer_freq=config.profiling.min_kmer_freq
    )
    
    # Process samples
    if config.processing.sequential or len(fastq_files) <= 2:
        logging.info("Using sequential processing")
        failed_samples = []
        for filepath, sample_name in fastq_files:
            try:
                profiler.profile_sample(
                    filepath, 
                    sample_name, 
                    memory_efficient=config.profiling.memory_efficient
                )
            except Exception as e:
                logging.error(f"Failed to process {sample_name}: {e}")
                failed_samples.append(sample_name)
    else:
        logging.info("Using parallel processing")
        try:
            profiles, failed_samples = profiler.process_samples_parallel(
                fastq_files,
                n_processes=config.processing.n_processes,
                show_progress=config.processing.show_progress,
                memory_efficient=config.profiling.memory_efficient
            )
        except Exception as e:
            logging.error(f"Parallel processing failed: {e}, falling back to sequential")
            failed_samples = []
            for filepath, sample_name in fastq_files:
                try:
                    profiler.profile_sample(
                        filepath, 
                        sample_name, 
                        memory_efficient=config.profiling.memory_efficient
                    )
                except Exception as e:
                    logging.error(f"Failed to process {sample_name}: {e}")
                    failed_samples.append(sample_name)
    
    if not profiler.profiles:
        raise RuntimeError("No samples were successfully processed")
    
    # Analyze similarities
    analyzer = SimilarityAnalyzer(
        profiler.profiles, 
        memory_efficient=config.analysis.memory_efficient
    )
    distance_matrix = analyzer.compute_distance_matrix(metric=config.analysis.distance_metric)
    
    # Perform dimensionality reduction
    pca_result, pca = analyzer.perform_pca()
    mds_result = analyzer.perform_mds()
    
    print(f"✅ Successfully processed {len(profiler.profiles)} samples")
    if failed_samples:
        print(f"⚠️  Failed to process {len(failed_samples)} samples")
    
    return profiler, analyzer, distance_matrix, pca_result, pca, mds_result, failed_samples


def run_phase2_analysis(config: MetaGrouperConfig, distance_matrix, sample_names, 
                       pca_result, pca, output_path, args):
    """Run Phase 2: Metadata analysis (if available and requested)."""
    if not args.metadata or not PHASE2_AVAILABLE:
        return None, {}
    
    print("\n📊 Phase 2: Metadata Analysis")
    print("=" * 60)
    
    try:
        # Initialize metadata analyzer
        meta_analyzer = MetadataAnalyzer(distance_matrix, sample_names)
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
        meta_visualizer = MetadataVisualizer(sample_names, meta_analyzer.metadata)
        
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
            metadata_results_df.to_csv(output_path / "permanova_results.csv", index=False)
        
        print("✅ Phase 2 analysis completed successfully")
        return metadata_results_df, cluster_results
        
    except Exception as e:
        logging.error(f"Phase 2 analysis failed: {e}")
        print(f"❌ Phase 2 analysis failed: {e}")
        return None, {}


def run_phase3_analysis(config: MetaGrouperConfig, distance_matrix, sample_names,
                       metadata_results_df, output_path, args):
    """Run Phase 3: Assembly recommendations (if available)."""
    if not PHASE3_AVAILABLE:
        return None
    
    print("\n🔧 Phase 3: Assembly Strategy Recommendations")
    print("=" * 60)
    
    try:
        # Initialize assembly recommender
        recommender = AssemblyRecommender(distance_matrix, sample_names)
        
        # Configure thresholds
        recommender.strategy_engine.similarity_threshold_medium = args.similarity_threshold
        recommender.strategy_engine.min_group_size = args.min_group_size
        recommender.strategy_engine.max_group_size = args.max_group_size
        
        # Generate recommendations
        tools = ["megahit", "spades", "flye"] if "all" in args.assembly_tools else args.assembly_tools
        
        assembly_recommendation = recommender.generate_recommendations(
            metadata_results=metadata_results_df,
            metadata=getattr(locals().get('meta_analyzer'), 'metadata', None)
        )
        
        # Filter assembly commands
        filtered_commands = {
            tool: commands
            for tool, commands in assembly_recommendation.assembly_commands.items()
            if tool in tools
        }
        assembly_recommendation.assembly_commands = filtered_commands
        
        # Save recommendations
        save_recommendations(assembly_recommendation, output_path / "assembly_recommendations")
        
        # Create visualization
        visualize_assembly_strategy(
            assembly_recommendation,
            distance_matrix,
            sample_names,
            output_path / "assembly_strategy_overview.png",
        )
        
        print(f"✅ Phase 3 analysis completed successfully")
        print(f"   Strategy: {assembly_recommendation.strategy.title()}")
        print(f"   Confidence: {assembly_recommendation.overall_confidence:.2f}")
        print(f"   Groups: {len(assembly_recommendation.groups)}")
        
        return assembly_recommendation
        
    except Exception as e:
        logging.error(f"Phase 3 analysis failed: {e}")
        print(f"❌ Phase 3 analysis failed: {e}")
        return None


def generate_visualizations(visualizer, distance_matrix, pca_result, pca, mds_result, output_path):
    """Generate all Phase 1 visualizations."""
    print("\n📈 Generating visualizations...")
    
    visualizer.plot_distance_heatmap(distance_matrix, output_path / "distance_heatmap.png")
    visualizer.plot_pca(pca_result, pca, output_path / "pca_plot.png")
    visualizer.plot_mds(mds_result, output_path / "mds_plot.png")
    
    # Generate additional comprehensive plots
    visualizer.plot_sample_overview(distance_matrix, pca_result, mds_result, pca, 
                                   output_path / "sample_overview.png")
    visualizer.plot_clustering_dendrogram(distance_matrix, output_path / "dendrogram.png")
    visualizer.plot_distance_distribution(distance_matrix, output_path / "distance_distribution.png")


def print_summary(config: MetaGrouperConfig, profiler, failed_samples, 
                 metadata_results_df, assembly_recommendation):
    """Print analysis summary."""
    print("\n" + "="*60)
    print("🎉 MetaGrouper Analysis Complete!")
    print("="*60)
    
    print(f"📊 Processed: {len(profiler.profiles)} samples")
    print(f"🧬 K-mer size: {config.profiling.k_size}")
    print(f"📏 Distance metric: {config.analysis.distance_metric}")
    print(f"📁 Results saved to: {config.output.output_dir}")
    
    if failed_samples:
        print(f"⚠️  Failed samples: {len(failed_samples)}")
    
    print(f"\n📋 Phase 1 Files:")
    print(f"   • distance_heatmap.png: Sample distance matrix")
    print(f"   • pca_plot.png: PCA analysis")
    print(f"   • mds_plot.png: MDS analysis")
    print(f"   • sample_overview.png: Comprehensive overview")
    print(f"   • dendrogram.png: Hierarchical clustering")
    print(f"   • distance_distribution.png: Distance statistics")
    
    if metadata_results_df is not None and not metadata_results_df.empty:
        print(f"\n📋 Phase 2 Files:")
        print(f"   • analysis_report.md: Comprehensive report")
        print(f"   • variable_importance.png: Variable ranking")
        print(f"   • permanova_results.csv: Statistical results")
        print(f"   • pca_by_*.png: PCA colored by variables")
        print(f"   • clustering_*.png: Clustering results")
    
    if assembly_recommendation is not None:
        print(f"\n📋 Phase 3 Files:")
        print(f"   • assembly_recommendations/: Detailed recommendations")
        print(f"   • assembly_strategy.md: Human-readable strategy")
        print(f"   • run_*_assemblies.sh: Executable scripts")
        print(f"   • assembly_strategy_overview.png: Strategy visualization")


def main():
    """Main analysis workflow."""
    args = parse_arguments()
    
    # Initialize configuration
    config = MetaGrouperConfig(args.config if args.config else None)
    config.update_from_args(args)
    config.validate()
    
    # Setup logging
    setup_logging(config.output.verbose, config.output.log_file)
    
    # Print configuration summary
    if config.output.verbose:
        print(config.get_summary())
        print()
    
    # Find and validate input files
    logging.info(f"Searching for FASTQ files in {args.input_dir}")
    fastq_files = find_fastq_files(args.input_dir)
    
    if not fastq_files:
        logging.error(f"No FASTQ files found in {args.input_dir}")
        sys.exit(1)
    
    logging.info(f"Found {len(fastq_files)} samples")
    
    # Create output directory
    output_path = Path(config.output.output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Save configuration
    config.save_to_file(output_path / "analysis_config.json")
    
    try:
        # Phase 1: K-mer profiling and similarity analysis
        profiler, analyzer, distance_matrix, pca_result, pca, mds_result, failed_samples = \
            run_phase1_analysis(config, fastq_files)
        
        # Generate Phase 1 visualizations
        visualizer = Visualizer(profiler.sample_names)
        generate_visualizations(visualizer, distance_matrix, pca_result, pca, mds_result, output_path)
        
        # Save Phase 1 results
        save_results(profiler.profiles, distance_matrix, profiler.sample_names, config.output.output_dir)
        
        # Phase 2: Metadata analysis (optional)
        metadata_results_df, cluster_results = run_phase2_analysis(
            config, distance_matrix, profiler.sample_names, pca_result, pca, output_path, args
        )
        
        # Phase 3: Assembly recommendations (optional)
        assembly_recommendation = run_phase3_analysis(
            config, distance_matrix, profiler.sample_names, metadata_results_df, output_path, args
        )
        
        # Print summary
        print_summary(config, profiler, failed_samples, metadata_results_df, assembly_recommendation)
        
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        print(f"\n❌ Analysis failed: {e}")
        if config.output.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()