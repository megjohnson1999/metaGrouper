#!/usr/bin/env python3
"""
MetaGrouper: K-mer-based Analysis for Optimal Metagenomic Assembly Grouping

This is the main entry point for MetaGrouper, implementing all three phases:
Phase 1: K-mer profiling and similarity analysis
Phase 2: Metadata analysis and variable testing
Phase 3: Assembly strategy recommendations
"""

import sys
import argparse
import logging
import time
from pathlib import Path
import numpy as np
import pandas as pd

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))

from metagrouper import (
    find_fastq_files,
    setup_logging,
    save_results,
    KmerProfiler,
    SimilarityAnalyzer,
    Visualizer
)

# Phase 1 imports
try:
    from metagrouper.sketch_profiler import StreamingKmerProfiler
    from metagrouper.sparse_analyzer import SparseSimilarityAnalyzer
    PHASE1_AVAILABLE = True
except ImportError as e:
    PHASE1_AVAILABLE = False
    logging.error(f"Phase 1 not available: {e}")

# Phase 2 imports
try:
    from metadata_analyzer import (
        MetadataAnalyzer,
        MetadataVisualizer,
        generate_summary_report,
    )
    PHASE2_AVAILABLE = True
except ImportError as e:
    PHASE2_AVAILABLE = False
    logging.warning(f"Phase 2 metadata analysis not available: {e}")

# Phase 3 imports
try:
    from assembly_recommender import (
        AssemblyRecommender,
        save_recommendations,
        visualize_assembly_strategy,
    )
    PHASE3_AVAILABLE = True
except ImportError as e:
    PHASE3_AVAILABLE = False
    logging.warning(f"Phase 3 assembly recommendations not available: {e}")


def create_parser():
    """Create command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="MetaGrouper: K-mer profiling and metadata-driven assembly grouping",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis (Phase 1 only)
  python metagrouper.py /path/to/fastq/files -o results/
  
  # Full analysis with metadata (All phases)
  python metagrouper.py /path/to/fastq/files -m metadata.csv -o results/
  
  # Memory-efficient for large datasets
  python metagrouper.py /path/to/fastq/files -o results/ --use-sketching --sketch-size 5000
  
  # Complete workflow with assembly recommendations
  python metagrouper.py /path/to/fastq/files --metadata samples_metadata.csv \\
    --output results/ --assembly-tools megahit spades --similarity-threshold 0.25
        """
    )
    
    # Required arguments
    parser.add_argument("input_dir", help="Directory containing FASTQ files")
    parser.add_argument("-o", "--output", default="metagrouper_output", 
                       help="Output directory (default: metagrouper_output)")
    
    # K-mer Analysis (Phase 1)
    parser.add_argument("-k", "--kmer-size", type=int, default=21, 
                       help="K-mer size (default: 21)")
    parser.add_argument("--max-reads", type=int, 
                       help="Maximum reads per sample (for testing)")
    parser.add_argument("--distance-metric", default="braycurtis",
                       choices=["braycurtis", "jaccard", "cosine", "euclidean"],
                       help="Distance metric (default: braycurtis)")
    
    # Phase 1 optimization arguments
    parser.add_argument("--use-sketching", action="store_true",
                       help="Use streaming k-mer sketches for memory efficiency")
    parser.add_argument("--sketch-size", type=int, default=1000,
                       help="K-mer sketch size (default: 1000)")
    parser.add_argument("--sampling-method", choices=["reservoir", "frequency", "adaptive"],
                       default="frequency", help="Sketching method (default: frequency)")
    
    # Metadata Analysis (Phase 2)
    parser.add_argument("-m", "--metadata", help="Metadata file (CSV/TSV) for Phase 2 analysis")
    parser.add_argument("--sample-id-column", default="sample_id",
                       help="Column name for sample IDs in metadata (default: sample_id)")
    parser.add_argument("--variables", nargs="+",
                       help="Specific metadata variables to analyze (default: all)")
    parser.add_argument("--permutations", type=int, default=999,
                       help="Number of permutations for PERMANOVA (default: 999)")
    parser.add_argument("--cluster-range", nargs=2, type=int, default=[2, 8],
                       help="Range for number of clusters to test (default: 2 8)")
    
    # Assembly Recommendations (Phase 3)
    parser.add_argument("--assembly-tools", nargs="+",
                       choices=["megahit", "spades", "flye", "all"],
                       default=["megahit", "spades"],
                       help="Assembly tools to generate commands for (default: megahit spades)")
    parser.add_argument("--similarity-threshold", type=float, default=0.30,
                       help="Distance threshold for grouping samples (default: 0.30)")
    parser.add_argument("--min-group-size", type=int, default=2,
                       help="Minimum samples per assembly group (default: 2)")
    parser.add_argument("--max-group-size", type=int, default=10,
                       help="Maximum samples per assembly group (default: 10)")
    
    # Processing arguments
    parser.add_argument("--processes", type=int,
                       help="Number of parallel processes (default: CPU count)")
    parser.add_argument("--sequential", action="store_true",
                       help="Use sequential processing")
    
    # Other arguments
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    
    return parser


def run_analysis(args):
    """Run the main analysis workflow."""
    # Setup logging
    setup_logging(args.verbose)
    
    print("🧬 MetaGrouper: K-mer Analysis for Optimal Metagenomic Assembly Grouping")
    print("=" * 80)
    
    # Validate input
    if not Path(args.input_dir).exists():
        logging.error(f"Input directory not found: {args.input_dir}")
        return False
    
    # Find FASTQ files
    logging.info("Finding FASTQ files...")
    fastq_files = find_fastq_files(args.input_dir)
    if not fastq_files:
        logging.error("No FASTQ files found")
        return False
    
    logging.info(f"Found {len(fastq_files)} samples")
    
    # Create output directory
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Determine which phases to run
    run_phase2 = args.metadata and PHASE2_AVAILABLE
    run_phase3 = PHASE3_AVAILABLE  # Phase 3 can run without metadata
    
    print(f"📋 Analysis Plan:")
    print(f"   Phase 1: K-mer profiling and similarity ✅")
    print(f"   Phase 2: Metadata analysis {'✅' if run_phase2 else '❌ (no metadata provided)' if not args.metadata else '❌ (dependencies missing)'}")
    print(f"   Phase 3: Assembly recommendations {'✅' if run_phase3 else '❌ (dependencies missing)'}")
    print()
    
    # =============================================================================
    # PHASE 1: K-mer Profiling and Similarity Analysis
    # =============================================================================
    print("🧬 Phase 1: K-mer Profiling and Similarity Analysis")
    print("-" * 60)
    
    # Choose profiler based on options
    if args.use_sketching and PHASE1_AVAILABLE:
        print(f"🚀 Using optimized sketching approach")
        print(f"   Sketch size: {args.sketch_size}")
        print(f"   Sampling method: {args.sampling_method}")
        
        profiler = StreamingKmerProfiler(
            k=args.kmer_size,
            max_reads=args.max_reads,
            sketch_size=args.sketch_size,
            sampling_method=args.sampling_method
        )
    else:
        print(f"📊 Using traditional k-mer profiling")
        profiler = KmerProfiler(
            k=args.kmer_size,
            max_reads=args.max_reads
        )
    
    # Process samples
    print(f"🔬 Processing {len(fastq_files)} samples...")
    start_time = time.time()
    
    profiles = {}
    failed_samples = []
    
    for i, (filepath, sample_name) in enumerate(fastq_files):
        if i % 10 == 0 or i == len(fastq_files) - 1:
            print(f"   Progress: {i+1}/{len(fastq_files)} samples")
        
        try:
            profile = profiler.profile_sample(filepath, sample_name)
            profiles[sample_name] = profile
        except Exception as e:
            logging.error(f"Failed to process {sample_name}: {e}")
            failed_samples.append(sample_name)
    
    processing_time = time.time() - start_time
    success_count = len(profiles)
    
    print(f"✅ Processed {success_count}/{len(fastq_files)} samples in {processing_time:.1f}s")
    if failed_samples:
        print(f"❌ Failed samples: {', '.join(failed_samples)}")
    
    if not profiles:
        logging.error("No samples processed successfully")
        return False
    
    # Memory usage report for sketching
    if isinstance(profiler, StreamingKmerProfiler):
        memory_info = profiler.estimate_memory_usage()
        print(f"💾 Memory usage: {memory_info['sketch_memory_mb']:.1f} MB")
        print(f"💾 Memory reduction: {memory_info['memory_reduction']:.1f}x vs full profiles")
    
    # Compute similarity/distance matrix
    print(f"\n🔗 Computing similarity matrix...")
    start_time = time.time()
    
    if args.use_sketching and PHASE1_AVAILABLE:
        # Use sparse similarity analysis
        sparse_threshold = getattr(args, 'sparse_threshold', 0.1)  # Default sparse threshold
        analyzer = SparseSimilarityAnalyzer(similarity_threshold=sparse_threshold)
        similarity_matrix, sample_names = analyzer.compute_similarities(profiles, method=args.distance_metric)
        
        # Convert sparse to dense for downstream analysis
        distance_matrix = 1 - similarity_matrix.toarray()
        
        # Get sparse statistics
        sparse_stats = analyzer.compute_summary_statistics()
        similarity_time = time.time() - start_time
        
        print(f"✅ Sparse similarity computed in {similarity_time:.1f}s")
        print(f"📊 Sparsity: {sparse_stats['sparsity']:.1%}")
        print(f"📊 Significant pairs: {sparse_stats['n_significant_pairs']:,}")
        
        # Save sparse matrix
        analyzer.save_similarity_matrix(str(output_path / "similarity_matrix.npz"), format='npz')
        
    else:
        # Traditional dense similarity 
        analyzer = SimilarityAnalyzer(profiles)
        distance_matrix = analyzer.compute_distance_matrix(metric=args.distance_metric)
        sample_names = list(profiles.keys())
        
        similarity_time = time.time() - start_time
        print(f"✅ Dense similarity computed in {similarity_time:.1f}s")
        
        # Save traditional results
        save_results(profiles, distance_matrix, sample_names, args.output)
    
    # Perform dimensionality reduction for visualization
    pca_result, pca = analyzer.perform_pca() if hasattr(analyzer, 'perform_pca') else (None, None)
    if pca_result is None:
        # Fallback PCA using sklearn
        try:
            from sklearn.decomposition import PCA
            from sklearn.manifold import MDS
            
            pca = PCA(n_components=min(2, len(sample_names)-1))
            if args.use_sketching:
                # Use original profiles for PCA
                profile_matrix = np.array([list(profiles[name].values()) for name in sample_names])
            else:
                profile_matrix = 1 - distance_matrix
            pca_result = pca.fit_transform(profile_matrix)
        except Exception as e:
            logging.warning(f"PCA failed: {e}")
            pca_result, pca = None, None
    
    # Generate Phase 1 visualizations
    print(f"\n📈 Generating Phase 1 visualizations...")
    try:
        visualizer = Visualizer(sample_names)
        
        # Basic distance heatmap
        if len(sample_names) <= 100:  # Only for reasonably sized matrices
            visualizer.plot_distance_heatmap(distance_matrix, output_path / "distance_heatmap.png")
        
        # PCA plot
        if pca_result is not None:
            visualizer.plot_pca(pca_result, pca, output_path / "pca_plot.png")
        
        print(f"✅ Phase 1 visualizations saved")
        
    except Exception as e:
        logging.warning(f"Phase 1 visualization failed: {e}")
    
    # =============================================================================
    # PHASE 2: Metadata Analysis (if requested)
    # =============================================================================
    metadata_results_df = pd.DataFrame()
    cluster_results = {}
    
    if run_phase2:
        print(f"\n📊 Phase 2: Metadata Analysis")
        print("-" * 60)
        
        try:
            # Initialize metadata analyzer
            meta_analyzer = MetadataAnalyzer(distance_matrix, sample_names)
            meta_analyzer.load_metadata(args.metadata, args.sample_id_column)
            
            # Analyze variables
            metadata_results_df = meta_analyzer.analyze_variables(
                variables=args.variables, n_permutations=args.permutations
            )
            
            print(f"✅ Analyzed {len(metadata_results_df)} metadata variables")
            
            # Identify clusters
            cluster_results = meta_analyzer.identify_clusters(
                n_clusters_range=tuple(args.cluster_range)
            )
            
            print(f"✅ Clustering analysis complete")
            
            # Generate Phase 2 visualizations
            meta_visualizer = MetadataVisualizer(sample_names, meta_analyzer.metadata)
            
            # Variable importance plot
            if not metadata_results_df.empty:
                meta_visualizer.plot_variable_importance(
                    metadata_results_df, output_path / "variable_importance.png"
                )
            
            # Plot PCA colored by top variables
            if not metadata_results_df.empty and pca_result is not None:
                valid_results = metadata_results_df.dropna(subset=["r_squared"])
                top_variables = valid_results.head(3)["variable"].tolist()
                
                for var in top_variables:
                    safe_var_name = var.replace(" ", "_").replace("/", "_")
                    meta_visualizer.plot_samples_by_variable(
                        pca_result, var, output_path / f"pca_by_{safe_var_name}.png", pca
                    )
            
            # Plot clustering results
            for method, method_results in cluster_results.items():
                if "optimal" in method_results and pca_result is not None:
                    optimal = method_results["optimal"]
                    meta_visualizer.plot_clustering_results(
                        pca_result, optimal["labels"], method, optimal["n_clusters"],
                        output_path / f"clustering_{method}.png"
                    )
            
            # Generate summary report
            generate_summary_report(
                metadata_results_df, cluster_results, output_path / "analysis_report.md"
            )
            
            # Save metadata analysis results
            if not metadata_results_df.empty:
                metadata_results_df.to_csv(output_path / "permanova_results.csv", index=False)
            
            print(f"✅ Phase 2 analysis completed successfully")
            
        except Exception as e:
            logging.error(f"Phase 2 analysis failed: {e}")
            print(f"❌ Phase 2 analysis failed: {e}")
            run_phase2 = False
    
    # =============================================================================
    # PHASE 3: Assembly Strategy Recommendations
    # =============================================================================
    assembly_recommendation = None
    
    if run_phase3:
        print(f"\n🔧 Phase 3: Assembly Strategy Recommendations")
        print("-" * 60)
        
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
                metadata_results=metadata_results_df if run_phase2 else None,
                metadata=getattr(meta_analyzer, 'metadata', None) if run_phase2 else None
            )
            
            # Filter assembly commands to requested tools
            filtered_commands = {
                tool: commands
                for tool, commands in assembly_recommendation.assembly_commands.items()
                if tool in tools
            }
            assembly_recommendation.assembly_commands = filtered_commands
            
            # Save recommendations
            save_recommendations(assembly_recommendation, output_path / "assembly_recommendations")
            
            # Create visualization
            if pca_result is not None:
                visualize_assembly_strategy(
                    assembly_recommendation, distance_matrix, sample_names,
                    output_path / "assembly_strategy_overview.png"
                )
            
            print(f"✅ Phase 3 analysis completed successfully")
            print(f"   Strategy: {assembly_recommendation.strategy.title()}")
            print(f"   Confidence: {assembly_recommendation.overall_confidence:.1%}")
            print(f"   Groups: {len(assembly_recommendation.groups)}")
            
        except Exception as e:
            logging.error(f"Phase 3 analysis failed: {e}")
            print(f"❌ Phase 3 analysis failed: {e}")
    
    # =============================================================================
    # FINAL SUMMARY
    # =============================================================================
    total_time = processing_time + similarity_time
    print(f"\n" + "="*80)
    print(f"🎉 MetaGrouper Analysis Complete!")
    print(f"="*80)
    
    print(f"📊 Phase 1 Results:")
    print(f"   • Samples processed: {success_count}/{len(fastq_files)}")
    print(f"   • Processing time: {total_time:.1f}s")
    print(f"   • K-mer size: {args.kmer_size}")
    print(f"   • Distance metric: {args.distance_metric}")
    
    if args.use_sketching and isinstance(profiler, StreamingKmerProfiler):
        print(f"   • Memory efficiency: {memory_info['memory_reduction']:.1f}x improvement")
        print(f"   • Sparsity benefit: {sparse_stats['sparsity']:.1%} memory saved")
    
    if run_phase2:
        print(f"\n📊 Phase 2 Results:")
        significant_vars = metadata_results_df[metadata_results_df['p_value'] < 0.05] if not metadata_results_df.empty else pd.DataFrame()
        print(f"   • Variables analyzed: {len(metadata_results_df)}")
        print(f"   • Significant associations: {len(significant_vars)}")
        if len(significant_vars) > 0:
            top_var = significant_vars.iloc[0]
            print(f"   • Top variable: {top_var['variable']} (R² = {top_var['r_squared']:.3f})")
    
    if run_phase3 and assembly_recommendation:
        print(f"\n🔧 Phase 3 Results:")
        print(f"   • Assembly strategy: {assembly_recommendation.strategy.title()}")
        print(f"   • Confidence: {assembly_recommendation.overall_confidence:.1%}")
        print(f"   • Assembly groups: {len(assembly_recommendation.groups)}")
        print(f"   • Assembly tools: {', '.join(assembly_recommendation.assembly_commands.keys())}")
    
    print(f"\n📁 Results saved to: {args.output}")
    print(f"📚 Check the output directory for detailed results and visualizations")
    
    return True


def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    if not PHASE1_AVAILABLE and args.use_sketching:
        print("❌ Phase 1 optimizations not available but requested")
        print("   Install Phase 1 components or remove --use-sketching flag")
        return 1
    
    success = run_analysis(args)
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())