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
import multiprocessing
from pathlib import Path
import numpy as np
import pandas as pd

# Add the package and phases to the path
sys.path.insert(0, str(Path(__file__).parent / "metagrouper_package"))
sys.path.insert(0, str(Path(__file__).parent / "phases"))

from metagrouper import (
    find_fastq_files,
    setup_logging,
    save_results,
    SourmashProfiler,
    SimilarityAnalyzer,
    Visualizer
)

# Required sourmash imports
import sourmash

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

# Phase 4 imports
try:
    from metagrouper.interactive_visualizer import InteractiveVisualizer
    PHASE4_AVAILABLE = True
except ImportError as e:
    PHASE4_AVAILABLE = False
    logging.warning(f"Phase 4 interactive visualizations not available: {e}")


def get_memory_usage():
    """Get current memory usage in MB."""
    try:
        import psutil
        process = psutil.Process()
        return process.memory_info().rss / 1024 / 1024
    except ImportError:
        return None


def log_memory_usage(stage_name, start_memory=None):
    """Log memory usage at different stages."""
    current_memory = get_memory_usage()
    if current_memory is not None:
        if start_memory is not None:
            diff = current_memory - start_memory
            print(f"💾 {stage_name}: {current_memory:.1f} MB (+{diff:.1f} MB)")
        else:
            print(f"💾 {stage_name}: {current_memory:.1f} MB")
    return current_memory


def get_process_count(args):
    """Get the number of processes to use, defaulting to CPU count if not specified."""
    if args.processes is not None:
        return args.processes
    elif args.sequential:
        return 1
    else:
        # Auto-detect CPU count
        cpu_count = multiprocessing.cpu_count()
        print(f"🔧 Auto-detected {cpu_count} CPU cores, using all cores for processing")
        return cpu_count


def get_dense_distance_matrix(distance_matrix, sparse_similarity_matrix):
    """Convert sparse similarity to dense distance matrix only when needed."""
    if distance_matrix is not None:
        return distance_matrix
    elif sparse_similarity_matrix is not None:
        start_mem = log_memory_usage("Before sparse-to-dense conversion")
        print(f"💾 Converting sparse to dense matrix for downstream analysis...")
        dense_matrix = 1 - sparse_similarity_matrix.toarray()
        log_memory_usage("After sparse-to-dense conversion", start_mem)
        return dense_matrix
    else:
        raise ValueError("No distance matrix available")


def generate_visualizations(args, sample_names, distance_matrix, sparse_similarity_matrix, 
                          pca_result, pca, output_path, metadata_for_viz=None):
    """Generate both static and interactive visualizations efficiently."""
    print(f"\n📈 Generating visualizations...")
    
    try:
        # Generate interactive visualizations if requested
        if (args.interactive or args.interactive_only) and PHASE4_AVAILABLE:
            print(f"🌐 Generating interactive HTML visualizations...")
            interactive_viz = InteractiveVisualizer(sample_names, metadata_for_viz)
            
            if pca_result is not None:
                # Interactive PCA plot
                interactive_viz.create_interactive_pca(
                    pca_result, pca, 
                    output_path / "interactive_pca.html",
                    title=f"{args.html_title} - PCA Analysis"
                )
                
                # Interactive heatmap (for reasonable dataset sizes)
                if len(sample_names) <= 100:
                    dense_distance_matrix = get_dense_distance_matrix(distance_matrix, sparse_similarity_matrix)
                    interactive_viz.create_interactive_heatmap(
                        dense_distance_matrix, 
                        output_path / "interactive_heatmap.html",
                        title=f"{args.html_title} - Distance Heatmap"
                    )
                
                # Unified dashboard
                dense_distance_matrix = get_dense_distance_matrix(distance_matrix, sparse_similarity_matrix)
                interactive_viz.create_dashboard(
                    pca_result, pca, dense_distance_matrix,
                    output_path / "interactive_dashboard.html", 
                    title=f"{args.html_title} - Interactive Dashboard"
                )
                
                print(f"✅ Interactive visualizations saved")
        
        # Generate static visualizations (unless interactive-only is specified)
        if not args.interactive_only:
            visualizer = Visualizer(sample_names)
            
            # Basic distance heatmap
            if len(sample_names) <= 100:  # Only for reasonably sized matrices
                dense_distance_matrix = get_dense_distance_matrix(distance_matrix, sparse_similarity_matrix)
                visualizer.plot_distance_heatmap(dense_distance_matrix, output_path / "distance_heatmap.png")
            
            # PCA plot
            if pca_result is not None:
                visualizer.plot_pca(pca_result, pca, output_path / "pca_plot.png")
            
            print(f"✅ Static visualizations saved")
        
    except Exception as e:
        logging.warning(f"Visualization generation failed: {e}")
        import traceback
        traceback.print_exc()


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
  
  # High sensitivity analysis (robust to PCR bias)
  python metagrouper.py /path/to/fastq/files --scaled 100 --additional-k-sizes 31 51 -o results/
  
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
    
    # Sourmash k-mer profiling arguments (high sensitivity defaults)
    parser.add_argument("--scaled", type=int, default=100,
                       help="Sourmash scaled parameter (1 in N k-mers kept, default: 100 for high sensitivity)")
    parser.add_argument("--num-hashes", type=int, default=0,
                       help="Number of hashes (0 for scaled mode, default: 0)")
    parser.add_argument("--track-abundance", action="store_true", default=False,
                       help="Track k-mer abundances in sourmash signatures (more sensitive to PCR bias)")
    parser.add_argument("--no-track-abundance", action="store_true",
                       help="Explicitly disable k-mer abundance tracking (default behavior)")
    parser.add_argument("--additional-k-sizes", nargs="+", type=int,
                       help="Additional k-mer sizes for multi-scale analysis (e.g., --additional-k-sizes 31 51)")
    parser.add_argument("--save-signatures", action="store_true",
                       help="Save sourmash signatures for future use")
    
    # Metadata Analysis (Phase 2)
    parser.add_argument("-m", "--metadata", help="Metadata file (CSV/TSV) for Phase 2 analysis")
    parser.add_argument("--sample-id-column", default="sample_id",
                       help="Column name for sample IDs in metadata (default: sample_id, common alternatives: Sample_ID, sample, accession, run_id)")
    parser.add_argument("--variables", nargs="+",
                       help="Specific metadata variables to analyze (default: all)")
    parser.add_argument("--auto-filter-variables", action="store_true",
                       help="Automatically filter metadata variables to focus on biologically relevant ones")
    parser.add_argument("--exclude-variables", nargs="+",
                       help="Metadata variables to exclude from analysis")
    parser.add_argument("--permutations", type=int, default=999,
                       help="Number of permutations for PERMANOVA (default: 999)")
    parser.add_argument("--cluster-range", nargs=2, type=int, default=[2, 8],
                       help="Range for number of clusters to test (default: 2 8)")
    
    # Assembly Recommendations (Phase 3)
    parser.add_argument("--assembly-tools", nargs="+",
                       choices=["megahit", "spades", "flye", "all"],
                       default=["megahit", "spades"],
                       help="Assembly tools to generate commands for (default: megahit spades)")
    parser.add_argument("--similarity-threshold", type=float, default=0.45,
                       help="Distance threshold for grouping samples (default: 0.45)")
    parser.add_argument("--min-group-size", type=int, default=2,
                       help="Minimum samples per assembly group (default: 2)")
    parser.add_argument("--max-group-size", type=int, default=20,
                       help="Maximum samples per assembly group (default: 20)")
    
    # Processing arguments
    parser.add_argument("--processes", type=int,
                       help="Number of parallel processes (default: CPU count)")
    parser.add_argument("--sequential", action="store_true",
                       help="Use sequential processing")
    
    # Interactive visualization (Phase 4)
    parser.add_argument("--interactive", action="store_true",
                       help="Generate interactive HTML visualizations")
    parser.add_argument("--interactive-only", action="store_true", 
                       help="Generate only interactive plots, skip static plots")
    parser.add_argument("--comprehensive-report", action="store_true",
                       help="Generate comprehensive interactive HTML report with explanations")
    parser.add_argument("--html-title", default="MetaGrouper Analysis",
                       help="Title for interactive HTML reports")
    
    # Phase selection arguments
    parser.add_argument("--phases", nargs="+", type=int, choices=[1, 2, 3, 4],
                       help="Run specific phases only (e.g., --phases 2 3 to run only Phase 2 and 3)")
    parser.add_argument("--skip-phases", nargs="+", type=int, choices=[1, 2, 3, 4],
                       help="Skip specific phases (e.g., --skip-phases 1 to skip k-mer profiling)")
    parser.add_argument("--load-from", 
                       help="Load Phase 1 results from previous run directory (skip k-mer profiling)")
    
    # Other arguments
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    
    return parser


def run_analysis(args):
    """Run the main analysis workflow."""
    # Setup logging
    setup_logging(args.verbose)
    
    print("🧬 MetaGrouper: K-mer Analysis for Optimal Metagenomic Assembly Grouping")
    print("🔧 VERSION: 2024-07-11 Phase Selection Enhanced (DEBUGGING)")
    print("=" * 80)
    
    # Initial memory usage
    initial_memory = log_memory_usage("Initial memory usage")
    
    # Determine which phases to run based on arguments
    print(f"🔍 DEBUG: args.phases = {getattr(args, 'phases', 'NOT FOUND')}")
    print(f"🔍 DEBUG: args.skip_phases = {getattr(args, 'skip_phases', 'NOT FOUND')}")
    print(f"🔍 DEBUG: args.load_from = {getattr(args, 'load_from', 'NOT FOUND')}")
    print(f"🔍 DEBUG: ALL ARGS: {vars(args)}")
    
    phases_to_run = set([1, 2, 3, 4])  # Default: all phases
    print(f"🔍 DEBUG: Initial phases_to_run = {sorted(phases_to_run)}")
    
    if args.phases:
        # If specific phases requested, run only those
        phases_to_run = set(args.phases)
        print(f"🔍 DEBUG: After --phases, phases_to_run = {sorted(phases_to_run)}")
        
    if args.skip_phases:
        # Remove skipped phases
        phases_to_run -= set(args.skip_phases)
        print(f"🔍 DEBUG: After --skip-phases, phases_to_run = {sorted(phases_to_run)}")
        
    if args.load_from:
        print(f"🔍 DEBUG: --load-from specified: {args.load_from}")
        # If loading from previous run, skip Phase 1 by default (unless explicitly requested)
        if args.phases is None and args.skip_phases is None:
            # Default behavior when using --load-from: skip Phase 1
            phases_to_run.discard(1)
            print(f"🔍 DEBUG: Default --load-from behavior, phases_to_run = {sorted(phases_to_run)}")
        elif 1 in phases_to_run and args.phases and 1 in args.phases:
            # User explicitly requested Phase 1 even with --load-from, warn them
            print("⚠️  Warning: --load-from specified but Phase 1 is requested. Will run Phase 1 anyway.")
        elif 1 in phases_to_run:
            # Phase 1 is in the list but --load-from is specified, remove it
            phases_to_run.discard(1)
            print(f"🔍 DEBUG: Removed Phase 1 due to --load-from, phases_to_run = {sorted(phases_to_run)}")
    else:
        print(f"🔍 DEBUG: No --load-from specified")
        
    # Validate phase selection
    if not phases_to_run:
        logging.error("No phases selected to run")
        return False
        
    print(f"📋 Phases to run: {sorted(phases_to_run)}")
    
    # Create output directory
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Initialize variables for later phases
    profiles = None
    distance_matrix = None
    sample_names = None
    fastq_files = None
    
    # Handle loading from previous run
    if args.load_from:
        print(f"\n📂 Loading Phase 1 results from: {args.load_from}")
        load_path = Path(args.load_from)
        
        try:
            # Load k-mer profiles
            import pickle
            import json
            
            with open(load_path / "kmer_profiles.pkl", "rb") as f:
                profiles = pickle.load(f)
            
            # Load distance matrix
            distance_matrix = np.load(load_path / "distance_matrix.npy")
            
            # Load sample names
            with open(load_path / "sample_names.json", "r") as f:
                sample_names = json.load(f)
                
            print(f"✅ Loaded results for {len(sample_names)} samples")
            
        except FileNotFoundError as e:
            logging.error(f"Could not load saved results: {e}")
            return False
    
    # Validate input for Phase 1
    if 1 in phases_to_run:
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
    
    # Get process count (auto-detect if not specified)
    process_count = get_process_count(args)
    
    # Determine which phases to run (updated to use phase selection)
    run_phase2 = 2 in phases_to_run and args.metadata and PHASE2_AVAILABLE
    run_phase3 = 3 in phases_to_run and PHASE3_AVAILABLE  # Phase 3 can run without metadata
    
    print(f"📋 Analysis Plan:")
    print(f"   Phase 1: K-mer profiling and similarity {'✅' if 1 in phases_to_run else '⏭️  (skipped)'}")
    print(f"   Phase 2: Metadata analysis {'✅' if run_phase2 else '⏭️  (skipped)' if 2 not in phases_to_run else '❌ (no metadata provided)' if not args.metadata else '❌ (dependencies missing)'}")
    print(f"   Phase 3: Assembly recommendations {'✅' if run_phase3 else '⏭️  (skipped)' if 3 not in phases_to_run else '❌ (dependencies missing)'}")
    run_phase4 = 4 in phases_to_run and (args.comprehensive_report or args.interactive)
    print(f"   Phase 4: Interactive visualizations {'✅' if run_phase4 else '⏭️  (skipped)' if 4 not in phases_to_run else '❌ (not requested)'}")
    print(f"   Processing: {process_count} CPU cores")
    if args.load_from:
        print(f"   Loading Phase 1 results from: {args.load_from}")
    print()
    
    # =============================================================================
    # PHASE 1: K-mer Profiling and Similarity Analysis
    # =============================================================================
    if 1 in phases_to_run:
        print("🧬 Phase 1: K-mer Profiling and Similarity Analysis")
        print("-" * 60)
    
    # Handle track_abundance logic (default to False for robustness to PCR bias)
    track_abundance = args.track_abundance and not getattr(args, 'no_track_abundance', False)
    
    # Always use sourmash for k-mer profiling
    print(f"⚡ Using sourmash for fast MinHash k-mer sketching")
    print(f"   K-mer size: {args.kmer_size}")
    print(f"   Scaled: {args.scaled} (higher sensitivity than default 1000)")
    print(f"   Track abundance: {track_abundance} (presence/absence mode more robust to PCR bias)")
    if hasattr(args, 'additional_k_sizes') and args.additional_k_sizes:
        print(f"   Additional k-mer sizes: {args.additional_k_sizes} (multi-scale analysis)")
    elif args.scaled <= 100:
        print(f"   Multi-scale analysis: k=21,31,51 (auto-enabled for high sensitivity)")
    
    profiler = SourmashProfiler(
        k=args.kmer_size,
        scaled=args.scaled,
        num_hashes=args.num_hashes,
        processes=process_count,
        track_abundance=track_abundance,
        additional_k_sizes=getattr(args, 'additional_k_sizes', None)
    )
    
    # Process samples
    print(f"🔬 Processing {len(fastq_files)} samples using {process_count} processes...")
    start_time = time.time()
    
    # Process samples with sourmash
    signatures = profiler.process_samples_parallel(fastq_files)
    
    # Use primary k-mer size for similarity analysis
    similarity_matrix = profiler.compute_similarity_matrix(signatures, use_k_size=args.kmer_size)
    
    profiles, sample_names = profiler.export_to_metagrouper_format(
        signatures, 
        similarity_matrix,
        use_k_size=args.kmer_size
    )
    failed_samples = []
    
    # Save signatures if requested
    if args.save_signatures:
        sig_path = output_path / "signatures.sig"
        profiler.save_signatures(signatures, str(sig_path))
        print(f"💾 Saved signatures to {sig_path}")
    
    processing_time = time.time() - start_time
    success_count = len(profiles)
    
    print(f"✅ Processed {success_count}/{len(fastq_files)} samples in {processing_time:.1f}s")
    if failed_samples:
        print(f"❌ Failed samples: {', '.join(failed_samples)}")
    
    # Memory usage after processing
    log_memory_usage("After k-mer profiling", initial_memory)
    
    if not profiles:
        logging.error("No samples processed successfully")
        return False
    
    # Memory usage report
    print(f"💾 Using sourmash MinHash sketches for efficient memory usage")
    
    # Convert similarity to distance matrix
    print(f"\n🔗 Converting similarity to distance matrix...")
    start_time = time.time()
    
    # Convert similarity to distance matrix
    distance_matrix = 1 - similarity_matrix
    
    similarity_time = time.time() - start_time
    print(f"✅ Similarity matrix converted to distances in {similarity_time:.1f}s")
    log_memory_usage("After similarity matrix conversion")
    
    # Save results
    save_results(profiles, distance_matrix, sample_names, args.output)
    
    # No sparse matrix handling needed with sourmash
    
    # Perform dimensionality reduction for visualization
    pca_result, pca = None, None
    try:
        from sklearn.decomposition import PCA
        
        pca = PCA(n_components=min(2, len(sample_names)-1))
        profile_matrix = 1 - distance_matrix
        pca_result = pca.fit_transform(profile_matrix)
    except Exception as e:
        logging.warning(f"PCA failed: {e}")
        pca_result, pca = None, None
    
    # Load metadata if available for visualizations
    metadata_for_viz = None
    if args.metadata:
        try:
            metadata_for_viz = pd.read_csv(args.metadata, sep=None, engine='python')
        except Exception as e:
            logging.warning(f"Could not load metadata for visualization: {e}")
    
    # Generate visualizations
    generate_visualizations(args, sample_names, distance_matrix, None, 
                          pca_result, pca, output_path, metadata_for_viz)
    
    # =============================================================================
    # PHASE 2: Metadata Analysis (if requested)
    # =============================================================================
    metadata_results_df = pd.DataFrame()
    cluster_results = {}
    
    if run_phase2:
        print(f"\n📊 Phase 2: Metadata Analysis")
        print("-" * 60)
        
        try:
            # Initialize metadata analyzer with distance matrix
            meta_analyzer = MetadataAnalyzer(distance_matrix, sample_names)
            meta_analyzer.load_metadata(args.metadata, args.sample_id_column)
            
            # Analyze variables (PERMANOVA)
            metadata_results_df = meta_analyzer.analyze_variables(
                variables=args.variables, 
                n_permutations=args.permutations,
                auto_filter=args.auto_filter_variables,
                exclude_variables=args.exclude_variables
            )
            
            print(f"✅ Analyzed {len(metadata_results_df)} metadata variables")
            
            # Identify clusters (separate try-catch to not lose PERMANOVA results)
            cluster_results = {}
            try:
                cluster_results = meta_analyzer.identify_clusters(
                    n_clusters_range=tuple(args.cluster_range)
                )
                print(f"✅ Clustering analysis complete")
            except Exception as cluster_error:
                logging.warning(f"Clustering analysis failed: {cluster_error}")
                print(f"⚠️  Clustering analysis failed: {cluster_error}")
                print(f"✅ PERMANOVA results still available")
            
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
            
            # Save variable filtering report if filtering was applied
            if hasattr(meta_analyzer, 'filtering_report_data'):
                report_text = meta_analyzer.generate_filtering_report(
                    meta_analyzer.filtering_report_data['exclusion_reasons'],
                    meta_analyzer.filtering_report_data['included_variables'],
                    output_path / "variable_filtering_report.md"
                )
                if args.auto_filter_variables:
                    print(f"📋 Variable filtering report saved to: variable_filtering_report.md")
            
            print(f"✅ Phase 2 analysis completed successfully")
            
            # Display statistical testing summary
            if not metadata_results_df.empty:
                print(f"\n📈 Statistical Testing Summary (PERMANOVA):")
                print("-" * 45)
                
                total_variables = len(metadata_results_df)
                significant_vars = metadata_results_df[metadata_results_df['p_value'] < 0.05]
                num_significant = len(significant_vars)
                
                # Get top variable info
                top_variable = metadata_results_df.iloc[0] if len(metadata_results_df) > 0 else None
                
                print(f"📊 Tested {total_variables} metadata variables for association with sample composition")
                
                if num_significant > 0:
                    print(f"✅ Found {num_significant} significant variable{'s' if num_significant != 1 else ''} (p < 0.05)")
                    if top_variable is not None:
                        r_squared_pct = top_variable['r_squared'] * 100
                        print(f"🏆 Top variable '{top_variable['variable']}' explains {r_squared_pct:.1f}% of variation (R² = {top_variable['r_squared']:.3f})")
                        
                        # Biological interpretation
                        if top_variable['r_squared'] > 0.3:
                            interpretation = "strong biological association"
                        elif top_variable['r_squared'] > 0.15:
                            interpretation = "moderate biological association"
                        else:
                            interpretation = "weak but detectable association"
                        print(f"💡 This indicates a {interpretation} between this variable and microbial composition")
                else:
                    print(f"⚠️  No variables showed significant association (all p ≥ 0.05)")
                    print(f"💭 Consider similarity-based grouping instead of metadata-based grouping")
                
                print(f"📄 Detailed results saved to: permanova_results.csv")
                print()
            
        except Exception as e:
            logging.error(f"Phase 2 analysis failed: {e}")
            print(f"❌ Phase 2 analysis failed: {e}")
            run_phase2 = False
    
    # =============================================================================
    # PHASE 3: Assembly Strategy Recommendations
    # =============================================================================
    assembly_recommendation = None
    grouping_recommendations = None
    
    if run_phase3:
        print(f"\n🔧 Phase 3: Assembly Strategy Recommendations")
        print("-" * 60)
        
        try:
            # Initialize assembly recommender with distance matrix
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
            
            # Generate and display intelligent grouping recommendations
            if run_phase2 and metadata_results_df is not None:
                print(f"\n📊 Metadata Grouping Recommendations:")
                print("-" * 40)
                
                grouping_recommendations = recommender.generate_metadata_grouping_recommendations(
                    metadata_results_df, 
                    getattr(meta_analyzer, 'metadata', None)
                )
                
                if grouping_recommendations:
                    print(f"💡 Top assembly grouping strategies based on your metadata:")
                    print()
                    
                    for i, rec in enumerate(grouping_recommendations[:5], 1):  # Show top 5
                        confidence_emoji = "🟢" if rec['confidence'] > 0.7 else "🟡" if rec['confidence'] > 0.4 else "🔴"
                        strategy_emoji = "👥" if rec['strategy'] == "grouped_coassembly" else "🔄"
                        
                        print(f"{i}. {confidence_emoji} {strategy_emoji} {rec['recommendation_text']}")
                        print(f"   📈 {rec['rationale']}")
                        if rec['benefits']:
                            print(f"   ✅ Benefits: {', '.join(rec['benefits'])}")
                        if rec['challenges']:
                            print(f"   ⚠️  Challenges: {', '.join(rec['challenges'])}")
                        print()
                    
                    # Save detailed recommendations to file
                    grouping_file = output_path / "metadata_grouping_recommendations.json"
                    import json
                    with open(grouping_file, 'w') as f:
                        json.dump(grouping_recommendations, f, indent=2)
                    print(f"💾 Detailed recommendations saved to: {grouping_file}")
                    print()
                else:
                    print(f"   ℹ️  No strong metadata-based grouping strategies found")
                    print(f"   📝 Consider individual assembly or similarity-based grouping")
                    print()
            
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
    
    # Memory efficiency already reported above
    
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
    
    # =============================================================================
    # PHASE 4 (Enhanced): Comprehensive Interactive Report
    # =============================================================================
    
    if (args.comprehensive_report or args.interactive) and 4 in phases_to_run:
        try:
            print(f"\n🌟 Generating comprehensive interactive report...")
            
            from interactive_report_generator import create_interactive_report
            
            # Collect all data for the report
            kmer_data_dict = None
            if 'profiles' in locals() and profiles:
                print(f"🔬 Found k-mer profiles for report: {len(profiles)} samples")
                print(f"📋 Profile sample names: {list(profiles.keys())[:3]}...")
                kmer_data_dict = {'profiles': profiles}
            else:
                print(f"❌ No k-mer profiles found for interactive report")
            
            # Generate the comprehensive report
            metadata_for_report = None
            if args.metadata:
                print(f"🔍 Loading metadata from: {args.metadata}")
                try:
                    # Check if file exists
                    if not Path(args.metadata).exists():
                        print(f"❌ Metadata file not found: {args.metadata}")
                        logging.error(f"Metadata file not found: {args.metadata}")
                    else:
                        # Load metadata directly for visualization (even if Phase 2 didn't run)
                        metadata_for_report = pd.read_csv(args.metadata, sep=None, engine='python')
                        print(f"✅ Loaded metadata: {len(metadata_for_report)} rows, {len(metadata_for_report.columns)} columns")
                        print(f"📋 Metadata columns: {list(metadata_for_report.columns)}")
                        print(f"📝 First few sample IDs: {metadata_for_report.iloc[:3, 0].tolist() if len(metadata_for_report) > 0 else 'None'}")
                        logging.info(f"Loaded metadata for report: {len(metadata_for_report)} rows, {len(metadata_for_report.columns)} columns")
                except Exception as e:
                    print(f"❌ Error loading metadata: {e}")
                    logging.warning(f"Could not load metadata for report: {e}")
            else:
                print("⚠️  No metadata file specified - interactive coloring will not be available")
                    
            # If Phase 2 ran, use its processed metadata instead
            if run_phase2 and 'meta_analyzer' in locals():
                metadata_for_report = meta_analyzer.metadata
                print(f"✅ Using Phase 2 processed metadata for report")
                logging.info(f"Using Phase 2 processed metadata for report")
            
            report_path = create_interactive_report(
                distance_matrix=distance_matrix,
                sample_names=sample_names,
                output_dir=str(output_path),
                metadata=metadata_for_report,
                permanova_results=metadata_results_df if run_phase2 and 'metadata_results_df' in locals() and not metadata_results_df.empty else None,
                assembly_recommendation=assembly_recommendation if run_phase3 else None,
                kmer_data=kmer_data_dict,
                grouping_recommendations=grouping_recommendations if run_phase3 and grouping_recommendations else None,
                title=args.html_title,
                sample_id_column=args.sample_id_column,
                analyzed_variables=args.variables
            )
            
            print(f"✅ Comprehensive interactive report generated!")
            print(f"   📄 {report_path}")
            print(f"   🎯 Open in browser for interactive exploration")
            
        except Exception as e:
            logging.warning(f"Comprehensive report generation failed: {e}")
            print(f"⚠️  Could not generate comprehensive report: {e}")
    
    print(f"\n📁 Results saved to: {args.output}")
    print(f"📚 Check the output directory for detailed results and visualizations")
    
    return True


def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    
    success = run_analysis(args)
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())