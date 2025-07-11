#!/usr/bin/env python3
"""
Enhanced MetaGrouper with phase selection support.

This version allows running specific phases independently.
"""

import sys
import argparse
import logging
import time
import multiprocessing
from pathlib import Path
import numpy as np
import pandas as pd
import pickle
import json

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

# Import functions from original metagrouper
sys.path.insert(0, str(Path(__file__).parent))
import metagrouper as mg_main

# Get functions we need
try:
    get_memory_usage = mg_main.get_memory_usage
    log_memory_usage = mg_main.log_memory_usage  
    get_process_count = mg_main.get_process_count
    generate_visualizations = mg_main.generate_visualizations
    create_parser = mg_main.create_parser
    PHASE2_AVAILABLE = mg_main.PHASE2_AVAILABLE
    PHASE3_AVAILABLE = mg_main.PHASE3_AVAILABLE  
    PHASE4_AVAILABLE = mg_main.PHASE4_AVAILABLE
except AttributeError:
    # Fallback implementations
    def get_memory_usage():
        try:
            import psutil
            return psutil.Process().memory_info().rss / 1024 / 1024
        except:
            return None
    
    def log_memory_usage(stage_name, start_memory=None):
        memory = get_memory_usage()
        if memory:
            print(f"💾 {stage_name}: {memory:.1f} MB")
        return memory
    
    def get_process_count(args):
        if hasattr(args, 'processes') and args.processes:
            return args.processes
        elif hasattr(args, 'sequential') and args.sequential:
            return 1
        else:
            return multiprocessing.cpu_count()
    
    PHASE2_AVAILABLE = True
    PHASE3_AVAILABLE = True
    PHASE4_AVAILABLE = True

# Phase 2 imports
if PHASE2_AVAILABLE:
    from metadata_analyzer import (
        MetadataAnalyzer,
        MetadataVisualizer,
        generate_summary_report,
    )

# Phase 3 imports
if PHASE3_AVAILABLE:
    from assembly_recommender import (
        AssemblyRecommender,
        save_recommendations,
        visualize_assembly_strategy,
    )

# Phase 4 imports
try:
    from metagrouper.interactive_visualizer import InteractiveVisualizer
    from phases.interactive_report_generator import create_interactive_report
    PHASE4_AVAILABLE = True
except ImportError:
    PHASE4_AVAILABLE = False


def create_enhanced_parser():
    """Create enhanced command-line argument parser with phase selection."""
    try:
        parser = create_parser()
    except:
        # Create a simple parser if we can't import the original
        parser = argparse.ArgumentParser(description="MetaGrouper with phase selection")
        parser.add_argument("input_dir", help="Directory containing FASTQ files")
        parser.add_argument("-o", "--output", default="metagrouper_output", help="Output directory")
        parser.add_argument("-m", "--metadata", help="Metadata file")
        parser.add_argument("--sample-id-column", default="database_ID", help="Sample ID column")
        parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    
    # Update the description
    parser.description = """
MetaGrouper: K-mer profiling and metadata-driven assembly grouping with phase selection support.

This enhanced version allows running specific phases independently:
- Phase 1: K-mer profiling and similarity analysis
- Phase 2: Metadata analysis and PERMANOVA 
- Phase 3: Assembly strategy recommendations
- Phase 4: Interactive visualizations

You can run specific phases using --phases or skip phases using --skip-phases.
Load previous results using --load-from to skip Phase 1.
"""
    
    # Update examples in epilog
    parser.epilog = """
Examples:
  # Run all phases (default)
  python metagrouper_phases.py /path/to/fastq/files -m metadata.csv -o results/
  
  # Run only Phase 1 (k-mer profiling)
  python metagrouper_phases.py /path/to/fastq/files -o results/ --phases 1
  
  # Run only Phase 2 and 3 using previous results
  python metagrouper_phases.py dummy_path -m metadata.csv -o new_results/ --load-from results/ --phases 2 3
  
  # Skip Phase 1 (same as above)
  python metagrouper_phases.py dummy_path -m metadata.csv -o new_results/ --load-from results/ --skip-phases 1
  
  # Run Phase 2 with different metadata on existing profiles
  python metagrouper_phases.py dummy_path -m new_metadata.csv -o reanalysis/ --load-from results/ --phases 2
"""
    
    # Add phase selection arguments (these work with any parser)
    parser.add_argument("--phases", nargs="+", type=int, choices=[1, 2, 3, 4],
                       help="Run specific phases only (e.g., --phases 2 3)")
    parser.add_argument("--skip-phases", nargs="+", type=int, choices=[1, 2, 3, 4],
                       help="Skip specific phases (e.g., --skip-phases 1)")
    parser.add_argument("--load-from", 
                       help="Load Phase 1 results from previous run directory")
    
    return parser


def run_phase1(args, output_path):
    """Run Phase 1: K-mer profiling and similarity analysis."""
    print("🧬 Phase 1: K-mer Profiling and Similarity Analysis")
    print("-" * 60)
    
    # Find FASTQ files
    logging.info("Finding FASTQ files...")
    fastq_files = find_fastq_files(args.input_dir)
    if not fastq_files:
        logging.error("No FASTQ files found")
        return None, None, None
    
    logging.info(f"Found {len(fastq_files)} samples")
    
    # Get process count
    process_count = get_process_count(args)
    
    # Handle track_abundance logic
    track_abundance = args.track_abundance and not getattr(args, 'no_track_abundance', False)
    
    # Setup profiler
    print(f"⚡ Using sourmash for fast MinHash k-mer sketching")
    print(f"   K-mer size: {args.kmer_size}")
    print(f"   Scaled: {args.scaled}")
    print(f"   Track abundance: {track_abundance}")
    
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
    
    signatures = profiler.process_samples_parallel(fastq_files)
    similarity_matrix = profiler.compute_similarity_matrix(signatures, use_k_size=args.kmer_size)
    
    profiles, sample_names = profiler.export_to_metagrouper_format(
        signatures, 
        similarity_matrix,
        use_k_size=args.kmer_size
    )
    
    # Save signatures if requested
    if args.save_signatures:
        sig_path = output_path / "signatures.sig"
        profiler.save_signatures(signatures, str(sig_path))
        print(f"💾 Saved signatures to {sig_path}")
    
    processing_time = time.time() - start_time
    print(f"✅ Processed {len(profiles)} samples in {processing_time:.1f}s")
    
    # Convert similarity to distance
    distance_matrix = 1 - similarity_matrix
    np.fill_diagonal(distance_matrix, 0)
    
    # Save results
    save_results(profiles, distance_matrix, sample_names, str(output_path))
    
    # Also save k-mer data for Phase 4
    kmer_data_dict = {}
    for sample_name in sample_names:
        kmer_data_dict[sample_name] = profiles[sample_name]
    
    with open(output_path / "kmer_data.pkl", "wb") as f:
        pickle.dump(kmer_data_dict, f)
    
    return profiles, distance_matrix, sample_names


def run_phase2(args, distance_matrix, sample_names, output_path):
    """Run Phase 2: Metadata analysis."""
    if not PHASE2_AVAILABLE:
        print("❌ Phase 2 dependencies not available")
        return None, {}
        
    print(f"\n📊 Phase 2: Metadata Analysis")
    print("-" * 60)
    
    if not args.metadata:
        print("❌ No metadata file provided")
        return None, {}
    
    try:
        # Initialize metadata analyzer
        meta_analyzer = MetadataAnalyzer(distance_matrix, sample_names)
        meta_analyzer.load_metadata(args.metadata, args.sample_id_column)
        
        # Analyze variables (disable auto-filter to avoid the bug)
        metadata_results_df = meta_analyzer.analyze_variables(
            variables=args.variables,
            n_permutations=args.permutations,
            auto_filter=False  # Disable to avoid the filtering bug
        )
        
        if not metadata_results_df.empty:
            print(f"✅ Analyzed {len(metadata_results_df)} metadata variables")
            
            # Save results
            metadata_results_df.to_csv(output_path / "permanova_results.csv", index=False)
            
            # Generate visualizations
            meta_visualizer = MetadataVisualizer(sample_names, meta_analyzer.metadata)
            
            # Variable importance plot
            meta_visualizer.plot_variable_importance(
                metadata_results_df, output_path / "variable_importance.png"
            )
            
            # Generate summary report
            generate_summary_report(
                metadata_results_df, {}, output_path / "analysis_report.md"
            )
            
            print(f"✅ Phase 2 analysis completed successfully")
            
            return metadata_results_df, meta_analyzer
        else:
            print(f"⚠️  No variables were analyzed")
            return pd.DataFrame(), meta_analyzer
            
    except Exception as e:
        logging.error(f"Phase 2 analysis failed: {e}")
        print(f"❌ Phase 2 analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None, {}


def run_phase3(args, distance_matrix, sample_names, metadata_results_df, meta_analyzer, output_path):
    """Run Phase 3: Assembly recommendations."""
    if not PHASE3_AVAILABLE:
        print("❌ Phase 3 dependencies not available")
        return None, None
        
    print(f"\n🔧 Phase 3: Assembly Strategy Recommendations")
    print("-" * 60)
    
    try:
        # Initialize assembly recommender
        recommender = AssemblyRecommender(distance_matrix, sample_names)
        
        # Configure thresholds
        recommender.strategy_engine.similarity_threshold_medium = args.similarity_threshold
        
        # Generate recommendations
        tools = ["megahit", "spades", "flye"] if "all" in args.assembly_tools else args.assembly_tools
        
        metadata = getattr(meta_analyzer, 'metadata', None) if meta_analyzer else None
        
        assembly_recommendation = recommender.generate_recommendations(
            metadata_results=metadata_results_df,
            metadata=metadata
        )
        
        # Filter assembly commands
        filtered_commands = {
            tool: commands
            for tool, commands in assembly_recommendation.assembly_commands.items()
            if tool in tools
        }
        assembly_recommendation.assembly_commands = filtered_commands
        
        # Save recommendations
        save_recommendations(assembly_recommendation, output_path)
        
        # Visualize strategy
        visualize_assembly_strategy(
            assembly_recommendation, 
            distance_matrix,
            sample_names,
            output_path / "assembly_strategy.png"
        )
        
        print(f"✅ Phase 3 analysis completed successfully")
        
        return assembly_recommendation, None
        
    except Exception as e:
        logging.error(f"Phase 3 analysis failed: {e}")
        print(f"❌ Phase 3 analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None


def run_phase4(args, distance_matrix, sample_names, metadata_results_df, assembly_recommendation, output_path):
    """Run Phase 4: Interactive visualizations."""
    if not PHASE4_AVAILABLE:
        print("❌ Phase 4 dependencies not available")
        return
        
    if not args.comprehensive_report:
        return
        
    print(f"\n🌟 Phase 4: Generating comprehensive interactive report...")
    
    try:
        # Load k-mer data if available
        kmer_data = None
        if (output_path / "kmer_data.pkl").exists():
            with open(output_path / "kmer_data.pkl", "rb") as f:
                kmer_data = pickle.load(f)
        
        # Load metadata for visualization
        metadata_for_report = None
        if args.metadata:
            metadata_for_report = pd.read_csv(args.metadata, sep=None, engine='python')
        
        report_path = create_interactive_report(
            distance_matrix=distance_matrix,
            sample_names=sample_names,
            output_dir=str(output_path),
            metadata=metadata_for_report,
            permanova_results=metadata_results_df if metadata_results_df is not None else None,
            assembly_recommendation=assembly_recommendation,
            kmer_data=kmer_data,
            title=args.html_title,
            sample_id_column=args.sample_id_column
        )
        
        print(f"✅ Comprehensive interactive report generated!")
        print(f"   📄 {report_path}")
        
    except Exception as e:
        logging.warning(f"Phase 4 report generation failed: {e}")
        print(f"⚠️  Could not generate comprehensive report: {e}")


def run_analysis(args):
    """Run the main analysis workflow with phase selection."""
    # Setup logging
    setup_logging(args.verbose)
    
    print("🧬 MetaGrouper: K-mer Analysis with Phase Selection Support")
    print("=" * 80)
    
    # Initial memory usage
    initial_memory = log_memory_usage("Initial memory usage")
    
    # Determine which phases to run
    phases_to_run = set([1, 2, 3, 4])  # Default: all phases
    
    if args.phases:
        phases_to_run = set(args.phases)
        
    if args.skip_phases:
        phases_to_run -= set(args.skip_phases)
        
    if args.load_from:
        phases_to_run.discard(1)
        
    if not phases_to_run:
        logging.error("No phases selected to run")
        return False
        
    print(f"📋 Phases to run: {sorted(phases_to_run)}")
    
    # Create output directory
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Initialize variables
    profiles = None
    distance_matrix = None
    sample_names = None
    metadata_results_df = None
    meta_analyzer = None
    assembly_recommendation = None
    
    # Load from previous run if specified
    if args.load_from and 1 not in phases_to_run:
        print(f"\n📂 Loading Phase 1 results from: {args.load_from}")
        load_path = Path(args.load_from)
        
        try:
            with open(load_path / "kmer_profiles.pkl", "rb") as f:
                profiles = pickle.load(f)
            
            distance_matrix = np.load(load_path / "distance_matrix.npy")
            
            with open(load_path / "sample_names.json", "r") as f:
                sample_names = json.load(f)
                
            print(f"✅ Loaded results for {len(sample_names)} samples")
            
        except FileNotFoundError as e:
            logging.error(f"Could not load saved results: {e}")
            print(f"❌ Error: Required files not found in {args.load_from}")
            print(f"   Expected: kmer_profiles.pkl, distance_matrix.npy, sample_names.json")
            return False
    
    # Run Phase 1 if requested
    if 1 in phases_to_run:
        profiles, distance_matrix, sample_names = run_phase1(args, output_path)
        if profiles is None:
            return False
    
    # Run Phase 2 if requested
    if 2 in phases_to_run:
        if distance_matrix is None or sample_names is None:
            print("❌ Phase 2 requires Phase 1 results (distance matrix and sample names)")
            return False
        metadata_results_df, meta_analyzer = run_phase2(args, distance_matrix, sample_names, output_path)
    
    # Run Phase 3 if requested
    if 3 in phases_to_run:
        if distance_matrix is None or sample_names is None:
            print("❌ Phase 3 requires Phase 1 results (distance matrix and sample names)")
            return False
        assembly_recommendation, _ = run_phase3(args, distance_matrix, sample_names, 
                                              metadata_results_df, meta_analyzer, output_path)
    
    # Run Phase 4 if requested
    if 4 in phases_to_run:
        if distance_matrix is None or sample_names is None:
            print("❌ Phase 4 requires Phase 1 results (distance matrix and sample names)")
            return False
        run_phase4(args, distance_matrix, sample_names, metadata_results_df, 
                  assembly_recommendation, output_path)
    
    print(f"\n📁 Results saved to: {args.output}")
    print(f"📚 Analysis complete!")
    
    return True


def main():
    """Main entry point."""
    parser = create_enhanced_parser()
    args = parser.parse_args()
    
    success = run_analysis(args)
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())