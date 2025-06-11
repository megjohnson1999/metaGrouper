#!/usr/bin/env python3
"""
MetaGrouper Command Line Interface

Enhanced CLI with better user experience, error handling, and progress reporting.
"""

import argparse
import sys
import os
import logging
from pathlib import Path
from typing import List, Optional
import time
from datetime import datetime
import json

def create_enhanced_parser() -> argparse.ArgumentParser:
    """Create enhanced argument parser with better organization and help."""
    
    parser = argparse.ArgumentParser(
        prog='metagrouper',
        description='MetaGrouper: K-mer-based analysis for optimal metagenomic assembly grouping',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic k-mer analysis
  metagrouper samples/ -o results/

  # With metadata analysis
  metagrouper samples/ -m metadata.csv -o results/

  # Full analysis with assembly recommendations
  metagrouper samples/ -m metadata.csv --assembly-tools all -o results/

  # Quick test with reduced parameters
  metagrouper samples/ --kmer-size 15 --max-reads 1000 --permutations 99

For more information, see the tutorial: https://github.com/user/metagrouper/blob/main/TUTORIAL.md
        """
    )
    
    # Required arguments
    parser.add_argument(
        'input_dir',
        help='Directory containing FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz)'
    )
    
    # Output and general options
    general_group = parser.add_argument_group('General Options')
    general_group.add_argument(
        '-o', '--output',
        default='metagrouper_output',
        help='Output directory (default: metagrouper_output)'
    )
    general_group.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    general_group.add_argument(
        '--version',
        action='version',
        version='MetaGrouper 1.0.0'
    )
    
    # Phase 1: K-mer analysis options
    kmer_group = parser.add_argument_group('K-mer Analysis (Phase 1)')
    kmer_group.add_argument(
        '-k', '--kmer-size',
        type=int,
        default=21,
        choices=range(11, 32, 2),  # Odd numbers only
        help='K-mer size (must be odd, 11-31, default: 21)'
    )
    kmer_group.add_argument(
        '--max-reads',
        type=int,
        help='Maximum reads per sample (for testing/large datasets)'
    )
    kmer_group.add_argument(
        '--distance-metric',
        choices=['braycurtis', 'jaccard', 'cosine', 'euclidean'],
        default='braycurtis',
        help='Distance metric for sample comparison (default: braycurtis)'
    )
    
    # Phase 2: Metadata analysis options
    metadata_group = parser.add_argument_group('Metadata Analysis (Phase 2)')
    metadata_group.add_argument(
        '-m', '--metadata',
        help='Metadata file (CSV/TSV) for statistical analysis'
    )
    metadata_group.add_argument(
        '--sample-id-column',
        default='sample_id',
        help='Column name for sample IDs in metadata (default: sample_id)'
    )
    metadata_group.add_argument(
        '--variables',
        nargs='+',
        help='Specific metadata variables to analyze (default: all columns)'
    )
    metadata_group.add_argument(
        '--permutations',
        type=int,
        default=999,
        help='Number of permutations for PERMANOVA (default: 999, min: 99)'
    )
    metadata_group.add_argument(
        '--cluster-range',
        nargs=2,
        type=int,
        default=[2, 8],
        help='Range for number of clusters to test (default: 2 8)'
    )
    
    # Phase 3: Assembly recommendation options
    assembly_group = parser.add_argument_group('Assembly Recommendations (Phase 3)')
    assembly_group.add_argument(
        '--assembly-tools',
        nargs='+',
        choices=['megahit', 'spades', 'flye', 'all'],
        default=['megahit', 'spades'],
        help='Assembly tools to generate commands for (default: megahit spades)'
    )
    assembly_group.add_argument(
        '--similarity-threshold',
        type=float,
        default=0.30,
        help='Distance threshold for grouping samples (0.0-1.0, default: 0.30)'
    )
    assembly_group.add_argument(
        '--min-group-size',
        type=int,
        default=2,
        help='Minimum samples per assembly group (default: 2)'
    )
    assembly_group.add_argument(
        '--max-group-size',
        type=int,
        default=10,
        help='Maximum samples per assembly group (default: 10)'
    )
    
    # Performance options
    perf_group = parser.add_argument_group('Performance Options')
    perf_group.add_argument(
        '--quick',
        action='store_true',
        help='Quick analysis mode (k=15, max-reads=1000, permutations=99)'
    )
    perf_group.add_argument(
        '--high-quality',
        action='store_true',
        help='High-quality mode (k=25, permutations=9999)'
    )
    
    return parser

def validate_arguments(args) -> List[str]:
    """Validate command line arguments and return list of errors."""
    errors = []
    
    # Check input directory
    if not os.path.exists(args.input_dir):
        errors.append(f"Input directory does not exist: {args.input_dir}")
    elif not os.path.isdir(args.input_dir):
        errors.append(f"Input path is not a directory: {args.input_dir}")
    
    # Check metadata file if provided
    if args.metadata:
        if not os.path.exists(args.metadata):
            errors.append(f"Metadata file does not exist: {args.metadata}")
        elif not args.metadata.lower().endswith(('.csv', '.tsv', '.txt')):
            errors.append(f"Metadata file must be CSV or TSV: {args.metadata}")
    
    # Validate k-mer size
    if args.kmer_size % 2 == 0:
        errors.append(f"K-mer size must be odd: {args.kmer_size}")
    
    # Validate similarity threshold
    if not 0.0 <= args.similarity_threshold <= 1.0:
        errors.append(f"Similarity threshold must be between 0.0 and 1.0: {args.similarity_threshold}")
    
    # Validate group sizes
    if args.min_group_size < 2:
        errors.append(f"Minimum group size must be at least 2: {args.min_group_size}")
    
    if args.max_group_size < args.min_group_size:
        errors.append(f"Maximum group size must be >= minimum group size")
    
    # Validate permutations
    if args.permutations < 99:
        errors.append(f"Number of permutations must be at least 99: {args.permutations}")
    
    # Validate cluster range
    if args.cluster_range[0] >= args.cluster_range[1]:
        errors.append(f"Cluster range invalid: start must be < end")
    
    if args.cluster_range[0] < 2:
        errors.append(f"Minimum cluster number must be at least 2")
    
    return errors

def apply_preset_modes(args):
    """Apply quick or high-quality mode presets."""
    if args.quick:
        print("🚀 Quick mode enabled: optimized for speed")
        args.kmer_size = 15
        args.max_reads = 1000
        args.permutations = 99
        args.similarity_threshold = 0.35  # More permissive for faster grouping
    
    if args.high_quality:
        print("🎯 High-quality mode enabled: optimized for accuracy")
        args.kmer_size = 25
        args.permutations = 9999
        args.similarity_threshold = 0.20  # More stringent for better groups

def setup_enhanced_logging(verbose: bool, output_dir: str):
    """Set up enhanced logging with progress indicators."""
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(exist_ok=True)
    
    # Set up logging
    log_level = logging.DEBUG if verbose else logging.INFO
    
    # Create formatters
    detailed_formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s'
    )
    simple_formatter = logging.Formatter('%(message)s')
    
    # File handler (detailed)
    file_handler = logging.FileHandler(Path(output_dir) / 'metagrouper.log')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(detailed_formatter)
    
    # Console handler (simple for non-verbose, detailed for verbose)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    
    if verbose:
        console_handler.setFormatter(detailed_formatter)
    else:
        console_handler.setFormatter(simple_formatter)
    
    # Configure root logger
    logging.basicConfig(
        level=logging.DEBUG,
        handlers=[file_handler, console_handler],
        force=True
    )

def print_analysis_plan(args):
    """Print a summary of the planned analysis."""
    print("\n" + "="*60)
    print("📋 METAGROUPER ANALYSIS PLAN")
    print("="*60)
    
    print(f"📁 Input directory: {args.input_dir}")
    print(f"📁 Output directory: {args.output}")
    
    print(f"\n🧬 Phase 1 - K-mer Analysis:")
    print(f"   • K-mer size: {args.kmer_size}")
    print(f"   • Distance metric: {args.distance_metric}")
    if args.max_reads:
        print(f"   • Max reads per sample: {args.max_reads:,}")
    
    if args.metadata:
        print(f"\n📊 Phase 2 - Metadata Analysis:")
        print(f"   • Metadata file: {args.metadata}")
        print(f"   • Sample ID column: {args.sample_id_column}")
        print(f"   • PERMANOVA permutations: {args.permutations:,}")
        if args.variables:
            print(f"   • Variables: {', '.join(args.variables)}")
        else:
            print(f"   • Variables: All columns")
    
    print(f"\n🛠️ Phase 3 - Assembly Recommendations:")
    if 'all' in args.assembly_tools:
        tools = ['MEGAHIT', 'SPAdes', 'Flye']
    else:
        tools = [tool.upper() for tool in args.assembly_tools]
    print(f"   • Assembly tools: {', '.join(tools)}")
    print(f"   • Similarity threshold: {args.similarity_threshold}")
    print(f"   • Group size range: {args.min_group_size}-{args.max_group_size}")
    
    print("\n" + "="*60)

def check_dependencies():
    """Check if required dependencies are available."""
    missing = []
    
    try:
        import numpy
    except ImportError:
        missing.append('numpy')
    
    try:
        import pandas
    except ImportError:
        missing.append('pandas')
    
    try:
        import sklearn
    except ImportError:
        missing.append('scikit-learn')
    
    try:
        import matplotlib
    except ImportError:
        missing.append('matplotlib')
    
    try:
        import seaborn
    except ImportError:
        missing.append('seaborn')
    
    if missing:
        print(f"❌ Missing required dependencies: {', '.join(missing)}")
        print("Install with: pip install " + " ".join(missing))
        return False
    
    return True

def find_fastq_files(input_dir: str) -> List[str]:
    """Find FASTQ files in input directory."""
    fastq_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    files = []
    
    input_path = Path(input_dir)
    for ext in fastq_extensions:
        files.extend(list(input_path.glob(f'*{ext}')))
        files.extend(list(input_path.glob(f'**/*{ext}')))
    
    return sorted([str(f) for f in files])

def print_progress_summary(phase: str, status: str, details: str = ""):
    """Print formatted progress updates."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    
    if status == "start":
        print(f"\n🔄 [{timestamp}] Starting {phase}...")
        if details:
            print(f"   {details}")
    elif status == "complete":
        print(f"✅ [{timestamp}] {phase} completed successfully!")
        if details:
            print(f"   {details}")
    elif status == "error":
        print(f"❌ [{timestamp}] {phase} failed!")
        if details:
            print(f"   Error: {details}")

def save_run_summary(args, start_time: float, output_dir: str, success: bool):
    """Save a summary of the analysis run."""
    
    end_time = time.time()
    duration = end_time - start_time
    
    summary = {
        'timestamp': datetime.now().isoformat(),
        'command': ' '.join(sys.argv),
        'parameters': {
            'input_dir': args.input_dir,
            'output_dir': args.output,
            'kmer_size': args.kmer_size,
            'distance_metric': args.distance_metric,
            'metadata_file': args.metadata,
            'assembly_tools': args.assembly_tools,
            'similarity_threshold': args.similarity_threshold,
            'permutations': args.permutations if args.metadata else None
        },
        'runtime_seconds': duration,
        'runtime_formatted': f"{duration:.1f}s" if duration < 60 else f"{duration/60:.1f}m",
        'success': success,
        'version': '1.0.0'
    }
    
    with open(Path(output_dir) / 'run_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

def print_final_summary(args, start_time: float, fastq_files: List[str], success: bool):
    """Print final analysis summary."""
    
    duration = time.time() - start_time
    
    print("\n" + "="*60)
    if success:
        print("🎉 METAGROUPER ANALYSIS COMPLETE!")
    else:
        print("❌ METAGROUPER ANALYSIS FAILED!")
    print("="*60)
    
    print(f"📊 Processed {len(fastq_files)} FASTQ files")
    print(f"⏱️  Total runtime: {duration:.1f}s" if duration < 60 else f"{duration/60:.1f}m")
    print(f"📁 Results saved to: {args.output}")
    
    if success:
        print(f"\n📋 Generated outputs:")
        
        # Always generated (Phase 1)
        print(f"   Phase 1 (K-mer Analysis):")
        print(f"   • distance_heatmap.png - Sample similarity visualization")
        print(f"   • pca_plot.png - Principal component analysis")
        print(f"   • distance_matrix.csv - Quantitative similarity data")
        
        # Phase 2 outputs
        if args.metadata:
            print(f"   Phase 2 (Metadata Analysis):")
            print(f"   • analysis_report.md - Statistical analysis summary")
            print(f"   • variable_importance.png - PERMANOVA results")
            print(f"   • permanova_results.csv - Statistical test results")
        
        # Phase 3 outputs
        print(f"   Phase 3 (Assembly Recommendations):")
        print(f"   • assembly_recommendations/ - Detailed recommendations")
        print(f"   • assembly_strategy.md - Human-readable strategy")
        print(f"   • run_*_assemblies.sh - Executable assembly scripts")
        
        print(f"\n🚀 Next steps:")
        print(f"   1. Review assembly_strategy.md for recommendations")
        print(f"   2. Run generated assembly scripts")
        print(f"   3. Validate results with your biological knowledge")
    
    print("\n" + "="*60)

def enhanced_main():
    """Enhanced main function with better user experience."""
    
    # Parse arguments
    parser = create_enhanced_parser()
    args = parser.parse_args()
    
    start_time = time.time()
    
    # Apply preset modes
    apply_preset_modes(args)
    
    # Validate arguments
    errors = validate_arguments(args)
    if errors:
        print("❌ Argument validation failed:")
        for error in errors:
            print(f"   • {error}")
        sys.exit(1)
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Set up logging
    setup_enhanced_logging(args.verbose, args.output)
    
    # Find FASTQ files
    fastq_files = find_fastq_files(args.input_dir)
    if not fastq_files:
        print(f"❌ No FASTQ files found in {args.input_dir}")
        print("   Supported formats: .fastq, .fq, .fastq.gz, .fq.gz")
        sys.exit(1)
    
    # Print analysis plan
    print_analysis_plan(args)
    
    print(f"\n📁 Found {len(fastq_files)} FASTQ files")
    if args.verbose:
        for f in fastq_files[:5]:  # Show first 5
            print(f"   • {f}")
        if len(fastq_files) > 5:
            print(f"   • ... and {len(fastq_files) - 5} more")
    
    # Confirmation prompt for large analyses
    if len(fastq_files) > 20 or (args.permutations > 999 and args.metadata):
        print(f"\n⚠️  Large analysis detected:")
        print(f"   • {len(fastq_files)} samples")
        if args.metadata and args.permutations > 999:
            print(f"   • {args.permutations:,} permutations")
        print(f"   • Estimated runtime: {estimate_runtime(args, len(fastq_files))}")
        
        response = input("\nContinue? [y/N]: ").strip().lower()
        if response not in ['y', 'yes']:
            print("Analysis cancelled by user.")
            sys.exit(0)
    
    # Import and run main analysis
    try:
        # Import here to avoid import errors if dependencies missing
        from metagrouper import main as run_metagrouper
        
        print_progress_summary("MetaGrouper Analysis", "start")
        
        # Run the analysis
        success = run_metagrouper_with_progress(args)
        
        if success:
            print_progress_summary("MetaGrouper Analysis", "complete")
        else:
            print_progress_summary("MetaGrouper Analysis", "error")
        
    except Exception as e:
        print_progress_summary("MetaGrouper Analysis", "error", str(e))
        logging.error(f"Analysis failed: {e}", exc_info=True)
        success = False
    
    # Save run summary
    save_run_summary(args, start_time, args.output, success)
    
    # Print final summary
    print_final_summary(args, start_time, fastq_files, success)
    
    sys.exit(0 if success else 1)

def estimate_runtime(args, n_samples: int) -> str:
    """Estimate analysis runtime based on parameters."""
    
    # Base time per sample (seconds)
    base_time = 10 if args.kmer_size <= 17 else 20 if args.kmer_size <= 21 else 30
    
    if args.max_reads:
        base_time *= min(1.0, args.max_reads / 10000)  # Scale with read count
    
    # Add metadata analysis time
    if args.metadata:
        meta_time = args.permutations * n_samples * 0.01  # Rough estimate
        base_time += meta_time
    
    total_time = base_time * n_samples
    
    if total_time < 60:
        return f"{total_time:.0f} seconds"
    elif total_time < 3600:
        return f"{total_time/60:.1f} minutes"
    else:
        return f"{total_time/3600:.1f} hours"

def run_metagrouper_with_progress(args) -> bool:
    """Run MetaGrouper with progress reporting."""
    
    try:
        # This would be replaced with actual MetaGrouper main function
        # For now, just import and run the existing main
        import sys
        sys.argv = ['metagrouper.py', args.input_dir, '--output', args.output]
        
        if args.metadata:
            sys.argv.extend(['--metadata', args.metadata])
        if args.kmer_size != 21:
            sys.argv.extend(['--kmer-size', str(args.kmer_size)])
        if args.max_reads:
            sys.argv.extend(['--max-reads', str(args.max_reads)])
        if args.verbose:
            sys.argv.append('--verbose')
        
        # Add other arguments as needed...
        
        # Import and run the main function
        from metagrouper import main
        main()
        
        return True
        
    except Exception as e:
        logging.error(f"MetaGrouper execution failed: {e}")
        return False

if __name__ == '__main__':
    enhanced_main()