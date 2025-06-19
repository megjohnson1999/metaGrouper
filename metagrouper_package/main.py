#!/usr/bin/env python3
"""
MetaGrouper: Clean main entry point with Phase 1 integration.
"""

import sys
import argparse
import logging
import time
from pathlib import Path

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))

from metagrouper import (
    find_fastq_files,
    setup_logging,
    save_results,
)

# Phase 1 imports
try:
    from metagrouper.sketch_profiler import StreamingKmerProfiler
    from metagrouper.sparse_analyzer import SparseSimilarityAnalyzer
    from metagrouper.profiler import KmerProfiler
    PHASE1_AVAILABLE = True
except ImportError as e:
    PHASE1_AVAILABLE = False
    logging.error(f"Phase 1 not available: {e}")


def create_parser():
    """Create command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="MetaGrouper: K-mer analysis for metagenomic assembly grouping",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  python main_clean.py /path/to/fastq/files -o results/
  
  # Memory-efficient with sketching (Phase 1)
  python main_clean.py /path/to/fastq/files -o results/ --use-sketching --sketch-size 2000
  
  # Large dataset with sparse analysis
  python main_clean.py /path/to/fastq/files -o results/ --use-sketching --similarity-threshold 0.15
        """
    )
    
    # Required arguments
    parser.add_argument("input_dir", help="Directory containing FASTQ files")
    parser.add_argument("-o", "--output", default="metagrouper_output", 
                       help="Output directory (default: metagrouper_output)")
    
    # Core k-mer arguments
    parser.add_argument("-k", "--kmer-size", type=int, default=21, 
                       help="K-mer size (default: 21)")
    parser.add_argument("--max-reads", type=int, 
                       help="Maximum reads per sample (for testing)")
    parser.add_argument("--min-kmer-freq", type=int, default=1,
                       help="Minimum k-mer frequency (default: 1)")
    
    # Phase 1 arguments
    parser.add_argument("--use-sketching", action="store_true",
                       help="Use streaming k-mer sketches for memory efficiency")
    parser.add_argument("--sketch-size", type=int, default=1000,
                       help="K-mer sketch size (default: 1000)")
    parser.add_argument("--sampling-method", choices=["reservoir", "frequency", "adaptive"],
                       default="frequency", help="Sketching method (default: frequency)")
    parser.add_argument("--similarity-threshold", type=float, default=0.1,
                       help="Sparse similarity threshold (default: 0.1)")
    
    # Processing arguments
    parser.add_argument("--processes", type=int,
                       help="Number of parallel processes (default: CPU count)")
    parser.add_argument("--sequential", action="store_true",
                       help="Use sequential processing")
    
    # Analysis arguments
    parser.add_argument("--distance-metric", default="jaccard",
                       choices=["jaccard", "weighted_jaccard", "cosine"],
                       help="Similarity metric (default: jaccard)")
    
    # Other arguments
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    
    return parser


def run_analysis(args):
    """Run the main analysis workflow."""
    # Setup logging
    setup_logging(args.verbose)
    
    print("🧬 MetaGrouper Analysis")
    print("=" * 50)
    
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
    
    # Choose profiler based on options
    if args.use_sketching and PHASE1_AVAILABLE:
        print(f"🚀 Using Phase 1 optimizations (sketching + sparse analysis)")
        print(f"   Sketch size: {args.sketch_size}")
        print(f"   Sampling method: {args.sampling_method}")
        print(f"   Similarity threshold: {args.similarity_threshold}")
        
        profiler = StreamingKmerProfiler(
            k=args.kmer_size,
            max_reads=args.max_reads,
            min_kmer_freq=args.min_kmer_freq,
            sketch_size=args.sketch_size,
            sampling_method=args.sampling_method
        )
    else:
        print(f"📊 Using traditional k-mer profiling")
        profiler = KmerProfiler(
            k=args.kmer_size,
            max_reads=args.max_reads,
            min_kmer_freq=args.min_kmer_freq
        )
    
    # Process samples
    print(f"\n🔬 Processing {len(fastq_files)} samples...")
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
    
    # Similarity analysis
    print(f"\n🔗 Computing similarity matrix...")
    start_time = time.time()
    
    if args.use_sketching and PHASE1_AVAILABLE:
        # Use sparse similarity analysis
        analyzer = SparseSimilarityAnalyzer(similarity_threshold=args.similarity_threshold)
        similarity_matrix, sample_names = analyzer.compute_similarities(profiles, method=args.distance_metric)
        
        # Get sparse statistics
        sparse_stats = analyzer.compute_summary_statistics()
        similarity_time = time.time() - start_time
        
        print(f"✅ Sparse similarity computed in {similarity_time:.1f}s")
        print(f"📊 Sparsity: {sparse_stats['sparsity']:.1%}")
        print(f"📊 Significant pairs: {sparse_stats['n_significant_pairs']:,}")
        
        # Save sparse matrix
        analyzer.save_similarity_matrix(str(output_path / "similarity_matrix.npz"), format='npz')
        
    else:
        # Traditional dense similarity (for small datasets)
        from metagrouper.analyzer import SimilarityAnalyzer
        analyzer = SimilarityAnalyzer(profiles)
        distance_matrix = analyzer.compute_distance_matrix(metric=args.distance_metric)
        similarity_matrix = 1 - distance_matrix
        sample_names = list(profiles.keys())
        
        similarity_time = time.time() - start_time
        print(f"✅ Dense similarity computed in {similarity_time:.1f}s")
        
        # Save results
        save_results(profiles, distance_matrix, sample_names, args.output)
    
    # Generate basic visualization
    try:
        from metagrouper.visualizer import Visualizer
        visualizer = Visualizer()
        
        if args.use_sketching:
            # Convert sparse to dense for visualization (small matrices only)
            if len(sample_names) <= 100:
                dense_similarity = similarity_matrix.toarray()
                visualizer.plot_distance_heatmap(1 - dense_similarity, sample_names, 
                                                str(output_path / "similarity_heatmap.png"))
        else:
            visualizer.plot_distance_heatmap(distance_matrix, sample_names,
                                           str(output_path / "distance_heatmap.png"))
        
        print(f"📈 Visualization saved to {output_path}")
        
    except Exception as e:
        logging.warning(f"Visualization failed: {e}")
    
    # Summary
    total_time = processing_time + similarity_time
    print(f"\n🎯 Analysis Complete!")
    print(f"   Total time: {total_time:.1f}s")
    print(f"   Samples processed: {success_count}")
    print(f"   Results saved to: {args.output}")
    
    if args.use_sketching:
        print(f"   Memory efficiency: {memory_info['memory_reduction']:.1f}x improvement")
        print(f"   Sparsity benefit: {sparse_stats['sparsity']:.1%} memory saved")
    
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