#!/usr/bin/env python3
"""
Demonstration of the new modular MetaGrouper architecture.

This script showcases the improved modular structure, configuration management,
and enhanced functionality introduced in Phase 3.
"""

import tempfile
import logging
import shutil
import sys
from pathlib import Path

# Add the package to the path
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
from metagrouper.utils import (
    summarize_fastq_files,
    estimate_memory_usage,
    get_system_info,
    check_dependencies
)


def create_demo_dataset():
    """Create a demonstration dataset with diverse samples."""
    temp_dir = Path(tempfile.mkdtemp())
    print(f"📁 Creating demo dataset in: {temp_dir}")
    
    # Define sample compositions (different microbial profiles)
    samples = {
        "gut_healthy": {
            "sequences": [
                "ATGAAAGCGTGGCTGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
                "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
                "TTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGG"
            ],
            "type": "paired"
        },
        "gut_diseased": {
            "sequences": [
                "CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "TATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTAT",
                "AAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGG"
            ],
            "type": "paired"
        },
        "soil_agricultural": {
            "sequences": [
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG",
                "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
                "GAATTCGAATTCGAATTCGAATTCGAATTCGAATTCGAATTCGAATTCGAATTC"
            ],
            "type": "single"
        },
        "marine_surface": {
            "sequences": [
                "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG",
                "ATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATA",
                "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
            ],
            "type": "paired"
        },
        "freshwater": {
            "sequences": [
                "TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATT",
                "CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC",
                "AGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTA"
            ],
            "type": "single"
        }
    }
    
    # Create FASTQ files
    for sample_name, data in samples.items():
        if data["type"] == "paired":
            # Create paired-end files
            for direction in ["R1", "R2"]:
                filename = f"{sample_name}_{direction}.fastq"
                filepath = temp_dir / filename
                with open(filepath, "w") as f:
                    for i, seq in enumerate(data["sequences"]):
                        # Add some variation between R1 and R2
                        if direction == "R2":
                            # Simulate reverse complement for some reads
                            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                            seq = ''.join(complement.get(base, base) for base in reversed(seq))
                        
                        f.write(f"@{sample_name}_{direction}_read_{i}\n")
                        f.write(f"{seq}\n")
                        f.write("+\n")
                        f.write("~" * len(seq) + "\n")
        else:
            # Create single-end file
            filename = f"{sample_name}.fastq"
            filepath = temp_dir / filename
            with open(filepath, "w") as f:
                for i, seq in enumerate(data["sequences"]):
                    f.write(f"@{sample_name}_read_{i}\n")
                    f.write(f"{seq}\n")
                    f.write("+\n")
                    f.write("~" * len(seq) + "\n")
    
    return temp_dir


def demonstrate_configuration_system():
    """Demonstrate the new configuration system."""
    print("\n⚙️  Configuration System Demonstration")
    print("-" * 50)
    
    # Create default configuration
    print("1. Default Configuration:")
    config = MetaGrouperConfig()
    print("   ✓ Created with sensible defaults")
    
    # Show configuration summary
    print("\n2. Configuration Summary:")
    summary_lines = config.get_summary().split('\n')
    for line in summary_lines[:10]:  # Show first 10 lines
        print(f"   {line}")
    print("   ...")
    
    # Demonstrate configuration modification
    print("\n3. Configuration Modification:")
    config.profiling.k_size = 19
    config.profiling.min_kmer_freq = 2
    config.analysis.distance_metric = "jaccard"
    print("   ✓ Modified k-mer size, frequency threshold, and distance metric")
    
    # Demonstrate validation
    print("\n4. Configuration Validation:")
    try:
        config.validate()
        print("   ✓ Configuration is valid")
    except ValueError as e:
        print(f"   ❌ Configuration error: {e}")
    
    return config


def demonstrate_utility_functions():
    """Demonstrate the new utility functions."""
    print("\n🔧 Utility Functions Demonstration")
    print("-" * 50)
    
    # System information
    print("1. System Information:")
    sys_info = get_system_info()
    print(f"   Platform: {sys_info['platform']}")
    print(f"   Python: {sys_info['python_version']}")
    print(f"   CPU cores: {sys_info['cpu_count']}")
    
    # Dependency checking
    print("\n2. Dependency Check:")
    deps = check_dependencies()
    for dep, available in deps.items():
        status = "✓" if available else "❌"
        print(f"   {status} {dep}")
    
    # Memory estimation
    print("\n3. Memory Usage Estimation:")
    memory_est = estimate_memory_usage(5, 100, 21)
    print(f"   Estimated k-mers: {memory_est['estimated_kmers']:,}")
    print(f"   Profile storage: {memory_est['profile_storage_mb']:.1f} MB")
    print(f"   Distance matrix: {memory_est['distance_matrix_mb']:.1f} MB")
    print(f"   Total estimated: {memory_est['total_estimated_mb']:.1f} MB")


def demonstrate_enhanced_file_discovery(test_dir):
    """Demonstrate enhanced file discovery capabilities."""
    print("\n📂 Enhanced File Discovery")
    print("-" * 50)
    
    # Find files
    fastq_files = find_fastq_files(str(test_dir))
    
    print("1. Discovered Files:")
    for filepath, sample_name in fastq_files:
        if isinstance(filepath, list):
            print(f"   📌 Paired-end: {sample_name}")
            print(f"      R1: {Path(filepath[0]).name}")
            print(f"      R2: {Path(filepath[1]).name}")
        else:
            print(f"   📄 Single-end: {sample_name}")
            print(f"      File: {Path(filepath).name}")
    
    # File summary
    print("\n2. File Summary:")
    summary = summarize_fastq_files(fastq_files)
    print(f"   Total samples: {summary['total_samples']}")
    print(f"   Paired-end samples: {summary['paired_end_samples']}")
    print(f"   Single-end samples: {summary['single_end_samples']}")
    print(f"   Total files: {summary['total_files']}")
    print(f"   Total size: {summary['total_size_formatted']}")
    
    return fastq_files


def demonstrate_modular_processing(config, fastq_files):
    """Demonstrate the modular processing pipeline."""
    print("\n🧬 Modular Processing Pipeline")
    print("-" * 50)
    
    # Initialize profiler with configuration
    print("1. Initializing K-mer Profiler:")
    profiler = KmerProfiler(
        k=config.profiling.k_size,
        max_reads=config.profiling.max_reads,
        min_kmer_freq=config.profiling.min_kmer_freq
    )
    print(f"   ✓ K-mer size: {profiler.k}")
    print(f"   ✓ Min frequency: {profiler.min_kmer_freq}")
    
    # Process samples
    print("\n2. Processing Samples:")
    for filepath, sample_name in fastq_files:
        print(f"   Processing {sample_name}...")
        try:
            profile = profiler.profile_sample(
                filepath, 
                sample_name, 
                memory_efficient=config.profiling.memory_efficient
            )
            print(f"   ✓ Generated {len(profile)} unique k-mers")
        except Exception as e:
            print(f"   ❌ Error: {e}")
    
    print(f"\n   📊 Successfully profiled {len(profiler.profiles)} samples")
    return profiler


def demonstrate_modular_analysis(profiler, config):
    """Demonstrate the modular analysis pipeline."""
    print("\n📈 Modular Analysis Pipeline")
    print("-" * 50)
    
    # Initialize analyzer
    print("1. Initializing Similarity Analyzer:")
    analyzer = SimilarityAnalyzer(
        profiler.profiles, 
        memory_efficient=config.analysis.memory_efficient
    )
    print(f"   ✓ Loaded {len(analyzer.sample_names)} sample profiles")
    print(f"   ✓ Memory efficient mode: {analyzer.memory_efficient}")
    
    # Compute distance matrix
    print("\n2. Computing Distance Matrix:")
    distance_matrix = analyzer.compute_distance_matrix(metric=config.analysis.distance_metric)
    print(f"   ✓ Matrix shape: {distance_matrix.shape}")
    print(f"   ✓ Metric: {config.analysis.distance_metric}")
    print(f"   ✓ Mean distance: {distance_matrix.mean():.3f}")
    
    # Perform dimensionality reduction
    print("\n3. Dimensionality Reduction:")
    pca_result, pca = analyzer.perform_pca()
    print(f"   ✓ PCA completed: {pca_result.shape}")
    print(f"   ✓ Explained variance: {pca.explained_variance_ratio_[0]:.1%}, {pca.explained_variance_ratio_[1]:.1%}")
    
    mds_result = analyzer.perform_mds()
    print(f"   ✓ MDS completed: {mds_result.shape}")
    
    return analyzer, distance_matrix, pca_result, pca, mds_result


def demonstrate_enhanced_visualization(profiler, distance_matrix, pca_result, pca, mds_result, output_dir):
    """Demonstrate enhanced visualization capabilities."""
    print("\n🎨 Enhanced Visualization")
    print("-" * 50)
    
    # Initialize visualizer
    visualizer = Visualizer(profiler.sample_names)
    
    print("1. Generating Standard Plots:")
    visualizer.plot_distance_heatmap(distance_matrix, str(output_dir / "heatmap.png"))
    print("   ✓ Distance heatmap")
    
    visualizer.plot_pca(pca_result, pca, str(output_dir / "pca.png"))
    print("   ✓ PCA plot")
    
    visualizer.plot_mds(mds_result, str(output_dir / "mds.png"))
    print("   ✓ MDS plot")
    
    print("\n2. Generating Enhanced Plots:")
    visualizer.plot_sample_overview(
        distance_matrix, pca_result, mds_result, pca, 
        str(output_dir / "overview.png")
    )
    print("   ✓ Comprehensive overview")
    
    visualizer.plot_clustering_dendrogram(distance_matrix, str(output_dir / "dendrogram.png"))
    print("   ✓ Hierarchical clustering dendrogram")
    
    visualizer.plot_distance_distribution(distance_matrix, str(output_dir / "distances.png"))
    print("   ✓ Distance distribution")


def demonstrate_result_management(profiler, distance_matrix, output_dir, config):
    """Demonstrate result saving and configuration management."""
    print("\n💾 Result Management")
    print("-" * 50)
    
    # Save analysis results
    print("1. Saving Analysis Results:")
    save_results(
        profiler.profiles, 
        distance_matrix, 
        profiler.sample_names, 
        str(output_dir)
    )
    print("   ✓ K-mer profiles saved")
    print("   ✓ Distance matrix saved")
    print("   ✓ Sample names saved")
    
    # Save configuration
    print("\n2. Saving Configuration:")
    config.save_to_file(str(output_dir / "analysis_config.json"))
    print("   ✓ Configuration saved for reproducibility")
    
    # List all outputs
    print("\n3. Generated Files:")
    for file_path in sorted(output_dir.iterdir()):
        file_size = file_path.stat().st_size
        print(f"   📄 {file_path.name} ({file_size:,} bytes)")


def main():
    """Main demonstration workflow."""
    setup_logging(verbose=False)
    
    print("🎯 MetaGrouper Phase 3: Modular Architecture Demo")
    print("=" * 60)
    
    # Create demo dataset
    test_dir = create_demo_dataset()
    output_dir = Path(tempfile.mkdtemp()) / "modular_demo_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Demonstrate configuration system
        config = demonstrate_configuration_system()
        
        # Demonstrate utility functions
        demonstrate_utility_functions()
        
        # Demonstrate file discovery
        fastq_files = demonstrate_enhanced_file_discovery(test_dir)
        
        # Demonstrate modular processing
        profiler = demonstrate_modular_processing(config, fastq_files)
        
        # Demonstrate modular analysis
        analyzer, distance_matrix, pca_result, pca, mds_result = demonstrate_modular_analysis(profiler, config)
        
        # Demonstrate enhanced visualization
        demonstrate_enhanced_visualization(profiler, distance_matrix, pca_result, pca, mds_result, output_dir)
        
        # Demonstrate result management
        demonstrate_result_management(profiler, distance_matrix, output_dir, config)
        
        # Summary
        print("\n🎉 Phase 3 Modular Architecture Benefits")
        print("=" * 60)
        print("✅ Modular Design:")
        print("   • Separated concerns into focused modules")
        print("   • Improved code maintainability and testability")
        print("   • Easy to extend with new functionality")
        
        print("\n✅ Configuration Management:")
        print("   • Centralized configuration system")
        print("   • JSON-based configuration files")
        print("   • Environment variable support")
        print("   • Parameter validation and recommendations")
        
        print("\n✅ Enhanced Utilities:")
        print("   • System information gathering")
        print("   • Dependency checking")
        print("   • Memory usage estimation")
        print("   • Comprehensive file handling")
        
        print("\n✅ Improved Visualizations:")
        print("   • Additional plot types")
        print("   • Comprehensive overview plots")
        print("   • Better statistical summaries")
        print("   • Enhanced formatting and layouts")
        
        print("\n✅ Better User Experience:")
        print("   • Clear progress reporting")
        print("   • Informative error messages")
        print("   • Comprehensive result summaries")
        print("   • Backwards compatibility")
        
        print(f"\n📁 Demo outputs saved to: {output_dir}")
        
    finally:
        # Cleanup
        shutil.rmtree(test_dir)
        print(f"\n🧹 Cleaned up temporary files")


if __name__ == "__main__":
    main()