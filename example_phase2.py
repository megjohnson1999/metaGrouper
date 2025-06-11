#!/usr/bin/env python3
"""
Example usage script for MetaGrouper Phase 2 (Metadata Analysis)

This script demonstrates how to use MetaGrouper with metadata for
comprehensive analysis including PERMANOVA and clustering.
"""

import os
import pandas as pd
import numpy as np
import subprocess
import sys
from pathlib import Path

def generate_example_metadata():
    """Generate example metadata for the test samples."""
    
    # Create metadata for our 4 example samples
    metadata = {
        'sample_id': ['sample_A', 'sample_B', 'sample_C', 'sample_D'],
        'patient_id': ['P001', 'P001', 'P002', 'P002'],
        'timepoint': ['baseline', 'week4', 'baseline', 'week4'],
        'treatment': ['control', 'treated', 'control', 'treated'],
        'age': [25, 25, 34, 34],
        'bmi': [22.5, 23.1, 28.2, 27.8],
        'location': ['clinic_A', 'clinic_A', 'clinic_B', 'clinic_B'],
        'sequencing_depth': [1.2e6, 1.5e6, 0.9e6, 1.8e6]
    }
    
    df = pd.DataFrame(metadata)
    
    # Save metadata file
    metadata_file = Path("example_data") / "metadata.csv"
    df.to_csv(metadata_file, index=False)
    
    print(f"Generated example metadata:")
    print(df)
    print(f"\nSaved to: {metadata_file}")
    
    return str(metadata_file)

def run_phase2_example():
    """Run MetaGrouper with Phase 2 metadata analysis."""
    
    # First, ensure we have the example data from Phase 1
    example_data_dir = Path("example_data")
    if not example_data_dir.exists() or not any(example_data_dir.glob("*.fastq")):
        print("Example data not found. Running Phase 1 example first...")
        try:
            subprocess.run([sys.executable, "example_usage.py"], check=True)
        except subprocess.CalledProcessError:
            print("Failed to generate Phase 1 example data")
            return
    
    # Generate example metadata
    metadata_file = generate_example_metadata()
    
    # Run MetaGrouper with Phase 2
    print(f"\nRunning MetaGrouper with Phase 2 metadata analysis...")
    
    cmd = [
        sys.executable, "metagrouper.py",
        "example_data",
        "--output", "phase2_output",
        "--kmer-size", "15",  # Smaller k-mer for example
        "--max-reads", "500",  # Limit reads for faster processing
        "--metadata", metadata_file,
        "--sample-id-column", "sample_id",
        "--permutations", "199",  # Fewer permutations for faster testing
        "--cluster-range", "2", "4",
        "--verbose"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        print("STDOUT:")
        print(result.stdout)
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
        
        if result.returncode == 0:
            print("\nPhase 2 example completed successfully!")
            print("Check the 'phase2_output' directory for results.")
            print("\nGenerated files should include:")
            print("- analysis_report.md: Comprehensive analysis report")
            print("- variable_importance.png: Variable importance ranking")
            print("- permanova_results.csv: Statistical results")
            print("- pca_by_*.png: PCA plots colored by metadata")
            print("- clustering_*.png: Clustering results")
        else:
            print(f"Phase 2 example failed with return code {result.returncode}")
            
    except Exception as e:
        print(f"Error running Phase 2 example: {e}")

def demonstrate_custom_analysis():
    """Demonstrate custom metadata analysis."""
    print("\n" + "="*60)
    print("Custom Metadata Analysis Example")
    print("="*60)
    
    # Run with specific variables only
    metadata_file = Path("example_data") / "metadata.csv"
    if not metadata_file.exists():
        print("Metadata file not found. Run main example first.")
        return
    
    print("Running analysis for specific variables: patient_id, treatment")
    
    cmd = [
        sys.executable, "metagrouper.py",
        "example_data",
        "--output", "custom_analysis",
        "--kmer-size", "15",
        "--max-reads", "500",
        "--metadata", str(metadata_file),
        "--variables", "patient_id", "treatment",
        "--permutations", "99",  # Very fast for demo
        "--verbose"
    ]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("Custom analysis completed successfully!")
        print("Results saved to 'custom_analysis' directory.")
    except subprocess.CalledProcessError as e:
        print(f"Custom analysis failed: {e}")

if __name__ == '__main__':
    print("MetaGrouper Phase 2 Example")
    print("="*50)
    
    # Run comprehensive Phase 2 example
    run_phase2_example()
    
    # Run custom analysis example
    demonstrate_custom_analysis()
    
    print("\nPhase 2 examples complete!")
    print("\nNext steps:")
    print("1. Examine the analysis_report.md files for detailed results")
    print("2. Review the variable importance plots")
    print("3. Compare clustering results with metadata variables")
    print("4. Use the insights to plan your assembly strategy")