#!/usr/bin/env python3
"""
Example usage script for MetaGrouper Phase 1

This script demonstrates how to use MetaGrouper with example data.
"""

import os
import tempfile
import random
from pathlib import Path

def generate_example_fastq(filepath: str, num_reads: int = 1000, read_length: int = 150):
    """Generate a simple example FASTQ file for testing."""
    bases = ['A', 'T', 'C', 'G']
    
    with open(filepath, 'w') as f:
        for i in range(num_reads):
            # Generate random DNA sequence
            sequence = ''.join(random.choices(bases, k=read_length))
            quality = '~' * read_length  # High quality scores
            
            f.write(f"@read_{i+1}\n")
            f.write(f"{sequence}\n")
            f.write("+\n")
            f.write(f"{quality}\n")

def create_example_dataset():
    """Create example dataset with 4 samples."""
    # Create temporary directory for example data
    example_dir = Path("example_data")
    example_dir.mkdir(exist_ok=True)
    
    print("Generating example FASTQ files...")
    
    # Generate 4 example samples
    samples = [
        ("sample_A.fastq", 1000),
        ("sample_B.fastq", 1200),
        ("sample_C.fastq", 800),
        ("sample_D.fastq", 1500)
    ]
    
    for sample_name, num_reads in samples:
        filepath = example_dir / sample_name
        generate_example_fastq(str(filepath), num_reads)
        print(f"  Generated {sample_name} with {num_reads} reads")
    
    return str(example_dir)

def run_example():
    """Run MetaGrouper on example data."""
    import subprocess
    import sys
    
    # Create example dataset
    input_dir = create_example_dataset()
    
    # Run MetaGrouper
    print(f"\nRunning MetaGrouper on example data...")
    cmd = [
        sys.executable, "metagrouper.py",
        input_dir,
        "--output", "example_output",
        "--kmer-size", "15",  # Smaller k-mer for example
        "--max-reads", "500",  # Limit reads for faster processing
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
            print("\nExample completed successfully!")
            print("Check the 'example_output' directory for results.")
        else:
            print(f"Example failed with return code {result.returncode}")
            
    except Exception as e:
        print(f"Error running example: {e}")

if __name__ == '__main__':
    # Set random seed for reproducible results
    random.seed(42)
    run_example()