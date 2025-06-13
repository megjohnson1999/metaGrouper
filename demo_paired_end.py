#!/usr/bin/env python3
"""
Demonstration of MetaGrouper's improved paired-end support.

This script demonstrates the new paired-end read detection and processing
capabilities added in Phase 1 improvements.
"""

import tempfile
import logging
from pathlib import Path
from metagrouper import find_fastq_files, KmerProfiler

def create_demo_files():
    """Create demonstration FASTQ files."""
    temp_dir = Path(tempfile.mkdtemp())
    print(f"Creating demo files in: {temp_dir}")
    
    # Create paired-end samples
    samples = {
        "gut_sample_001": {
            "R1": ["ATCGATCGATCGATCGATCG", "GCTAGCTAGCTAGCTAGCTA", "TTCCGGAATTCCGGAATTCC"],
            "R2": ["CGTATCGATCGATCGATCGT", "TAGATGATGATGATGATGAT", "GGAATTCCGGAATTCCGGAA"]
        },
        "soil_sample_002": {
            "R1": ["AAAATTTCCCGGGTTTAAAA", "CCCCGGGGAAAATTTTCCCC"],
            "R2": ["TTTTAAAAGGGCCCAAATTT", "GGGGTTTTCCCCAAAAGGGG"]
        }
    }
    
    # Create single-end sample
    single_end = {
        "water_sample_003": ["ATGCATGCATGCATGCATGC", "CGATCGATCGATCGATCGAT"]
    }
    
    # Write paired-end files
    for sample_name, reads in samples.items():
        for direction in ["R1", "R2"]:
            filename = f"{sample_name}_{direction}.fastq"
            filepath = temp_dir / filename
            with open(filepath, "w") as f:
                for i, seq in enumerate(reads[direction]):
                    f.write(f"@{sample_name}_{direction}_read_{i}\n")
                    f.write(f"{seq}\n")
                    f.write("+\n")
                    f.write("~" * len(seq) + "\n")
    
    # Write single-end file
    for sample_name, reads in single_end.items():
        filename = f"{sample_name}.fastq"
        filepath = temp_dir / filename
        with open(filepath, "w") as f:
            for i, seq in enumerate(reads):
                f.write(f"@{sample_name}_read_{i}\n")
                f.write(f"{seq}\n")
                f.write("+\n")
                f.write("~" * len(seq) + "\n")
    
    return temp_dir

def main():
    """Demonstrate paired-end functionality."""
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    
    print("=== MetaGrouper Paired-End Support Demo ===\n")
    
    # Create demo files
    demo_dir = create_demo_files()
    
    # Demonstrate file detection
    print("1. File Detection:")
    print("-" * 50)
    files = find_fastq_files(str(demo_dir))
    
    for filepath, sample_name in files:
        if isinstance(filepath, list):
            print(f"✓ Paired-end sample: {sample_name}")
            print(f"  - R1: {Path(filepath[0]).name}")
            print(f"  - R2: {Path(filepath[1]).name}")
        else:
            print(f"✓ Single-end sample: {sample_name}")
            print(f"  - File: {Path(filepath).name}")
        print()
    
    # Demonstrate processing
    print("2. K-mer Profiling:")
    print("-" * 50)
    profiler = KmerProfiler(k=17, max_reads=10)
    
    for filepath, sample_name in files:
        print(f"Processing {sample_name}...")
        try:
            profile = profiler.profile_sample(filepath, sample_name)
            print(f"  ✓ Generated {len(profile)} unique k-mers")
            print(f"  ✓ Profile normalized (sum={sum(profile.values()):.3f})")
        except Exception as e:
            print(f"  ✗ Error: {e}")
        print()
    
    # Show benefits
    print("3. Benefits of Paired-End Support:")
    print("-" * 50)
    print("• Processes both R1 and R2 reads for paired-end samples")
    print("• Automatically detects multiple naming conventions (_R1/_R2, _1/_2, etc.)")
    print("• Gracefully handles orphaned mates (single R1 or R2 files)")
    print("• Maintains backward compatibility with single-end files")
    print("• Improved k-mer diversity from combined read information")
    print()
    
    print("4. Validation Features:")
    print("-" * 50)
    print("• FASTQ format validation with detailed error messages")
    print("• Quality score length validation")
    print("• DNA sequence validation (warns on non-standard bases)")
    print("• Robust error handling with context information")
    
    # Cleanup
    import shutil
    shutil.rmtree(demo_dir)
    print(f"\nDemo completed! (cleaned up {demo_dir})")

if __name__ == "__main__":
    main()