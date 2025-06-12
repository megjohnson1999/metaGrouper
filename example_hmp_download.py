#!/usr/bin/env python3
"""
Example: Download and preprocess HMP data for MetaGrouper analysis

This script demonstrates how to:
1. Download a subset of Human Microbiome Project (HMP) data
2. Preprocess the data using MetaGrouper's preprocessing pipeline
3. Run MetaGrouper analysis with metadata

Note: Requires SRA Toolkit to be installed.
Install with: conda install -c bioconda sra-tools
"""

import subprocess
import pandas as pd
import os
from pathlib import Path
import sys


def download_hmp_subset():
    """Download a small subset of HMP data for testing."""
    
    # Example HMP samples: gut vs oral
    samples = {
        # Gut samples (stool)
        "SRS011086": {"body_site": "stool", "subject": "159005608", "visit": "1"},
        "SRS011134": {"body_site": "stool", "subject": "159005725", "visit": "1"},
        "SRS011271": {"body_site": "stool", "subject": "159005903", "visit": "1"},
        "SRS011239": {"body_site": "stool", "subject": "159005873", "visit": "1"},
        "SRS011405": {"body_site": "stool", "subject": "159006033", "visit": "1"},
        
        # Oral samples (tongue dorsum)
        "SRS011090": {"body_site": "tongue_dorsum", "subject": "159005608", "visit": "1"},
        "SRS011138": {"body_site": "tongue_dorsum", "subject": "159005725", "visit": "1"},
        "SRS011275": {"body_site": "tongue_dorsum", "subject": "159005903", "visit": "1"},
        "SRS011243": {"body_site": "tongue_dorsum", "subject": "159005873", "visit": "1"},
        "SRS011409": {"body_site": "tongue_dorsum", "subject": "159006033", "visit": "1"},
    }
    
    print("📥 Downloading HMP samples...")
    print(f"   {len([s for s in samples.values() if s['body_site'] == 'stool'])} gut samples")
    print(f"   {len([s for s in samples.values() if s['body_site'] == 'tongue_dorsum'])} oral samples")
    
    # Create directories
    raw_dir = Path("hmp_raw_data")
    raw_dir.mkdir(exist_ok=True)
    
    # Download samples
    failed_downloads = []
    for sample_id, metadata in samples.items():
        print(f"   Downloading {sample_id} ({metadata['body_site']})...")
        
        try:
            # First, prefetch the SRA file
            print(f"     Prefetching {sample_id}...")
            prefetch_cmd = ["prefetch", sample_id]
            subprocess.run(prefetch_cmd, capture_output=True, text=True, check=True)
            
            # Then extract FASTQ from the prefetched SRA file
            print(f"     Extracting FASTQ for {sample_id}...")
            fastq_cmd = [
                "fastq-dump", 
                "--split-files",  # Split paired-end reads
                "--outdir", str(raw_dir),
                sample_id
            ]
            
            result = subprocess.run(fastq_cmd, capture_output=True, text=True, check=True)
            print(f"     ✅ {sample_id} downloaded")
            
        except subprocess.CalledProcessError as e:
            print(f"     ❌ {sample_id} failed: {e}")
            failed_downloads.append(sample_id)
            continue
        except FileNotFoundError:
            print("❌ fastq-dump not found. Install SRA toolkit:")
            print("   conda install -c bioconda sra-tools")
            sys.exit(1)
    
    # Remove failed samples from metadata
    for failed in failed_downloads:
        del samples[failed]
    
    # Create metadata file
    metadata_df = pd.DataFrame([
        {"sample_id": sample_id, **metadata} 
        for sample_id, metadata in samples.items()
    ])
    
    metadata_file = "hmp_metadata.csv"
    metadata_df.to_csv(metadata_file, index=False)
    
    print(f"✅ Downloaded {len(samples)} samples to {raw_dir}")
    print(f"📋 Metadata saved to {metadata_file}")
    
    return str(raw_dir), metadata_file


def preprocess_data(raw_dir: str):
    """Preprocess the downloaded data."""
    
    print("\n🔄 Preprocessing data...")
    
    # Check if preprocess.py is available
    if not Path("preprocess.py").exists():
        print("❌ preprocess.py not found in current directory")
        print("   Make sure you're running from the MetaGrouper directory")
        sys.exit(1)
    
    # Run preprocessing
    cmd = [
        "python", "preprocess.py",
        raw_dir,
        "-o", "hmp_preprocessed",
        "--quick",  # Quick mode for testing
        "-v"  # Verbose
    ]
    
    print(f"   Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True)
        print("✅ Preprocessing completed")
        return "hmp_preprocessed"
    except subprocess.CalledProcessError as e:
        print(f"❌ Preprocessing failed: {e}")
        sys.exit(1)


def run_metagrouper(preprocessed_dir: str, metadata_file: str):
    """Run MetaGrouper analysis."""
    
    print("\n🧬 Running MetaGrouper analysis...")
    
    cmd = [
        "python", "metagrouper.py",
        preprocessed_dir,
        "-m", metadata_file,
        "-o", "hmp_results",
        "--kmer-size", "15",  # Smaller k-mer for quick analysis
        "--max-reads", "5000",  # Limit reads for testing
        "--permutations", "99",  # Fewer permutations for speed
        "-v"
    ]
    
    print(f"   Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True)
        print("✅ MetaGrouper analysis completed")
        print("\n📊 Results saved to hmp_results/")
        print("   📈 Check distance_heatmap.png for sample clustering")
        print("   📋 Check analysis_report.md for statistical results")
        
    except subprocess.CalledProcessError as e:
        print(f"❌ MetaGrouper analysis failed: {e}")
        sys.exit(1)


def main():
    """Run the complete HMP example workflow."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Download and analyze HMP data with MetaGrouper")
    parser.add_argument("--auto", action="store_true", help="Run automatically without confirmation")
    args = parser.parse_args()
    
    print("🧪 MetaGrouper HMP Example Workflow")
    print("=" * 50)
    
    print("\nThis example will:")
    print("1. Download 10 HMP samples (5 gut + 5 oral)")
    print("2. Preprocess the data (quality trimming)")
    print("3. Run MetaGrouper analysis with metadata")
    print("4. Generate visualizations and reports")
    
    if not args.auto:
        # Ask for confirmation
        response = input("\nContinue? [y/N]: ").strip().lower()
        if response not in ['y', 'yes']:
            print("Cancelled.")
            sys.exit(0)
    
    try:
        # Step 1: Download data
        raw_dir, metadata_file = download_hmp_subset()
        
        # Step 2: Preprocess
        preprocessed_dir = preprocess_data(raw_dir)
        
        # Step 3: Run MetaGrouper
        run_metagrouper(preprocessed_dir, metadata_file)
        
        print("\n🎉 Complete workflow finished successfully!")
        print("\n📋 Expected results:")
        print("   • Gut and oral samples should cluster separately")
        print("   • body_site should be the top-ranked variable")
        print("   • Grouped assembly strategy recommended")
        
    except KeyboardInterrupt:
        print("\n⚠️  Workflow interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Workflow failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()