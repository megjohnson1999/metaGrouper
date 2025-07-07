#!/usr/bin/env python3
"""
MetaGrouper Complete Demonstration

This script provides a comprehensive demonstration of MetaGrouper's capabilities
across all three phases with realistic example data and complete workflow.
"""

import os
import sys
import tempfile
import shutil
import random
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import time


def create_demo_fastq_files(output_dir: Path, n_samples: int = 6):
    """Create realistic demo FASTQ files with controlled similarity patterns."""

    print(f"🧬 Generating {n_samples} demo FASTQ files...")

    # Define base sequences for different "species"
    base_sequences = {
        "species_A": "ATGCGATCGATCGATCGATCGAATGCGATCGATCGATCGATCG",
        "species_B": "ATGCGATCGATCGATCGATCGCCTGCGATCGATCGATCGATCG",
        "species_C": "ATGCGATCGATCGATCGATCGTTTGCGATCGATCGATCGATCG",
        "common": "ATGCGATCGATCGATCGATCGAAATGCGATCGATCGATCGATCG",
    }

    # Define sample compositions (which species are present)
    sample_compositions = {
        "gut_baseline_001": {"species_A": 0.4, "species_B": 0.3, "common": 0.3},
        "gut_baseline_002": {"species_A": 0.3, "species_B": 0.4, "common": 0.3},
        "gut_treated_001": {"species_A": 0.1, "species_B": 0.6, "common": 0.3},
        "gut_treated_002": {"species_A": 0.1, "species_B": 0.7, "common": 0.2},
        "skin_baseline_001": {"species_C": 0.7, "common": 0.3},
        "skin_baseline_002": {"species_C": 0.6, "common": 0.4},
    }

    output_dir.mkdir(exist_ok=True)

    for sample_name, composition in sample_compositions.items():
        fastq_file = output_dir / f"{sample_name}.fastq"

        with open(fastq_file, "w") as f:
            read_count = 0

            for species, proportion in composition.items():
                n_reads = int(800 * proportion)  # Total ~800 reads per sample
                base_seq = base_sequences[species]

                for i in range(n_reads):
                    # Add some variation to the base sequence
                    seq = base_seq
                    if random.random() < 0.1:  # 10% chance of mutation
                        pos = random.randint(0, len(seq) - 1)
                        new_base = random.choice(["A", "T", "C", "G"])
                        seq = seq[:pos] + new_base + seq[pos + 1 :]

                    # Extend sequence to realistic length
                    while len(seq) < 150:
                        seq += random.choice(["A", "T", "C", "G"])

                    seq = seq[:150]  # Trim to 150bp

                    f.write(f"@{sample_name}_read_{read_count}\n")
                    f.write(f"{seq}\n")
                    f.write("+\n")
                    f.write("~" * len(seq) + "\n")  # High quality scores

                    read_count += 1

        print(f"  ✓ Created {sample_name}.fastq ({read_count} reads)")

    return list(sample_compositions.keys())


def create_demo_metadata(output_dir: Path, sample_names: list):
    """Create realistic metadata file for demo samples."""

    print("📋 Generating demo metadata...")

    metadata = {
        "sample_id": sample_names,
        "patient_id": ["P001", "P002", "P003", "P004", "P005", "P006"],
        "body_site": ["gut", "gut", "gut", "gut", "skin", "skin"],
        "timepoint": [
            "baseline",
            "baseline",
            "treated",
            "treated",
            "baseline",
            "baseline",
        ],
        "treatment_group": [
            "control",
            "control",
            "antibiotic",
            "antibiotic",
            "control",
            "control",
        ],
        "age_group": ["adult", "adult", "adult", "adult", "elderly", "elderly"],
        "sequencing_batch": [
            "batch1",
            "batch1",
            "batch1",
            "batch2",
            "batch2",
            "batch2",
        ],
        "reads_millions": [0.8, 0.8, 0.8, 0.8, 0.8, 0.8],
        "gc_content": [0.48, 0.47, 0.52, 0.51, 0.45, 0.46],
    }

    df = pd.DataFrame(metadata)
    metadata_file = output_dir / "sample_metadata.csv"
    df.to_csv(metadata_file, index=False)

    print(f"  ✓ Created sample_metadata.csv")
    print("     Expected patterns:")
    print("     - body_site: Strong effect (gut vs skin)")
    print("     - timepoint: Strong effect (baseline vs treated)")
    print("     - treatment_group: Moderate effect (antibiotic treatment)")
    print("     - sequencing_batch: Possible technical effect")

    return str(metadata_file)


def run_metagrouper_demo(fastq_dir: Path, metadata_file: str, output_dir: Path):
    """Run MetaGrouper with demo data."""

    print(f"\n🔬 Running MetaGrouper comprehensive analysis...")

    cmd = [
        sys.executable,
        "metagrouper.py",
        str(fastq_dir),
        "--metadata",
        metadata_file,
        "--output",
        str(output_dir),
        "--kmer-size",
        "17",  # Smaller for demo speed
        "--max-reads",
        "500",  # Limit for demo speed
        "--permutations",
        "199",  # Reduced for demo speed
        "--assembly-tools",
        "megahit",
        "spades",
        "--similarity-threshold",
        "0.25",
        "--verbose",
    ]

    start_time = time.time()

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        runtime = time.time() - start_time

        if result.returncode == 0:
            print(f"  ✅ Analysis completed successfully in {runtime:.1f}s")
            return True, result.stdout
        else:
            print(f"  ❌ Analysis failed with return code {result.returncode}")
            print(f"  Error: {result.stderr}")
            return False, result.stderr

    except subprocess.TimeoutExpired:
        print(f"  ⏰ Analysis timed out after 5 minutes")
        return False, "Timeout"
    except Exception as e:
        print(f"  💥 Analysis failed with exception: {e}")
        return False, str(e)


def analyze_demo_results(output_dir: Path):
    """Analyze and display demo results."""

    print(f"\n📊 Analyzing demo results...")

    # Check what files were generated
    files_generated = []
    for file_path in output_dir.rglob("*"):
        if file_path.is_file():
            rel_path = file_path.relative_to(output_dir)
            files_generated.append(str(rel_path))

    print(f"  📁 Generated {len(files_generated)} output files")

    # Display key results
    key_files = {
        "distance_matrix.csv": "Sample similarity matrix",
        "permanova_results.csv": "Statistical test results",
        "analysis_report.md": "Comprehensive analysis report",
        "assembly_recommendations/assembly_strategy.md": "Assembly recommendations",
    }

    print(f"\n📋 Key Results Summary:")
    print("=" * 50)

    for file_name, description in key_files.items():
        file_path = output_dir / file_name
        if file_path.exists():
            print(f"\n📄 {description}:")

            if file_name == "distance_matrix.csv":
                try:
                    df = pd.read_csv(file_path, index_col=0)
                    print(
                        f"   • Sample distances (0=identical, 1=completely different):"
                    )
                    print(
                        f"   • Average distance: {df.values[np.triu_indices_from(df.values, k=1)].mean():.3f}"
                    )
                    print(
                        f"   • Most similar pair: {df.values[np.triu_indices_from(df.values, k=1)].min():.3f}"
                    )
                    print(
                        f"   • Most different pair: {df.values[np.triu_indices_from(df.values, k=1)].max():.3f}"
                    )
                except:
                    print(f"   • Could not parse distance matrix")

            elif file_name == "permanova_results.csv":
                try:
                    df = pd.read_csv(file_path)
                    df = df.dropna(subset=["r_squared"]).sort_values(
                        "r_squared", ascending=False
                    )
                    print(f"   • Variable importance (R² = variation explained):")
                    for _, row in df.head(3).iterrows():
                        significance = (
                            "***"
                            if row["p_value"] < 0.001
                            else (
                                "**"
                                if row["p_value"] < 0.01
                                else "*" if row["p_value"] < 0.05 else ""
                            )
                        )
                        print(
                            f"     - {row['variable']}: R²={row['r_squared']:.3f}, p={row['p_value']:.3f}{significance}"
                        )
                except:
                    print(f"   • Could not parse PERMANOVA results")

            elif file_name == "analysis_report.md":
                try:
                    with open(file_path, "r") as f:
                        content = f.read()
                        # Extract recommendations section
                        if "## Recommendations" in content:
                            recs_start = content.find("## Recommendations")
                            recs_section = content[recs_start : recs_start + 500]
                            lines = recs_section.split("\n")[:8]
                            print(f"   • " + "\n   • ".join(lines[1:]))  # Skip header
                except:
                    print(f"   • Could not parse analysis report")

            elif file_name == "assembly_recommendations/assembly_strategy.md":
                try:
                    with open(file_path, "r") as f:
                        content = f.read()
                        # Extract strategy info
                        if "## Overall Strategy:" in content:
                            for line in content.split("\n"):
                                if "Overall Strategy:" in line:
                                    print(f"   • {line.replace('##', '').strip()}")
                                elif "**Confidence Score:**" in line:
                                    print(f"   • {line.replace('**', '').strip()}")
                                elif "**Rationale:**" in line:
                                    print(f"   • {line.replace('**', '').strip()}")
                                    break
                except:
                    print(f"   • Could not parse assembly recommendations")
        else:
            print(f"   ❌ {description}: File not found")

    # Show assembly commands if available
    assembly_dir = output_dir / "assembly_recommendations"
    if assembly_dir.exists():
        scripts = list(assembly_dir.glob("run_*.sh"))
        if scripts:
            print(f"\n🛠️ Generated Assembly Scripts:")
            for script in scripts:
                print(f"   • {script.name}")
                # Show first command
                try:
                    with open(script, "r") as f:
                        lines = f.readlines()
                        for line in lines:
                            if line.strip() and not line.startswith("#"):
                                print(f"     First command: {line.strip()[:80]}...")
                                break
                except:
                    pass


def run_complete_demo():
    """Run the complete MetaGrouper demonstration."""

    print("🎯 METAGROUPER COMPLETE DEMONSTRATION")
    print("=" * 60)
    print("This demo shows MetaGrouper's full capabilities with realistic data.")
    print("Expected runtime: 2-5 minutes")
    print()

    # Create temporary directory for demo
    demo_dir = Path("metagrouper_demo")
    if demo_dir.exists():
        shutil.rmtree(demo_dir)
    demo_dir.mkdir()

    try:
        # Step 1: Generate demo data
        print("📋 STEP 1: Generating Demo Data")
        print("-" * 30)

        fastq_dir = demo_dir / "fastq_files"
        sample_names = create_demo_fastq_files(fastq_dir)
        metadata_file = create_demo_metadata(demo_dir, sample_names)

        # Step 2: Run MetaGrouper
        print(f"\n🔬 STEP 2: Running MetaGrouper Analysis")
        print("-" * 40)

        output_dir = demo_dir / "results"
        success, log_output = run_metagrouper_demo(fastq_dir, metadata_file, output_dir)

        if not success:
            print(f"Demo failed during analysis step. Check logs for details.")
            return False

        # Step 3: Analyze results
        print(f"\n📊 STEP 3: Analyzing Results")
        print("-" * 30)

        analyze_demo_results(output_dir)

        # Step 4: Summary and next steps
        print(f"\n🎉 DEMO COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        print(f"📁 Demo files location: {demo_dir.absolute()}")
        print(f"📁 Results location: {output_dir.absolute()}")

        print(f"\n📋 What was demonstrated:")
        print(f"   ✅ Phase 1: K-mer profiling and similarity analysis")
        print(f"   ✅ Phase 2: Metadata analysis with PERMANOVA")
        print(f"   ✅ Phase 3: Assembly strategy recommendations")
        print(f"   ✅ Comprehensive reporting and visualization")

        print(f"\n🔍 Explore the results:")
        print(f"   • Open: {output_dir}/analysis_report.md")
        print(f"   • Open: {output_dir}/assembly_recommendations/assembly_strategy.md")
        print(f"   • View: {output_dir}/*.png (visualization plots)")
        print(f"   • Run: {output_dir}/assembly_recommendations/run_*.sh")

        print(f"\n📚 Next steps:")
        print(f"   1. Review the comprehensive reports")
        print(f"   2. Examine the visualizations")
        print(f"   3. Try with your own data")
        print(f"   4. Customize parameters for your use case")

        return True

    except Exception as e:
        print(f"💥 Demo failed with unexpected error: {e}")
        return False

    finally:
        # Optional: Clean up demo files
        response = input(f"\n🗑️  Remove demo files? [y/N]: ").strip().lower()
        if response in ["y", "yes"]:
            shutil.rmtree(demo_dir)
            print(f"   ✅ Demo files removed")
        else:
            print(f"   📁 Demo files kept at: {demo_dir.absolute()}")


def run_quick_test():
    """Run a very quick test to verify MetaGrouper is working."""

    print("⚡ QUICK FUNCTIONALITY TEST")
    print("=" * 30)

    try:
        # Test imports
        print("🔧 Testing imports...")
        import metagrouper
        import metadata_analyzer
        import assembly_recommender

        print("   ✅ All modules imported successfully")

        # Test basic functionality with tiny dataset
        print("🧪 Testing basic functionality...")

        temp_dir = Path("quick_test")
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
        temp_dir.mkdir()

        # Create minimal test data
        fastq_file = temp_dir / "test.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\nATCGATCGATCGATCGATCG\n+\n~~~~~~~~~~~~~~~~~~~~\n")
            f.write("@read2\nATCGATCGATCGATCGATCG\n+\n~~~~~~~~~~~~~~~~~~~~\n")

        # Test k-mer profiling with sourmash
        profiler = metagrouper.SourmashProfiler(k=15, processes=1)
        signature = profiler.sketch_sample(str(fastq_file), "test")
        profile = {str(h): 1 for h in signature.minhash.hashes}

        if len(profile) > 0:
            print("   ✅ K-mer profiling works")
        else:
            print("   ❌ K-mer profiling failed")
            return False

        shutil.rmtree(temp_dir)

        print("🎉 Quick test passed!")
        print("   MetaGrouper is ready to use!")

        return True

    except Exception as e:
        print(f"❌ Quick test failed: {e}")
        return False


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="MetaGrouper Demonstration")
    parser.add_argument(
        "--quick", action="store_true", help="Run quick functionality test only"
    )
    parser.add_argument(
        "--full", action="store_true", help="Run complete demonstration (default)"
    )

    args = parser.parse_args()

    if args.quick:
        success = run_quick_test()
    else:
        success = run_complete_demo()

    sys.exit(0 if success else 1)
