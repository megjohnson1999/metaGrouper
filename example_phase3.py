#!/usr/bin/env python3
"""
Example usage script for MetaGrouper Phase 3 (Assembly Recommendations)

This script demonstrates how to use MetaGrouper with all three phases
to get comprehensive assembly strategy recommendations.
"""

import os
import pandas as pd
import numpy as np
import subprocess
import sys
from pathlib import Path


def generate_realistic_metadata():
    """Generate more realistic metadata that creates meaningful groupings."""

    # Create metadata that should result in clear assembly recommendations
    metadata = {
        "sample_id": ["sample_A", "sample_B", "sample_C", "sample_D"],
        "patient_id": ["P001", "P001", "P002", "P003"],  # Two samples from same patient
        "timepoint": ["baseline", "week4", "baseline", "baseline"],
        "treatment_group": [
            "control",
            "control",
            "treatment",
            "treatment",
        ],  # Clear grouping
        "sampling_site": ["gut", "gut", "gut", "skin"],  # Different body sites
        "age_group": ["young", "young", "adult", "adult"],
        "antibiotic_exposure": ["none", "low", "high", "none"],
        "sequencing_batch": [
            "batch1",
            "batch1",
            "batch2",
            "batch2",
        ],  # Technical grouping
        "reads_millions": [1.2, 1.5, 0.9, 1.8],
        "gc_content": [0.45, 0.44, 0.52, 0.51],
    }

    df = pd.DataFrame(metadata)

    # Save metadata file
    metadata_file = Path("example_data") / "realistic_metadata.csv"
    df.to_csv(metadata_file, index=False)

    print(f"Generated realistic metadata for assembly recommendations:")
    print(df)
    print(f"\nExpected groupings:")
    print("- treatment_group: control vs treatment samples")
    print("- sampling_site: gut vs skin samples")
    print("- sequencing_batch: technical replicates")
    print(f"\nSaved to: {metadata_file}")

    return str(metadata_file)


def run_comprehensive_example():
    """Run MetaGrouper with all three phases for comprehensive analysis."""

    # Ensure we have example data
    example_data_dir = Path("example_data")
    if not example_data_dir.exists() or not any(example_data_dir.glob("*.fastq")):
        print("Example data not found. Running Phase 1 example first...")
        try:
            subprocess.run([sys.executable, "example_usage.py"], check=True)
        except subprocess.CalledProcessError:
            print("Failed to generate Phase 1 example data")
            return

    # Generate realistic metadata
    metadata_file = generate_realistic_metadata()

    # Run comprehensive MetaGrouper analysis
    print(f"\nRunning comprehensive MetaGrouper analysis (all phases)...")

    cmd = [
        sys.executable,
        "metagrouper.py",
        "example_data",
        "--output",
        "comprehensive_analysis",
        "--kmer-size",
        "15",
        "--max-reads",
        "500",
        "--metadata",
        metadata_file,
        "--sample-id-column",
        "sample_id",
        "--permutations",
        "199",  # Faster for demo
        "--cluster-range",
        "2",
        "4",
        "--assembly-tools",
        "megahit",
        "spades",
        "--similarity-threshold",
        "0.25",  # More stringent
        "--min-group-size",
        "2",
        "--max-group-size",
        "6",
        "--verbose",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        print("STDOUT:")
        print(result.stdout)
        if result.stderr:
            print("STDERR:")
            print(result.stderr)

        if result.returncode == 0:
            print("\nComprehensive analysis completed successfully!")
            print("Check the 'comprehensive_analysis' directory for all results.")

            # Show what was generated
            print("\nGenerated outputs:")
            output_dir = Path("comprehensive_analysis")
            if output_dir.exists():
                all_files = list(output_dir.rglob("*"))
                for file_path in sorted(all_files):
                    if file_path.is_file():
                        rel_path = file_path.relative_to(output_dir)
                        print(f"  - {rel_path}")
        else:
            print(f"Comprehensive analysis failed with return code {result.returncode}")

    except Exception as e:
        print(f"Error running comprehensive analysis: {e}")


def demonstrate_assembly_tools():
    """Demonstrate different assembly tool configurations."""
    print("\n" + "=" * 60)
    print("Assembly Tools Comparison Example")
    print("=" * 60)

    metadata_file = Path("example_data") / "realistic_metadata.csv"
    if not metadata_file.exists():
        print("Metadata file not found. Run main example first.")
        return

    # Test different assembly tool combinations
    tool_configs = [
        (["megahit"], "megahit_only"),
        (["spades"], "spades_only"),
        (["flye"], "flye_only"),
        (["all"], "all_tools"),
    ]

    for tools, output_suffix in tool_configs:
        print(f"\nTesting assembly tools: {tools}")

        cmd = (
            [
                sys.executable,
                "metagrouper.py",
                "example_data",
                "--output",
                f"assembly_test_{output_suffix}",
                "--kmer-size",
                "15",
                "--max-reads",
                "300",  # Faster
                "--metadata",
                str(metadata_file),
                "--assembly-tools",
            ]
            + tools
            + ["--permutations", "99", "--similarity-threshold", "0.30"]  # Very fast
        )

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"✓ {output_suffix} completed successfully")

            # Check what assembly scripts were generated
            scripts_dir = (
                Path(f"assembly_test_{output_suffix}") / "assembly_recommendations"
            )
            if scripts_dir.exists():
                scripts = list(scripts_dir.glob("run_*.sh"))
                print(f"  Generated {len(scripts)} assembly scripts")

        except subprocess.CalledProcessError as e:
            print(f"✗ {output_suffix} failed: {e}")


def demonstrate_threshold_sensitivity():
    """Demonstrate how different similarity thresholds affect recommendations."""
    print("\n" + "=" * 60)
    print("Similarity Threshold Sensitivity Analysis")
    print("=" * 60)

    metadata_file = Path("example_data") / "realistic_metadata.csv"
    if not metadata_file.exists():
        print("Metadata file not found. Run main example first.")
        return

    # Test different similarity thresholds
    thresholds = [0.10, 0.20, 0.30, 0.40]

    for threshold in thresholds:
        print(f"\nTesting similarity threshold: {threshold}")

        cmd = [
            sys.executable,
            "metagrouper.py",
            "example_data",
            "--output",
            f"threshold_test_{threshold:.2f}".replace(".", "p"),
            "--kmer-size",
            "15",
            "--max-reads",
            "300",
            "--metadata",
            str(metadata_file),
            "--assembly-tools",
            "megahit",
            "--permutations",
            "99",
            "--similarity-threshold",
            str(threshold),
        ]

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Quick analysis of results
            strategy_file = (
                Path(f"threshold_test_{threshold:.2f}".replace(".", "p"))
                / "assembly_recommendations"
                / "assembly_strategy.md"
            )
            if strategy_file.exists():
                with open(strategy_file, "r") as f:
                    content = f.read()
                    if "Overall Strategy:" in content:
                        strategy_line = [
                            line
                            for line in content.split("\n")
                            if "Overall Strategy:" in line
                        ][0]
                        print(f"  Threshold {threshold}: {strategy_line}")

        except subprocess.CalledProcessError as e:
            print(f"  Threshold {threshold} failed: {e}")


def analyze_recommendations():
    """Analyze the generated assembly recommendations."""
    print("\n" + "=" * 60)
    print("Assembly Recommendations Analysis")
    print("=" * 60)

    # Look for comprehensive analysis results
    results_dir = Path("comprehensive_analysis")
    if not results_dir.exists():
        print("No comprehensive analysis results found. Run main example first.")
        return

    # Read and display key recommendations
    strategy_file = results_dir / "assembly_recommendations" / "assembly_strategy.md"
    if strategy_file.exists():
        print("\n📋 ASSEMBLY STRATEGY SUMMARY:")
        print("=" * 40)
        with open(strategy_file, "r") as f:
            content = f.read()

            # Extract key sections
            lines = content.split("\n")
            in_summary = False
            in_groups = False

            for line in lines:
                if "## Overall Strategy:" in line:
                    print(line)
                elif "**Confidence Score:**" in line:
                    print(line)
                elif "**Rationale:**" in line:
                    print(line)
                elif "## Strategy Summary" in line:
                    in_summary = True
                    print(f"\n{line}")
                elif "## Recommended Assembly Groups" in line:
                    in_summary = False
                    in_groups = True
                    print(f"\n{line}")
                elif "## Assembly Commands" in line:
                    in_groups = False
                    print(f"\n{line}")
                    break
                elif in_summary or in_groups:
                    if line.strip():
                        print(line)

    # Show generated scripts
    scripts_dir = results_dir / "assembly_recommendations"
    if scripts_dir.exists():
        scripts = list(scripts_dir.glob("run_*.sh"))
        if scripts:
            print(f"\n🔧 GENERATED ASSEMBLY SCRIPTS:")
            print("=" * 40)
            for script in scripts:
                print(f"  - {script.name}")

                # Show first few commands
                with open(script, "r") as f:
                    lines = f.readlines()
                    command_lines = [
                        line.strip()
                        for line in lines
                        if line.strip() and not line.startswith("#")
                    ]
                    if command_lines:
                        print(f"    First command: {command_lines[0][:80]}...")


if __name__ == "__main__":
    print("MetaGrouper Phase 3 (Assembly Recommendations) Example")
    print("=" * 60)

    # Run comprehensive analysis
    run_comprehensive_example()

    # Demonstrate different configurations
    demonstrate_assembly_tools()

    # Test threshold sensitivity
    demonstrate_threshold_sensitivity()

    # Analyze final recommendations
    analyze_recommendations()

    print("\n" + "=" * 60)
    print("Phase 3 examples complete!")
    print("=" * 60)
    print("\nNext steps:")
    print("1. Review assembly_strategy.md for detailed recommendations")
    print("2. Examine assembly_strategy_overview.png for visual summary")
    print("3. Use the generated shell scripts to run assemblies")
    print("4. Compare results from different threshold settings")
    print("5. Validate recommendations with your biological knowledge")
