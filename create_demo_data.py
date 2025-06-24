#!/usr/bin/env python3
"""
Create demo data for MetaGrouper presentation.
"""

import os
import random
from pathlib import Path

def create_demo_fastq(filepath: str, sample_type: str, num_reads: int = 200):
    """Create realistic demo FASTQ files."""
    
    # Different sequence patterns for different sample types
    if sample_type == "gut":
        # Gut microbiome-like sequences (more diverse)
        base_sequences = [
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
            "TTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATT",
            "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
        ]
    elif sample_type == "oral":
        # Oral microbiome-like sequences (medium diversity)
        base_sequences = [
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        ]
    else:  # skin
        # Skin microbiome-like sequences (lower diversity, more similar)
        base_sequences = [
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        ]
    
    with open(filepath, 'w') as f:
        for i in range(num_reads):
            # Pick a base sequence and add some variation
            base_seq = random.choice(base_sequences)
            
            # Add some random mutations (5% of bases)
            seq_list = list(base_seq)
            for j in range(len(seq_list)):
                if random.random() < 0.05:  # 5% mutation rate
                    seq_list[j] = random.choice(['A', 'T', 'C', 'G'])
            
            seq = ''.join(seq_list)
            qual = 'I' * len(seq)
            
            f.write(f"@read_{i+1}\n{seq}\n+\n{qual}\n")

def create_metadata_file(output_path: str):
    """Create demo metadata file."""
    
    samples_data = [
        ("sample_gut_1", "gut", "patient_1", "baseline", "healthy"),
        ("sample_gut_2", "gut", "patient_2", "baseline", "healthy"), 
        ("sample_gut_3", "gut", "patient_3", "baseline", "IBD"),
        ("sample_oral_1", "oral", "patient_1", "baseline", "healthy"),
        ("sample_oral_2", "oral", "patient_2", "baseline", "healthy"),
        ("sample_skin_1", "skin", "patient_1", "baseline", "healthy"),
        ("sample_skin_2", "skin", "patient_2", "baseline", "healthy"),
        ("sample_skin_3", "skin", "patient_3", "baseline", "eczema"),
    ]
    
    with open(output_path, 'w') as f:
        f.write("sample_id,body_site,patient_id,timepoint,condition\n")
        for sample_id, body_site, patient_id, timepoint, condition in samples_data:
            f.write(f"{sample_id},{body_site},{patient_id},{timepoint},{condition}\n")

def main():
    """Create complete demo dataset."""
    
    # Create demo directory
    demo_dir = Path("presentation_demo")
    demo_dir.mkdir(exist_ok=True)
    
    print("🧬 Creating MetaGrouper presentation demo data...")
    
    # Sample types and their characteristics
    samples = [
        ("sample_gut_1", "gut"),
        ("sample_gut_2", "gut"), 
        ("sample_gut_3", "gut"),
        ("sample_oral_1", "oral"),
        ("sample_oral_2", "oral"),
        ("sample_skin_1", "skin"),
        ("sample_skin_2", "skin"),
        ("sample_skin_3", "skin"),
    ]
    
    # Create FASTQ files
    for sample_name, sample_type in samples:
        filepath = demo_dir / f"{sample_name}.fastq"
        create_demo_fastq(str(filepath), sample_type)
        print(f"✅ Created {filepath}")
    
    # Create metadata file
    metadata_path = demo_dir / "demo_metadata.csv"
    create_metadata_file(str(metadata_path))
    print(f"✅ Created {metadata_path}")
    
    # Create presentation script
    presentation_script = demo_dir / "run_demo.sh"
    with open(presentation_script, 'w') as f:
        f.write("""#!/bin/bash
# MetaGrouper Presentation Demo Script

echo "🧬 MetaGrouper Live Demo"
echo "======================="
echo
echo "📊 Dataset: 8 microbiome samples (gut, oral, skin)"
echo "⚡ Using sourmash for fast analysis"
echo

# Run MetaGrouper with sourmash
python3 ../metagrouper.py . \\
  --metadata demo_metadata.csv \\
  --output demo_results/ \\
  --use-sourmash \\
  --sourmash-scaled 100 \\
  --sourmash-save-sigs \\
  --comprehensive-report \\
  --processes 2 \\
  --verbose

echo
echo "🎉 Demo complete! Check demo_results/ for:"
echo "   • Interactive HTML report"  
echo "   • Assembly recommendations"
echo "   • Similarity visualizations"
""")
    
    presentation_script.chmod(0o755)
    print(f"✅ Created {presentation_script}")
    
    print(f"\n🎯 Demo ready in {demo_dir}/")
    print("To run presentation demo:")
    print(f"   cd {demo_dir}")
    print("   ./run_demo.sh")
    
    return True

if __name__ == "__main__":
    success = main()
    if success:
        print("\n🚀 PRESENTATION DEMO READY!")
    else:
        print("\n❌ Demo creation failed")