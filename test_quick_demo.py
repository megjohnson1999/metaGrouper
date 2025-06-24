#!/usr/bin/env python3
"""
Quick test to ensure sourmash integration works for presentation demo.
"""

import sys
import tempfile
import logging
from pathlib import Path

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent / "metagrouper_package"))

def create_demo_fastq(filename: str, num_reads: int = 100) -> None:
    """Create a demo FASTQ file with realistic sequences."""
    import random
    
    bases = ['A', 'T', 'C', 'G']
    with open(filename, 'w') as f:
        for i in range(num_reads):
            # Create 150bp reads
            seq = ''.join(random.choices(bases, k=150))
            qual = 'I' * 150
            f.write(f"@read_{i+1}\n{seq}\n+\n{qual}\n")

def test_presentation_demo():
    """Test the exact workflow for presentation."""
    print("🧪 Testing MetaGrouper for presentation demo...")
    
    # Setup logging to see what's happening
    logging.basicConfig(level=logging.INFO)
    
    try:
        from metagrouper.sourmash_profiler import SourmashProfiler
        print("✅ Successfully imported SourmashProfiler")
    except Exception as e:
        print(f"❌ Import failed: {e}")
        return False
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        print(f"📁 Working in: {tmpdir}")
        
        # Create demo files (simulating your dataset structure)
        demo_files = []
        for i in range(5):  # Small demo with 5 samples
            if i < 3:
                # Single-end samples
                filepath = tmpdir / f"sample_{i+1}.fastq"
                create_demo_fastq(str(filepath))
                demo_files.append((str(filepath), f"sample_{i+1}"))
            else:
                # Paired-end samples
                r1_path = tmpdir / f"sample_{i+1}_R1.fastq"
                r2_path = tmpdir / f"sample_{i+1}_R2.fastq"
                create_demo_fastq(str(r1_path))
                create_demo_fastq(str(r2_path))
                demo_files.append(([str(r1_path), str(r2_path)], f"sample_{i+1}"))
        
        print(f"📊 Created {len(demo_files)} demo samples")
        
        try:
            # Test sourmash profiler with presentation-ready parameters
            profiler = SourmashProfiler(
                k=21,
                scaled=100,  # Smaller for demo
                processes=1,
                track_abundance=False
            )
            print("✅ Created SourmashProfiler")
            
            # Process samples
            print("🔬 Processing samples...")
            signatures = profiler.process_samples_parallel(demo_files)
            print(f"✅ Processed {len(signatures)} samples successfully")
            
            # Test signature saving (the part that was failing)
            sig_path = tmpdir / "demo_signatures.sig"
            print("💾 Testing signature saving...")
            profiler.save_signatures(signatures, str(sig_path))
            print(f"✅ Saved signatures to {sig_path}")
            
            # Test similarity computation
            print("📊 Computing similarities...")
            similarity_matrix = profiler.compute_similarity_matrix(signatures)
            print(f"✅ Computed similarity matrix: {similarity_matrix.shape}")
            
            # Test export format
            print("📤 Testing export format...")
            profiles, sample_names = profiler.export_to_metagrouper_format(
                signatures, similarity_matrix
            )
            print(f"✅ Exported: {len(profiles)} profiles for {sample_names}")
            
            # Success summary
            print("\n" + "="*50)
            print("🎉 PRESENTATION DEMO READY!")
            print("✅ All core functions working")
            print("✅ Signature saving fixed")
            print("✅ Ready for live demo")
            print("="*50)
            
            return True
            
        except Exception as e:
            print(f"❌ Error during processing: {e}")
            import traceback
            traceback.print_exc()
            return False

def main():
    """Run the presentation readiness test."""
    print("🚀 MetaGrouper Presentation Readiness Test")
    print("=" * 50)
    
    # Check sourmash first
    try:
        import sourmash
        print("✅ Sourmash available")
    except ImportError:
        print("❌ Sourmash not available")
        return False
    
    success = test_presentation_demo()
    
    if success:
        print("\n🎯 READY FOR PRESENTATION!")
        print("💡 Demo command to show:")
        print("   python metagrouper.py demo_data/ --use-sourmash --comprehensive-report")
    else:
        print("\n❌ NOT READY - fix errors above first")
    
    return success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)