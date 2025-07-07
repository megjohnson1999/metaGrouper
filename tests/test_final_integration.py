#!/usr/bin/env python3
"""
Final integration test with realistic sequences.
"""

import sys
import tempfile
from pathlib import Path

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent / "metagrouper_package"))

def create_realistic_fastq(sample_name: str, seq_length: int = 150) -> str:
    """Create a realistic FASTQ file with longer sequences."""
    import random
    
    # Generate realistic DNA sequences
    bases = ['A', 'T', 'C', 'G']
    sequences = []
    
    for i in range(50):  # 50 reads per file
        seq = ''.join(random.choices(bases, k=seq_length))
        qual = 'I' * seq_length
        sequences.append(f"@read_{i+1}\n{seq}\n+\n{qual}")
    
    content = '\n'.join(sequences) + '\n'
    return content

def test_realistic_scenario():
    """Test with realistic MetaGrouper scenario."""
    print("🧪 Testing realistic scenario...")
    
    try:
        from metagrouper.sourmash_profiler import SourmashProfiler
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Create realistic test files
        single_file = tmpdir / "sample_single.fastq"
        pair_r1 = tmpdir / "sample_paired_R1.fastq"
        pair_r2 = tmpdir / "sample_paired_R2.fastq"
        
        with open(single_file, 'w') as f:
            f.write(create_realistic_fastq("single"))
        with open(pair_r1, 'w') as f:
            f.write(create_realistic_fastq("paired_r1"))
        with open(pair_r2, 'w') as f:
            f.write(create_realistic_fastq("paired_r2"))
        
        # MetaGrouper format exactly as it would be passed
        samples = [
            (str(single_file), "sample_single"),
            ([str(pair_r1), str(pair_r2)], "sample_paired"),
        ]
        
        try:
            # Test with sourmash profiler
            profiler = SourmashProfiler(
                k=21,
                scaled=100,  # Lower scaled for more hashes with small dataset
                processes=1,
                track_abundance=False
            )
            
            signatures = profiler.process_samples_parallel(samples)
            
            # Check results
            if len(signatures) != 2:
                print(f"❌ Expected 2 samples, got {len(signatures)}")
                return False
            
            for sample_name, sig in signatures.items():
                num_hashes = len(sig.minhash)
                print(f"✅ Sample {sample_name}: {num_hashes} hashes")
                
                if num_hashes == 0:
                    print(f"⚠️  Warning: No hashes for {sample_name} (sequences might be too short)")
            
            # Test similarity computation
            similarity_matrix = profiler.compute_similarity_matrix(signatures)
            print(f"✅ Similarity matrix shape: {similarity_matrix.shape}")
            
            # Test export format
            profiles, sample_names = profiler.export_to_metagrouper_format(
                signatures, similarity_matrix
            )
            print(f"✅ Exported profiles for: {sample_names}")
            
            return True
            
        except Exception as e:
            print(f"❌ Error: {e}")
            import traceback
            traceback.print_exc()
            return False

def main():
    """Run final integration test."""
    print("🚀 Final Integration Test for SourmashProfiler")
    print("=" * 50)
    
    success = test_realistic_scenario()
    
    print("\n" + "=" * 50)
    if success:
        print("🎉 Final integration test PASSED!")
        print("\n✅ SourmashProfiler is ready for MetaGrouper integration")
        print("   • Handles MetaGrouper's exact input format")
        print("   • Processes single-end and paired-end files")
        print("   • Computes similarity matrices")
        print("   • Exports to MetaGrouper format")
    else:
        print("❌ Final integration test FAILED!")
        print("   Please check the errors above before using in production.")
    
    return success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)