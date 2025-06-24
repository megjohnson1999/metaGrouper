#!/usr/bin/env python3
"""
Test sourmash profiler with real MetaGrouper data formats.
"""

import sys
import os
import tempfile
from pathlib import Path

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent / "metagrouper_package"))

def create_test_fastq(content: str, filename: str) -> str:
    """Create a test FASTQ file."""
    with open(filename, 'w') as f:
        f.write(content)
    return filename

def test_metagrouper_format():
    """Test with the exact format that MetaGrouper uses."""
    print("🧪 Testing with MetaGrouper format...")
    
    try:
        from metagrouper.sourmash_profiler import SourmashProfiler
        print("✅ Successfully imported SourmashProfiler")
    except ImportError as e:
        print(f"❌ Failed to import SourmashProfiler: {e}")
        return False
    
    # Create test data in temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Create test FASTQ files (making sure seq and quality lengths match)
        test_content = """@seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@seq2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
        
        # Create single-end files
        single_file1 = tmpdir / "sample1.fastq"
        single_file2 = tmpdir / "sample2.fastq"
        create_test_fastq(test_content, str(single_file1))
        create_test_fastq(test_content, str(single_file2))
        
        # Create paired-end files
        pair1_r1 = tmpdir / "sample3_R1.fastq"
        pair1_r2 = tmpdir / "sample3_R2.fastq"
        pair2_r1 = tmpdir / "sample4_R1.fastq"
        pair2_r2 = tmpdir / "sample4_R2.fastq"
        
        create_test_fastq(test_content, str(pair1_r1))
        create_test_fastq(test_content, str(pair1_r2))
        create_test_fastq(test_content, str(pair2_r1))
        create_test_fastq(test_content, str(pair2_r2))
        
        # Test MetaGrouper format: List[Tuple[Union[str, List[str]], str]]
        metagrouper_samples = [
            # Single-end samples: (str_path, sample_name)
            (str(single_file1), "sample1"),
            (str(single_file2), "sample2"),
            # Paired-end samples: ([path1, path2], sample_name)
            ([str(pair1_r1), str(pair1_r2)], "sample3"),
            ([str(pair2_r1), str(pair2_r2)], "sample4"),
        ]
        
        print(f"📊 Testing with {len(metagrouper_samples)} samples:")
        for i, sample in enumerate(metagrouper_samples):
            filepath, sample_name = sample
            if isinstance(filepath, list):
                print(f"  {i+1}. {sample_name}: paired-end {len(filepath)} files")
            else:
                print(f"  {i+1}. {sample_name}: single-end")
        
        try:
            # Initialize profiler
            profiler = SourmashProfiler(
                k=21,
                scaled=1000,
                processes=1,  # Use single process for testing
                track_abundance=False
            )
            print("✅ Successfully created SourmashProfiler")
            
            # Process samples with MetaGrouper format
            signatures = profiler.process_samples_parallel(metagrouper_samples)
            print(f"✅ Successfully processed {len(signatures)} samples")
            
            # Verify all samples were processed
            expected_samples = ["sample1", "sample2", "sample3", "sample4"]
            for sample_name in expected_samples:
                if sample_name in signatures:
                    sig = signatures[sample_name]
                    print(f"✅ Sample {sample_name}: {len(sig.minhash)} hashes")
                else:
                    print(f"❌ Missing sample: {sample_name}")
                    return False
            
            # Compute similarity matrix
            similarity_matrix = profiler.compute_similarity_matrix(signatures)
            print(f"✅ Successfully computed similarity matrix: {similarity_matrix.shape}")
            
            # Test export to MetaGrouper format
            profiles, sample_names = profiler.export_to_metagrouper_format(
                signatures, similarity_matrix
            )
            print(f"✅ Successfully exported: {len(profiles)} profiles, {len(sample_names)} sample names")
            
            return True
            
        except Exception as e:
            print(f"❌ Error during processing: {e}")
            import traceback
            traceback.print_exc()
            return False

def test_dict_format():
    """Test with dictionary format (backward compatibility)."""
    print("\n🧪 Testing dictionary format (backward compatibility)...")
    
    try:
        from metagrouper.sourmash_profiler import SourmashProfiler
    except ImportError as e:
        print(f"❌ Failed to import SourmashProfiler: {e}")
        return False
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Create test files (making sure seq and quality lengths match)
        test_content = """@seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
        
        file1 = tmpdir / "test1.fastq"
        file2 = tmpdir / "test2.fastq"
        create_test_fastq(test_content, str(file1))
        create_test_fastq(test_content, str(file2))
        
        # Dictionary format: {sample_name: filepath}
        dict_samples = {
            "test1": str(file1),
            "test2": str(file2),
        }
        
        try:
            profiler = SourmashProfiler(k=21, scaled=1000, processes=1)
            signatures = profiler.process_samples_parallel(dict_samples)
            
            if len(signatures) == 2 and "test1" in signatures and "test2" in signatures:
                print("✅ Dictionary format test passed")
                return True
            else:
                print(f"❌ Dictionary format test failed: got {list(signatures.keys())}")
                return False
                
        except Exception as e:
            print(f"❌ Dictionary format test error: {e}")
            return False

def main():
    """Run comprehensive tests."""
    print("🚀 Testing SourmashProfiler with Real MetaGrouper Formats")
    print("=" * 60)
    
    # Check sourmash availability
    try:
        import sourmash
        try:
            version = sourmash.__version__
        except AttributeError:
            version = "unknown"
        print(f"✅ Sourmash available: {version}")
    except ImportError:
        print("❌ Sourmash not available - please install: pip install sourmash")
        return False
    
    success = True
    
    # Test MetaGrouper format
    if not test_metagrouper_format():
        success = False
    
    # Test dictionary format
    if not test_dict_format():
        success = False
    
    print("\n" + "=" * 60)
    if success:
        print("🎉 All tests passed! SourmashProfiler handles all input formats correctly.")
        print("\n💡 The profiler now correctly handles:")
        print("   • Single-end files: (str_path, sample_name)")
        print("   • Paired-end files: ([path1, path2], sample_name)")
        print("   • Dictionary format: {sample_name: filepath}")
    else:
        print("❌ Some tests failed. Check the errors above.")
    
    return success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)