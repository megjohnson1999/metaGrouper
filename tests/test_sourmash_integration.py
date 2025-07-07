#!/usr/bin/env python3
"""
Test script for sourmash integration with MetaGrouper.
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

def test_sourmash_profiler():
    """Test the SourmashProfiler class."""
    print("🧪 Testing SourmashProfiler...")
    
    try:
        from metagrouper.sourmash_profiler import SourmashProfiler
        print("✅ Successfully imported SourmashProfiler")
    except ImportError as e:
        print(f"❌ Failed to import SourmashProfiler: {e}")
        return False
    
    # Create test data
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Create simple test FASTQ files
        test_sequences = {
            'sample1': """@seq1
ATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIII
@seq2
GCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIII
""",
            'sample2': """@seq1
ATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIII
@seq3
TTTTTTTTTTTTTTTTTTTTTTTT
+
IIIIIIIIIIIIIIIIIIIIIIII
"""
        }
        
        # Write test files
        test_files = {}
        for sample_name, content in test_sequences.items():
            filepath = tmpdir / f"{sample_name}.fastq"
            create_test_fastq(content, str(filepath))
            test_files[sample_name] = str(filepath)
        
        try:
            # Initialize profiler
            profiler = SourmashProfiler(
                k=21,
                scaled=1000,
                processes=1,
                track_abundance=False
            )
            print("✅ Successfully created SourmashProfiler")
            
            # Process samples
            signatures = profiler.process_samples_parallel(test_files)
            print(f"✅ Successfully processed {len(signatures)} samples")
            
            # Compute similarity matrix
            similarity_matrix = profiler.compute_similarity_matrix(signatures)
            print(f"✅ Successfully computed similarity matrix: {similarity_matrix.shape}")
            
            # Export to MetaGrouper format
            profiles, sample_names = profiler.export_to_metagrouper_format(
                signatures, similarity_matrix
            )
            print(f"✅ Successfully exported to MetaGrouper format: {len(profiles)} profiles")
            
            return True
            
        except Exception as e:
            print(f"❌ Error during processing: {e}")
            import traceback
            traceback.print_exc()
            return False

def test_cli_integration():
    """Test CLI integration."""
    print("\n🧪 Testing CLI integration...")
    
    # Test that the new arguments are available
    try:
        sys.path.insert(0, str(Path(__file__).parent / "metagrouper_package"))
        from metagrouper import create_parser
        
        parser = create_parser()
        
        # Test parsing sourmash arguments
        test_args = [
            "/fake/path",
            "--use-sourmash",
            "--sourmash-scaled", "2000",
            "--sourmash-track-abundance",
            "--sourmash-save-sigs"
        ]
        
        args = parser.parse_args(test_args)
        
        assert args.use_sourmash == True
        assert args.sourmash_scaled == 2000
        assert args.sourmash_track_abundance == True
        assert args.sourmash_save_sigs == True
        
        print("✅ CLI integration test passed")
        return True
        
    except Exception as e:
        print(f"❌ CLI integration test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("🚀 Testing MetaGrouper + Sourmash Integration")
    print("=" * 50)
    
    success = True
    
    # Test sourmash availability
    try:
        import sourmash
        print(f"✅ Sourmash available: {sourmash.__version__}")
    except ImportError:
        print("❌ Sourmash not available - please install: conda install sourmash")
        return False
    
    # Test profiler
    if not test_sourmash_profiler():
        success = False
    
    # Test CLI
    if not test_cli_integration():
        success = False
    
    print("\n" + "=" * 50)
    if success:
        print("🎉 All tests passed! Sourmash integration is working.")
        print("\n💡 To use sourmash with MetaGrouper, add --use-sourmash to your command:")
        print("   python metagrouper.py /path/to/fastq --use-sourmash --sourmash-scaled 1000")
    else:
        print("❌ Some tests failed. Please fix the issues above.")
    
    return success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)