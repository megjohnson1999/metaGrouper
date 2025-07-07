#!/usr/bin/env python3
"""
Simple CI test script to verify MetaGrouper functionality in CI environment.
"""

import sys
import tempfile
import shutil
from pathlib import Path


def test_imports():
    """Test that all imports work correctly."""
    print("🔧 Testing imports...")
    
    try:
        # Test modular imports
        sys.path.insert(0, 'metagrouper_package')
        from metagrouper import (
            SourmashProfiler, 
            SimilarityAnalyzer, 
            Visualizer,
            find_fastq_files,
            setup_logging,
            save_results,
            MetaGrouperConfig
        )
        print("   ✓ All modular imports successful")
    except ImportError as e:
        print(f"   ❌ Modular import failed: {e}")
        return False
    
    # Test legacy imports still work
    try:
        import metagrouper as legacy_mg
        print("   ✓ Legacy metagrouper imports work")
    except ImportError as e:
        print(f"   ⚠️  Legacy import failed (this is OK in CI): {e}")
    
    return True


def test_basic_functionality():
    """Test basic k-mer profiling functionality."""
    print("🧪 Testing basic functionality...")
    
    # Import after path setup
    sys.path.insert(0, 'metagrouper_package')
    from metagrouper import SourmashProfiler, find_fastq_files
    
    # Create minimal test data
    temp_dir = Path(tempfile.mkdtemp())
    try:
        # Create test FASTQ file
        fastq_file = temp_dir / 'test.fastq'
        with open(fastq_file, 'w') as f:
            f.write('@read1\nATCGATCGATCGATCGATCG\n+\n~~~~~~~~~~~~~~~~~~~~\n')
            f.write('@read2\nGCTAGCTAGCTAGCTAGCTA\n+\n~~~~~~~~~~~~~~~~~~~~\n')
            f.write('@read3\nTTCCGGAATTCCGGAATTCC\n+\n~~~~~~~~~~~~~~~~~~~~\n')
        
        # Test file discovery
        files = find_fastq_files(str(temp_dir))
        assert len(files) == 1, f'Expected 1 file, got {len(files)}'
        print('   ✓ File discovery works')
        
        # Test k-mer profiling (use num_hashes for small test data)
        profiler = SourmashProfiler(k=15, processes=1, num_hashes=100, scaled=0)
        signature = profiler.sketch_sample(str(fastq_file), 'test')
        profile = {str(h): 1 for h in signature.minhash.hashes}
        
        assert len(profile) > 0, 'Profile should not be empty'
        print(f'   ✓ K-mer profiling works ({len(profile)} k-mers)')
        
        return True
        
    finally:
        # Cleanup
        shutil.rmtree(temp_dir)


def test_compatibility_wrapper():
    """Test that the compatibility wrapper works."""
    print("🔄 Testing compatibility wrapper...")
    
    # Test import through wrapper
    import metagrouper
    
    # Check that key classes are available
    assert hasattr(metagrouper, 'SourmashProfiler'), "SourmashProfiler not available"
    assert hasattr(metagrouper, 'SimilarityAnalyzer'), "SimilarityAnalyzer not available"
    assert hasattr(metagrouper, 'find_fastq_files'), "find_fastq_files not available"
    
    print("   ✓ Compatibility wrapper works")
    return True


def main():
    """Run all CI tests."""
    print("🎯 MetaGrouper CI Test Suite")
    print("=" * 40)
    
    tests = [
        ("Import Test", test_imports),
        ("Functionality Test", test_basic_functionality),
        ("Compatibility Test", test_compatibility_wrapper),
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            print(f"\n📋 {test_name}")
            print("-" * len(test_name))
            result = test_func()
            if result:
                print(f"   ✅ {test_name} PASSED")
                passed += 1
            else:
                print(f"   ❌ {test_name} FAILED")
                failed += 1
        except Exception as e:
            print(f"   💥 {test_name} FAILED with exception: {e}")
            failed += 1
    
    print(f"\n🎯 CI Test Results")
    print("=" * 40)
    print(f"   ✅ Passed: {passed}")
    print(f"   ❌ Failed: {failed}")
    
    if failed == 0:
        print("   🎉 All tests passed!")
        return True
    else:
        print("   💥 Some tests failed!")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)