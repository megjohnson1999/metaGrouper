#!/usr/bin/env python3
"""
Minimal CI test script for MetaGrouper that tests only core functionality.
"""

import sys
import os
import tempfile
import shutil
from pathlib import Path


def test_modular_imports():
    """Test that modular imports work."""
    print("🔧 Testing modular imports...")
    
    try:
        # Add modular package to path
        sys.path.insert(0, 'metagrouper_package')
        
        # Test core imports
        from metagrouper.profiler import KmerProfiler
        from metagrouper.utils import find_fastq_files
        from metagrouper.analyzer import SimilarityAnalyzer
        
        print("   ✓ Core module imports work")
        
        # Test package-level imports
        from metagrouper import KmerProfiler as PkgKmerProfiler
        from metagrouper import find_fastq_files as pkg_find_files
        
        print("   ✓ Package-level imports work")
        return True
        
    except ImportError as e:
        print(f"   ❌ Import failed: {e}")
        return False


def test_basic_kmer_profiling():
    """Test basic k-mer profiling with minimal data."""
    print("🧪 Testing k-mer profiling...")
    
    try:
        # Import after setting path
        sys.path.insert(0, 'metagrouper_package')
        from metagrouper.profiler import KmerProfiler
        from metagrouper.utils import find_fastq_files
        
        # Create minimal test data
        temp_dir = Path(tempfile.mkdtemp())
        
        try:
            # Create test FASTQ
            fastq_file = temp_dir / 'test.fastq'
            with open(fastq_file, 'w') as f:
                f.write('@read1\nATCGATCGATCGATCGATCG\n+\n~~~~~~~~~~~~~~~~~~~~\n')
                f.write('@read2\nGCTAGCTAGCTAGCTAGCTA\n+\n~~~~~~~~~~~~~~~~~~~~\n')
            
            # Test file discovery
            files = find_fastq_files(str(temp_dir))
            if len(files) != 1:
                print(f"   ❌ Expected 1 file, got {len(files)}")
                return False
            print("   ✓ File discovery works")
            
            # Test k-mer profiling
            profiler = KmerProfiler(k=15, max_reads=3)
            profile = profiler.profile_sample(str(fastq_file), 'test')
            
            if len(profile) == 0:
                print("   ❌ Profile is empty")
                return False
            
            if abs(sum(profile.values()) - 1.0) > 1e-6:
                print("   ❌ Profile not normalized")
                return False
                
            print(f"   ✓ K-mer profiling works ({len(profile)} k-mers)")
            return True
            
        finally:
            shutil.rmtree(temp_dir)
            
    except Exception as e:
        print(f"   ❌ Test failed: {e}")
        return False


def test_compatibility_wrapper():
    """Test that the compatibility wrapper imports work."""
    print("🔄 Testing compatibility wrapper...")
    
    try:
        # Test that metagrouper.py can be imported
        old_path = sys.path[:]
        
        # Import the compatibility wrapper
        import metagrouper
        
        # Check that it has the expected attributes
        expected_attrs = ['KmerProfiler', 'SimilarityAnalyzer', 'find_fastq_files']
        
        for attr in expected_attrs:
            if not hasattr(metagrouper, attr):
                print(f"   ❌ Missing attribute: {attr}")
                return False
        
        print("   ✓ Compatibility wrapper works")
        return True
        
    except Exception as e:
        print(f"   ❌ Compatibility test failed: {e}")
        return False
    finally:
        sys.path = old_path


def main():
    """Run minimal CI tests."""
    print("🎯 MetaGrouper Minimal CI Test")
    print("=" * 35)
    
    tests = [
        ("Modular Imports", test_modular_imports),
        ("K-mer Profiling", test_basic_kmer_profiling),
        ("Compatibility", test_compatibility_wrapper),
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            print(f"\n📋 {test_name}")
            print("-" * len(test_name))
            
            if test_func():
                print(f"   ✅ {test_name} PASSED")
                passed += 1
            else:
                print(f"   ❌ {test_name} FAILED")
                failed += 1
                
        except Exception as e:
            print(f"   💥 {test_name} FAILED: {e}")
            failed += 1
    
    print(f"\n🎯 Results: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("🎉 All minimal tests passed!")
        return True
    else:
        print("💥 Some tests failed!")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)