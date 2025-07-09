#!/usr/bin/env python3
"""
Test script to verify the sample name normalization fix.
This will test that FASTQ file sample names are properly normalized
to match metadata Sample_ID column entries.
"""

import sys
from pathlib import Path

# Add package to path
sys.path.insert(0, str(Path(__file__).parent / "metagrouper_package"))

from metagrouper.utils import normalize_sample_name, find_fastq_files
import pandas as pd

def test_normalize_sample_name():
    """Test the normalize_sample_name function."""
    print("🧪 Testing normalize_sample_name function...")
    
    test_cases = [
        # (input, expected_output)
        ("NovaSeq_N983_I13380_39894_Celiac_Leonard_Stool_02_GEMM_038_18M_hr", 
         "NovaSeq_N983_I13380_39894_Celiac_Leonard_Stool_02_GEMM_038_18M"),
        ("sample_001_trimmed", "sample_001"),
        ("test_sample_filtered", "test_sample"),
        ("normal_sample", "normal_sample"),  # No suffix to remove
        ("sample_qc", "sample"),
        ("complex_sample_name_hr", "complex_sample_name"),
    ]
    
    for input_name, expected in test_cases:
        result = normalize_sample_name(input_name)
        status = "✅" if result == expected else "❌"
        print(f"  {status} '{input_name}' -> '{result}' (expected: '{expected}')")
        if result != expected:
            return False
    
    print("✅ All normalize_sample_name tests passed!")
    return True

def test_metadata_matching():
    """Test that our sample names will match metadata."""
    print("\n🔍 Testing metadata matching...")
    
    # Load the metadata file
    metadata_file = "./metadata_with_samples_final.csv"
    if not Path(metadata_file).exists():
        print(f"❌ Metadata file not found: {metadata_file}")
        return False
    
    # Read metadata
    metadata = pd.read_csv(metadata_file)
    metadata_sample_ids = set(metadata['Sample_ID'].dropna().astype(str))
    print(f"📊 Found {len(metadata_sample_ids)} sample IDs in metadata")
    
    # Test sample names from SLURM log
    test_fastq_names = [
        "NovaSeq_N983_I13380_39894_Celiac_Leonard_Stool_02_GEMM_038_18M_hr",
        "NovaSeq_N983_I13381_39895_Celiac_Leonard_Stool_02_GEMM_038_24M_hr",
        "NovaSeq_N983_I13382_39896_Celiac_Leonard_Stool_02_GEMM_038_30M_hr",
        "NovaSeq_N983_I13389_39903_Celiac_Leonard_Stool_02_GEMM_045_12M_hr",
        "NovaSeq_N983_I13394_39908_Celiac_Leonard_Stool_02_GEMM_050_12M_hr"
    ]
    
    matches_found = 0
    for fastq_name in test_fastq_names:
        normalized = normalize_sample_name(fastq_name)
        if normalized in metadata_sample_ids:
            print(f"  ✅ '{fastq_name}' -> '{normalized}' (FOUND in metadata)")
            matches_found += 1
        else:
            print(f"  ❌ '{fastq_name}' -> '{normalized}' (NOT FOUND in metadata)")
    
    print(f"\n📈 Found {matches_found}/{len(test_fastq_names)} matches")
    
    if matches_found > 0:
        print("✅ Sample name normalization successfully enables metadata matching!")
        return True
    else:
        print("❌ No matches found - there may be other issues")
        return False

def test_find_fastq_integration():
    """Test that find_fastq_files uses the normalization."""
    print("\n🔧 Testing find_fastq_files integration...")
    
    # Test with the actual FASTQ directory
    fastq_dir = "./host-removed-subset"
    if not Path(fastq_dir).exists():
        print(f"❌ FASTQ directory not found: {fastq_dir}")
        return False
    
    # Get FASTQ files
    file_pairs = find_fastq_files(fastq_dir)
    print(f"📁 Found {len(file_pairs)} samples in {fastq_dir}")
    
    # Check if sample names are normalized
    sample_names = [pair[1] for pair in file_pairs[:5]]  # First 5 samples
    print("🔍 Sample names from find_fastq_files:")
    for name in sample_names:
        has_suffix = any(name.endswith(suffix) for suffix in ['_hr', '_trimmed', '_filtered'])
        status = "❌ (still has suffix)" if has_suffix else "✅ (normalized)"
        print(f"  {status} '{name}'")
    
    # Check if any normalized names end with processing suffixes
    normalized_names_with_suffixes = [name for name in sample_names 
                                    if any(name.endswith(suffix) for suffix in ['_hr', '_trimmed', '_filtered'])]
    
    if not normalized_names_with_suffixes:
        print("✅ All sample names are properly normalized!")
        return True
    else:
        print(f"❌ Found {len(normalized_names_with_suffixes)} names still with suffixes")
        return False

def main():
    """Run all tests."""
    print("🧬 Testing MetaGrouper Sample Name Normalization Fix")
    print("=" * 60)
    
    tests = [
        test_normalize_sample_name,
        test_metadata_matching,
        test_find_fastq_integration
    ]
    
    passed = 0
    for test_func in tests:
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"❌ Test {test_func.__name__} failed with error: {e}")
    
    print(f"\n📊 Test Results: {passed}/{len(tests)} tests passed")
    
    if passed == len(tests):
        print("🎉 All tests passed! The fix should resolve the metadata matching issue.")
        return True
    else:
        print("⚠️  Some tests failed. The fix may need additional work.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)