#!/usr/bin/env python3
"""
Test script to verify smart metadata variable filtering functionality.
"""

import sys
from pathlib import Path
import pandas as pd

# Add packages to path
sys.path.insert(0, str(Path(__file__).parent / "phases"))

from metadata_analyzer import filter_metadata_variables, is_biological_id

def test_biological_id_detection():
    """Test the biological ID detection function."""
    print("🧪 Testing biological ID detection...")
    
    test_cases = [
        # (column_name, expected_result)
        ("PID", True),
        ("GEMM", True),
        ("Patient_ID", True),
        ("Study_ID", True),
        ("Subject_ID", True),
        ("Sample_ID", False),
        ("Collaborator.ID", False),
        ("Golay.Barcode", False),
        ("Plate", False),
        ("Well", False),
        ("Unnamed: 0", False),
    ]
    
    passed = 0
    for col_name, expected in test_cases:
        result = is_biological_id(col_name)
        status = "✅" if result == expected else "❌"
        print(f"  {status} '{col_name}' -> {result} (expected: {expected})")
        if result == expected:
            passed += 1
    
    print(f"✅ Biological ID detection: {passed}/{len(test_cases)} tests passed")
    return passed == len(test_cases)

def test_filtering_with_real_metadata():
    """Test filtering with the actual metadata file."""
    print("\n🔍 Testing smart filtering with real metadata...")
    
    metadata_file = "./metadata_with_samples_final.csv"
    if not Path(metadata_file).exists():
        print(f"❌ Metadata file not found: {metadata_file}")
        return False
    
    # Load metadata
    metadata = pd.read_csv(metadata_file)
    print(f"📊 Loaded metadata: {len(metadata)} rows, {len(metadata.columns)} columns")
    
    # Test auto-filtering
    filtered_vars, exclusion_reasons = filter_metadata_variables(
        metadata, 
        auto_filter=True,
        exclude_variables=None
    )
    
    print(f"\n📋 Auto-filtering results:")
    print(f"  Original variables: {len(metadata.columns)}")
    print(f"  Filtered variables: {len(filtered_vars)}")
    print(f"  Excluded variables: {len(exclusion_reasons)}")
    
    print(f"\n✅ Included variables (biological):")
    for var in filtered_vars:
        values = metadata[var].dropna()
        unique_count = len(values.unique())
        print(f"  - {var}: {unique_count} unique values")
    
    print(f"\n❌ Excluded variables (technical/low quality):")
    exclusion_groups = {}
    for var, reason in exclusion_reasons.items():
        if reason not in exclusion_groups:
            exclusion_groups[reason] = []
        exclusion_groups[reason].append(var)
    
    for reason, vars_list in exclusion_groups.items():
        print(f"  {reason}:")
        for var in vars_list:
            print(f"    - {var}")
    
    # Expected biological variables that should be included
    expected_biological = [
        'celiacs_group', 'case_control', 'Sex', 'Delivery Mode', 
        'HLA', 'Dx Status', 'month', 'Country', 'GEMM', 'PID'
    ]
    
    biological_found = sum(1 for var in expected_biological if var in filtered_vars)
    print(f"\n🎯 Expected biological variables found: {biological_found}/{len(expected_biological)}")
    
    # Expected technical variables that should be excluded
    expected_technical = [
        'Unnamed: 0', 'Plate', 'Well', 'Sample_ID', 'Golay.Barcode', 'Collaborator.ID'
    ]
    
    technical_excluded = sum(1 for var in expected_technical if var in exclusion_reasons)
    print(f"🗑️  Expected technical variables excluded: {technical_excluded}/{len(expected_technical)}")
    
    success = biological_found >= 8 and technical_excluded >= 5
    if success:
        print("✅ Smart filtering working correctly!")
    else:
        print("❌ Smart filtering needs adjustment")
    
    return success

def test_manual_filtering():
    """Test manual variable specification."""
    print("\n🎯 Testing manual variable specification...")
    
    metadata_file = "./metadata_with_samples_final.csv"
    if not Path(metadata_file).exists():
        print(f"❌ Metadata file not found: {metadata_file}")
        return False
    
    metadata = pd.read_csv(metadata_file)
    
    # Test manual inclusion
    manual_vars = ['celiacs_group', 'case_control', 'Sex', 'GEMM', 'month']
    filtered_vars, exclusion_reasons = filter_metadata_variables(
        metadata,
        auto_filter=False,
        include_variables=manual_vars
    )
    
    print(f"Manual specification:")
    print(f"  Requested: {manual_vars}")
    print(f"  Included: {filtered_vars}")
    print(f"  Match: {set(manual_vars) == set(filtered_vars)}")
    
    # Test manual exclusion
    exclude_vars = ['Plate', 'Well', 'Unnamed: 0']
    filtered_vars, exclusion_reasons = filter_metadata_variables(
        metadata,
        auto_filter=True,
        exclude_variables=exclude_vars
    )
    
    excluded_found = sum(1 for var in exclude_vars if var in exclusion_reasons)
    print(f"Manual exclusion:")
    print(f"  Requested exclusions: {exclude_vars}")
    print(f"  Actually excluded: {excluded_found}/{len(exclude_vars)}")
    
    success = excluded_found == len(exclude_vars)
    if success:
        print("✅ Manual filtering working correctly!")
    else:
        print("❌ Manual filtering needs adjustment")
    
    return success

def main():
    """Run all filtering tests."""
    print("🧬 Testing MetaGrouper Smart Variable Filtering")
    print("=" * 60)
    
    tests = [
        test_biological_id_detection,
        test_filtering_with_real_metadata,
        test_manual_filtering
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
        print("\n🎉 All tests passed! Smart filtering is ready for use.")
        print("\n🚀 Try these commands:")
        print("  # Auto-filter variables")
        print("  python metagrouper.py /path/to/fastq --metadata metadata_with_samples_final.csv --auto-filter-variables")
        print("  # Manual specification")
        print("  python metagrouper.py /path/to/fastq --metadata metadata_with_samples_final.csv --variables celiacs_group case_control Sex")
        print("  # Auto-filter with exclusions")
        print("  python metagrouper.py /path/to/fastq --metadata metadata_with_samples_final.csv --auto-filter-variables --exclude-variables Plate Well")
        return True
    else:
        print("⚠️  Some tests failed. The filtering system may need adjustments.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)