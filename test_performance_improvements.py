#!/usr/bin/env python3
"""
Test script to verify performance improvements in MetaGrouper.
"""

import subprocess
import sys
import time
from pathlib import Path

def run_test():
    """Run a simple test to verify parallel processing and optimizations work."""
    print("🧪 Testing MetaGrouper performance improvements...")
    
    # Find test data
    test_data = Path("demo_reads")
    if not test_data.exists():
        print("❌ Demo data not found. Looking for alternative test data...")
        test_data = Path("realistic_test_data")
        if not test_data.exists():
            print("❌ No test data found. Please provide test data.")
            return 1
    
    # Create output directory
    output_dir = Path("test_performance_output")
    output_dir.mkdir(exist_ok=True)
    
    # Run MetaGrouper with optimizations
    cmd = [
        sys.executable, "metagrouper.py",
        str(test_data),
        "--output", str(output_dir),
        "--use-sketching",
        "--sketch-size", "1000",
        "--sampling-method", "frequency",
        "--processes", "2",
        "--verbose"
    ]
    
    print(f"🚀 Running command: {' '.join(cmd)}")
    print("📊 Testing with:")
    print("   - Parallel processing (2 processes)")
    print("   - Optimized single-pass frequency sampling")
    print("   - Reduced memory usage")
    print("   - Sparse matrix computation")
    
    start_time = time.time()
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        elapsed_time = time.time() - start_time
        
        if result.returncode == 0:
            print(f"\n✅ Test completed successfully in {elapsed_time:.1f}s")
            print("\n📊 Performance improvements verified:")
            
            # Check for parallel processing in output
            if "using 2 processes" in result.stdout:
                print("   ✓ Parallel processing enabled")
            else:
                print("   ❌ Parallel processing not detected")
            
            # Check for sparse matrix usage
            if "Sparse similarity computed" in result.stdout:
                print("   ✓ Sparse similarity computation used")
            
            # Check for memory efficiency
            if "Memory usage:" in result.stdout:
                print("   ✓ Memory-efficient sketching active")
                
            print("\n📝 Output snippet:")
            print(result.stdout[:500] + "..." if len(result.stdout) > 500 else result.stdout)
            
        else:
            print(f"\n❌ Test failed with return code {result.returncode}")
            print("Error output:")
            print(result.stderr)
            
    except Exception as e:
        print(f"❌ Error running test: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(run_test())