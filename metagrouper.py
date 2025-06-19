#!/usr/bin/env python3
"""
MetaGrouper: K-mer-based Analysis for Optimal Metagenomic Assembly Grouping

This is a wrapper that calls the main metagrouper.py from the package directory.
"""

import sys
import subprocess
from pathlib import Path

def main():
    """Wrapper to call the main metagrouper.py script."""
    # Path to the actual metagrouper.py
    package_path = Path(__file__).parent / "metagrouper_package" / "metagrouper.py"
    
    if not package_path.exists():
        print("❌ MetaGrouper package not found")
        print("Please run from the metagrouper_package directory or install properly")
        return 1
    
    # Call the actual metagrouper.py with all arguments
    try:
        result = subprocess.run([sys.executable, str(package_path)] + sys.argv[1:])
        return result.returncode
    except Exception as e:
        print(f"❌ Failed to run MetaGrouper: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())