#!/usr/bin/env python3
"""
MetaGrouper: K-mer-based Analysis for Optimal Metagenomic Assembly Grouping

This is a compatibility wrapper that imports from the new modular architecture.
For the full-featured version, use metagrouper_package/main.py
"""

import sys
import warnings
from pathlib import Path

# Add the new package to the path
package_path = Path(__file__).parent / "metagrouper_package"
if package_path.exists():
    sys.path.insert(0, str(package_path))
    
    # Import from the new modular structure
    try:
        from metagrouper import *
        from metagrouper_package.main import main
        
        print("🔄 Using new modular MetaGrouper architecture")
        print("   For the latest features, use: python metagrouper_package/main.py")
        print()
        
    except ImportError as e:
        print(f"⚠️  Failed to import modular version: {e}")
        print("   Falling back to legacy implementation...")
        
        # Fall back to the legacy implementation
        exec(open("metagrouper_legacy.py").read())
        sys.exit(0)
else:
    print("⚠️  Modular package not found, using legacy implementation")
    # Fall back to the legacy implementation
    exec(open("metagrouper_legacy.py").read())
    sys.exit(0)


if __name__ == "__main__":
    main()