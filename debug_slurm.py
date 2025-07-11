#!/usr/bin/env python3
"""
Debug script to help identify slurm execution issues
"""

import sys
import os
import argparse
from pathlib import Path

def main():
    print("🔍 SLURM DEBUG SCRIPT")
    print("=" * 50)
    
    # Environment information
    print(f"🐍 Python executable: {sys.executable}")
    print(f"🐍 Python version: {sys.version}")
    print(f"📁 Current working directory: {os.getcwd()}")
    print(f"📁 Script location: {__file__}")
    print(f"📁 Script directory: {os.path.dirname(os.path.abspath(__file__))}")
    
    # Check if we're in the right directory
    expected_files = ['metagrouper.py', 'phases/', 'metagrouper_package/']
    print(f"\n📋 Checking for expected files:")
    for file in expected_files:
        exists = os.path.exists(file)
        print(f"   {file}: {'✅' if exists else '❌'}")
    
    # Python path
    print(f"\n🛤️  Python path:")
    for path in sys.path[:5]:  # Show first 5
        print(f"   {path}")
    
    # Environment variables
    print(f"\n🌍 Key environment variables:")
    for var in ['PATH', 'PYTHONPATH', 'SLURM_JOB_ID', 'SLURM_JOB_NAME']:
        value = os.environ.get(var, 'NOT SET')
        print(f"   {var}: {value}")
    
    # Command line arguments
    print(f"\n💬 Command line arguments:")
    print(f"   sys.argv: {sys.argv}")
    
    # Try importing critical modules
    print(f"\n📦 Module import test:")
    try:
        import argparse
        print(f"   argparse: ✅")
    except ImportError as e:
        print(f"   argparse: ❌ {e}")
    
    # Test argument parsing
    print(f"\n🧪 Argument parsing test:")
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help="Input directory")
    parser.add_argument("--phases", nargs="+", type=int, choices=[1, 2, 3, 4])
    parser.add_argument("--load-from", help="Load from directory")
    
    try:
        args = parser.parse_args()
        print(f"   ✅ Arguments parsed successfully")
        print(f"   📊 args.phases: {getattr(args, 'phases', 'None')}")
        print(f"   📁 args.load_from: {getattr(args, 'load_from', 'None')}")
        print(f"   📁 args.input_dir: {args.input_dir}")
    except Exception as e:
        print(f"   ❌ Argument parsing failed: {e}")

if __name__ == "__main__":
    main()