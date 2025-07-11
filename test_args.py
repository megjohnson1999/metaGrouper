#!/usr/bin/env python3
"""
Simple test script to debug argument parsing issues
"""

import argparse

def create_parser():
    """Create argument parser with phase selection options."""
    parser = argparse.ArgumentParser(description="Test argument parsing")
    
    parser.add_argument("input_dir", help="Input directory")
    parser.add_argument("-o", "--output", default="metagrouper_output", help="Output directory")
    
    # Phase selection arguments
    parser.add_argument("--phases", nargs="+", type=int, choices=[1, 2, 3, 4],
                       help="Run specific phases only (e.g., --phases 2 3)")
    parser.add_argument("--load-from", 
                       help="Load Phase 1 results from previous run directory")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    
    return parser

def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    print(f"🔍 TEST DEBUG: Immediately after parse_args()")
    print(f"🔍 TEST DEBUG: args.phases = {getattr(args, 'phases', 'NOT FOUND')}")
    print(f"🔍 TEST DEBUG: args.load_from = {getattr(args, 'load_from', 'NOT FOUND')}")
    print(f"🔍 TEST DEBUG: args.verbose = {getattr(args, 'verbose', 'NOT FOUND')}")
    print(f"🔍 TEST DEBUG: input_dir = {args.input_dir}")
    print(f"🔍 TEST DEBUG: output = {args.output}")
    print(f"🔍 TEST DEBUG: ALL ARGS: {vars(args)}")
    
    return 0

if __name__ == "__main__":
    main()