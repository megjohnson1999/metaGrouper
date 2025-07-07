#!/usr/bin/env python3
"""
MetaGrouper Preprocessing Helper

Automated preprocessing pipeline for metagenomic samples before k-mer analysis.
Handles quality trimming, adapter removal, host contamination removal, and formatting.
"""

import argparse
import subprocess
import sys
import os
import logging
from pathlib import Path
from typing import List, Tuple, Dict, Optional
import shutil
import json
from datetime import datetime
import glob


class PreprocessingPipeline:
    """Automated preprocessing pipeline for metagenomic samples."""
    
    def __init__(self, output_dir: str = "preprocessed_data", verbose: bool = False):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.verbose = verbose
        self.setup_logging()
        
        # Tool availability
        self.tools_available = self.check_dependencies()
        
    def setup_logging(self):
        """Set up logging configuration."""
        level = logging.DEBUG if self.verbose else logging.INFO
        logging.basicConfig(
            level=level,
            format="%(asctime)s - %(levelname)s - %(message)s",
            handlers=[
                logging.StreamHandler(sys.stdout),
                logging.FileHandler(self.output_dir / "preprocessing.log")
            ]
        )
        
    def check_dependencies(self) -> Dict[str, bool]:
        """Check availability of preprocessing tools."""
        tools = {
            "fastp": shutil.which("fastp") is not None,
            "bowtie2": shutil.which("bowtie2") is not None,
            "seqkit": shutil.which("seqkit") is not None,
        }
        
        logging.info("Tool availability check:")
        for tool, available in tools.items():
            status = "✅" if available else "❌"
            logging.info(f"  {tool}: {status}")
            
        return tools
    
    def find_fastq_files(self, input_dir: str) -> List[Tuple[str, str, Optional[str]]]:
        """
        Find FASTQ files and determine pairing.
        Returns list of (sample_name, read1_path, read2_path_or_None)
        """
        input_path = Path(input_dir)
        
        # Common FASTQ extensions
        extensions = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
        
        files = []
        for ext in extensions:
            files.extend(list(input_path.glob(ext)))
            files.extend(list(input_path.glob(f"**/{ext}")))
        
        # Group by sample name
        samples = {}
        for file_path in files:
            # Extract sample name (remove common suffixes)
            name = file_path.stem
            if name.endswith(('.fastq', '.fq')):
                name = name.rsplit('.', 1)[0]
            
            # Check for paired-end patterns
            if any(suffix in name for suffix in ['_1', '_R1', '.1']):
                sample_name = name.replace('_1', '').replace('_R1', '').replace('.1', '')
                if sample_name not in samples:
                    samples[sample_name] = [None, None]
                samples[sample_name][0] = str(file_path)
            elif any(suffix in name for suffix in ['_2', '_R2', '.2']):
                sample_name = name.replace('_2', '').replace('_R2', '').replace('.2', '')
                if sample_name not in samples:
                    samples[sample_name] = [None, None]
                samples[sample_name][1] = str(file_path)
            else:
                # Single-end or unpaired
                samples[name] = [str(file_path), None]
        
        # Convert to list format
        sample_list = []
        for sample_name, (r1, r2) in samples.items():
            if r1 is None and r2 is not None:
                # Only R2 found, treat as single-end
                sample_list.append((sample_name, r2, None))
            else:
                sample_list.append((sample_name, r1, r2))
        
        logging.info(f"Found {len(sample_list)} samples:")
        for sample_name, r1, r2 in sample_list:
            if r2:
                logging.info(f"  {sample_name}: PAIRED ({Path(r1).name}, {Path(r2).name})")
            else:
                logging.info(f"  {sample_name}: SINGLE ({Path(r1).name})")
        
        return sample_list
    
    def run_fastp(self, sample_name: str, read1: str, read2: Optional[str] = None,
                  min_length: int = 50, threads: int = 4) -> Tuple[str, Optional[str]]:
        """Run fastp for quality trimming and adapter removal."""
        
        if not self.tools_available["fastp"]:
            raise RuntimeError("fastp not found. Install with: conda install -c bioconda fastp")
        
        output_r1 = self.output_dir / f"{sample_name}_trimmed_1.fastq"
        
        if read2:
            # Paired-end
            output_r2 = self.output_dir / f"{sample_name}_trimmed_2.fastq"
            html_report = self.output_dir / f"{sample_name}_fastp.html"
            json_report = self.output_dir / f"{sample_name}_fastp.json"
            
            cmd = [
                "fastp",
                "-i", read1, "-I", read2,
                "-o", str(output_r1), "-O", str(output_r2),
                "--thread", str(threads),
                "--length_required", str(min_length),
                "--detect_adapter_for_pe",
                "--dedup",
                "--trim_poly_g",
                "--html", str(html_report),
                "--json", str(json_report)
            ]
            
            output_files = (str(output_r1), str(output_r2))
        else:
            # Single-end
            html_report = self.output_dir / f"{sample_name}_fastp.html"
            json_report = self.output_dir / f"{sample_name}_fastp.json"
            
            cmd = [
                "fastp",
                "-i", read1,
                "-o", str(output_r1),
                "--thread", str(threads),
                "--length_required", str(min_length),
                "--dedup",
                "--trim_poly_g",
                "--html", str(html_report),
                "--json", str(json_report)
            ]
            
            output_files = (str(output_r1), None)
        
        logging.info(f"Running fastp for {sample_name}...")
        if self.verbose:
            logging.debug(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if self.verbose:
                logging.debug(f"fastp stdout: {result.stdout}")
            logging.info(f"✅ fastp completed for {sample_name}")
            return output_files
        except subprocess.CalledProcessError as e:
            logging.error(f"❌ fastp failed for {sample_name}: {e.stderr}")
            raise
    
    def remove_host_contamination(self, sample_name: str, read1: str, read2: Optional[str] = None,
                                host_index: str = "human_genome", threads: int = 4) -> Tuple[str, Optional[str]]:
        """Remove host contamination using bowtie2."""
        
        if not self.tools_available["bowtie2"]:
            logging.warning("bowtie2 not found. Skipping host removal.")
            return read1, read2
        
        if not Path(f"{host_index}.1.bt2").exists():
            logging.warning(f"Host index {host_index} not found. Skipping host removal.")
            logging.info("To create human genome index:")
            logging.info("  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz")
            logging.info("  bowtie2-build GCF_000001405.39_GRCh38.p13_genomic.fna.gz human_genome")
            return read1, read2
        
        output_r1 = self.output_dir / f"{sample_name}_clean_1.fastq"
        
        if read2:
            # Paired-end
            output_r2 = self.output_dir / f"{sample_name}_clean_2.fastq"
            
            cmd = [
                "bowtie2",
                "-x", host_index,
                "-1", read1, "-2", read2,
                "--un-conc", str(self.output_dir / f"{sample_name}_clean_%.fastq"),
                "--threads", str(threads),
                "-S", "/dev/null"
            ]
            
            output_files = (str(output_r1), str(output_r2))
        else:
            # Single-end
            cmd = [
                "bowtie2",
                "-x", host_index,
                "-U", read1,
                "--un", str(output_r1),
                "--threads", str(threads),
                "-S", "/dev/null"
            ]
            
            output_files = (str(output_r1), None)
        
        logging.info(f"Removing host contamination for {sample_name}...")
        if self.verbose:
            logging.debug(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if self.verbose:
                logging.debug(f"bowtie2 stderr: {result.stderr}")
            logging.info(f"✅ Host removal completed for {sample_name}")
            return output_files
        except subprocess.CalledProcessError as e:
            logging.error(f"❌ Host removal failed for {sample_name}: {e.stderr}")
            raise
    
    def merge_paired_reads(self, sample_name: str, read1: str, read2: str) -> str:
        """Merge paired-end reads for MetaGrouper input."""
        output_file = self.output_dir / f"{sample_name}.fastq"
        
        logging.info(f"Merging paired reads for {sample_name}...")
        
        try:
            with open(output_file, 'w') as outf:
                # Copy read1
                with open(read1, 'r') as inf:
                    outf.write(inf.read())
                # Copy read2
                with open(read2, 'r') as inf:
                    outf.write(inf.read())
            
            logging.info(f"✅ Reads merged for {sample_name}")
            return str(output_file)
        except Exception as e:
            logging.error(f"❌ Read merging failed for {sample_name}: {e}")
            raise
    
    def process_sample(self, sample_name: str, read1: str, read2: Optional[str] = None,
                      host_index: Optional[str] = None, min_length: int = 50, threads: int = 4) -> str:
        """Process a single sample through the complete pipeline."""
        
        logging.info(f"🔄 Processing {sample_name}...")
        
        # Step 1: Quality trimming and adapter removal
        trimmed_r1, trimmed_r2 = self.run_fastp(sample_name, read1, read2, min_length, threads)
        
        # Step 2: Host contamination removal (if index provided)
        if host_index:
            clean_r1, clean_r2 = self.remove_host_contamination(
                sample_name, trimmed_r1, trimmed_r2, host_index, threads
            )
        else:
            clean_r1, clean_r2 = trimmed_r1, trimmed_r2
        
        # Step 3: Prepare final file for MetaGrouper
        if clean_r2:
            # Merge paired reads
            final_file = self.merge_paired_reads(sample_name, clean_r1, clean_r2)
        else:
            # Single-end, just rename
            final_file = self.output_dir / f"{sample_name}.fastq"
            shutil.copy2(clean_r1, final_file)
        
        logging.info(f"✅ {sample_name} processing complete: {final_file}")
        return str(final_file)
    
    def process_directory(self, input_dir: str, host_index: Optional[str] = None,
                         min_length: int = 50, threads: int = 4) -> Dict[str, str]:
        """Process all samples in a directory."""
        
        samples = self.find_fastq_files(input_dir)
        
        if not samples:
            raise ValueError(f"No FASTQ files found in {input_dir}")
        
        processed_files = {}
        
        logging.info(f"🚀 Starting preprocessing pipeline for {len(samples)} samples...")
        
        for i, (sample_name, read1, read2) in enumerate(samples, 1):
            logging.info(f"📊 Progress: {i}/{len(samples)}")
            
            try:
                final_file = self.process_sample(
                    sample_name, read1, read2, host_index, min_length, threads
                )
                processed_files[sample_name] = final_file
            except Exception as e:
                logging.error(f"❌ Failed to process {sample_name}: {e}")
                continue
        
        # Save processing summary
        summary = {
            "timestamp": datetime.now().isoformat(),
            "input_directory": str(input_dir),
            "output_directory": str(self.output_dir),
            "host_index": host_index,
            "min_length": min_length,
            "threads": threads,
            "processed_samples": processed_files,
            "tools_used": self.tools_available
        }
        
        with open(self.output_dir / "preprocessing_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
        
        logging.info(f"🎉 Preprocessing complete! {len(processed_files)}/{len(samples)} samples processed")
        logging.info(f"📁 Output directory: {self.output_dir}")
        logging.info(f"📋 Summary saved: {self.output_dir}/preprocessing_summary.json")
        
        return processed_files


def main():
    """Main function for preprocessing CLI."""
    parser = argparse.ArgumentParser(
        description="MetaGrouper Preprocessing Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic preprocessing
  python preprocess.py data/raw_fastq/ -o clean_data/

  # With host removal (human samples)
  python preprocess.py data/raw_fastq/ -o clean_data/ --host-index human_genome

  # Custom parameters
  python preprocess.py data/raw_fastq/ -o clean_data/ --min-length 75 --threads 8

  # Quick test mode
  python preprocess.py data/raw_fastq/ -o clean_data/ --min-length 30 --quick

Note: Install dependencies with:
  conda install -c bioconda fastp bowtie2 seqkit
        """
    )
    
    parser.add_argument("input_dir", help="Directory containing FASTQ files")
    parser.add_argument("-o", "--output", default="preprocessed_data", 
                       help="Output directory (default: preprocessed_data)")
    parser.add_argument("--host-index", help="Host genome index for contamination removal")
    parser.add_argument("--min-length", type=int, default=50,
                       help="Minimum read length after trimming (default: 50)")
    parser.add_argument("--threads", type=int, default=4,
                       help="Number of threads to use (default: 4)")
    parser.add_argument("--quick", action="store_true",
                       help="Quick mode: skip host removal, use min-length=30")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Verbose logging")
    
    args = parser.parse_args()
    
    # Quick mode adjustments
    if args.quick:
        args.min_length = 30
        args.host_index = None
        print("🚀 Quick mode enabled: min-length=30, no host removal")
    
    # Validate input directory
    if not Path(args.input_dir).exists():
        print(f"❌ Error: Input directory does not exist: {args.input_dir}")
        sys.exit(1)
    
    # Initialize pipeline
    pipeline = PreprocessingPipeline(output_dir=args.output, verbose=args.verbose)
    
    try:
        # Process samples
        processed_files = pipeline.process_directory(
            args.input_dir,
            host_index=args.host_index,
            min_length=args.min_length,
            threads=args.threads
        )
        
        # Print summary
        print(f"\n🎉 Preprocessing completed successfully!")
        print(f"📁 Processed files in: {args.output}")
        print(f"📊 Ready for MetaGrouper analysis:")
        print(f"   python metagrouper.py {args.output} -o metagrouper_results")
        
        if len(processed_files) < 10:
            print("\n📋 Processed samples:")
            for sample, file_path in processed_files.items():
                print(f"   {sample}: {Path(file_path).name}")
        
    except Exception as e:
        print(f"❌ Preprocessing failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()