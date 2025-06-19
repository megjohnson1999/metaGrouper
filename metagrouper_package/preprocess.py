#!/usr/bin/env python3
"""
MetaGrouper Preprocessing Pipeline

Built-in preprocessing pipeline for raw FASTQ files including:
- Quality trimming and adapter removal
- Deduplication  
- Host removal (optional)
- Format conversion

Requires: fastp (for quality control) and bowtie2 (for host removal)
"""

import sys
import argparse
import logging
import subprocess
import shutil
from pathlib import Path
from typing import List, Tuple, Optional
import tempfile
import os


def setup_logging(verbose: bool = False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def check_dependencies() -> Tuple[bool, List[str]]:
    """Check if required tools are available."""
    required_tools = {
        'fastp': 'Quality trimming and adapter removal',
        'gzip': 'File compression/decompression'
    }
    
    optional_tools = {
        'bowtie2': 'Host contamination removal'
    }
    
    missing_required = []
    missing_optional = []
    
    for tool, description in required_tools.items():
        if not shutil.which(tool):
            missing_required.append(f"{tool} ({description})")
    
    for tool, description in optional_tools.items():
        if not shutil.which(tool):
            missing_optional.append(f"{tool} ({description})")
    
    if missing_required:
        logging.error(f"Missing required tools: {', '.join(missing_required)}")
        logging.error("Install with: conda install -c bioconda fastp")
        return False, missing_required
    
    if missing_optional:
        logging.warning(f"Missing optional tools: {', '.join(missing_optional)}")
        logging.warning("Install with: conda install -c bioconda bowtie2")
    
    return True, []


def find_fastq_files(input_dir: str) -> List[Tuple[str, str, Optional[str]]]:
    """
    Find FASTQ files in input directory.
    
    Returns:
        List of (R1_path, sample_name, R2_path) tuples
        R2_path is None for single-end data
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        raise ValueError(f"Input directory not found: {input_dir}")
    
    # Find all FASTQ files
    fastq_patterns = ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']
    all_files = []
    
    for pattern in fastq_patterns:
        all_files.extend(input_path.rglob(pattern))
    
    if not all_files:
        raise ValueError(f"No FASTQ files found in {input_dir}")
    
    # Group paired-end files
    samples = {}
    
    for file_path in all_files:
        file_name = file_path.name
        
        # Extract sample name and read number
        sample_name = None
        read_num = None
        
        # Common patterns for paired-end reads
        patterns = [
            ('_R1', '_R2'),
            ('_1', '_2'),
            ('.R1', '.R2'),
            ('.1', '.2')
        ]
        
        for r1_pattern, r2_pattern in patterns:
            if r1_pattern in file_name:
                sample_name = file_name.replace(r1_pattern, '').split('.')[0]
                read_num = 1
                break
            elif r2_pattern in file_name:
                sample_name = file_name.replace(r2_pattern, '').split('.')[0]
                read_num = 2
                break
        
        # If no paired-end pattern found, treat as single-end
        if sample_name is None:
            sample_name = file_name.split('.')[0]
            read_num = 1
        
        if sample_name not in samples:
            samples[sample_name] = {'R1': None, 'R2': None}
        
        if read_num == 1:
            samples[sample_name]['R1'] = str(file_path)
        elif read_num == 2:
            samples[sample_name]['R2'] = str(file_path)
    
    # Convert to list format
    sample_list = []
    for sample_name, reads in samples.items():
        if reads['R1']:
            sample_list.append((reads['R1'], sample_name, reads['R2']))
        else:
            logging.warning(f"No R1 file found for sample {sample_name}")
    
    logging.info(f"Found {len(sample_list)} samples")
    
    # Check for paired vs single-end
    paired_samples = sum(1 for _, _, r2 in sample_list if r2 is not None)
    single_samples = len(sample_list) - paired_samples
    
    if paired_samples > 0:
        logging.info(f"Paired-end samples: {paired_samples}")
    if single_samples > 0:
        logging.info(f"Single-end samples: {single_samples}")
    
    return sample_list


def run_fastp(r1_path: str, sample_name: str, output_dir: str,
              r2_path: Optional[str] = None, min_length: int = 50,
              threads: int = 1, quick: bool = False) -> Tuple[str, Optional[str]]:
    """
    Run fastp for quality trimming and adapter removal.
    
    Returns:
        Tuple of (output_R1_path, output_R2_path)
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Output file paths
    if r2_path:
        out_r1 = output_path / f"{sample_name}_R1_trimmed.fastq.gz"
        out_r2 = output_path / f"{sample_name}_R2_trimmed.fastq.gz"
    else:
        out_r1 = output_path / f"{sample_name}_trimmed.fastq.gz"
        out_r2 = None
    
    # Build fastp command
    cmd = [
        'fastp',
        '-i', r1_path,
        '-o', str(out_r1),
        '--thread', str(threads),
        '--length_required', str(min_length),
        '--compression', '6'
    ]
    
    if r2_path:
        cmd.extend(['-I', r2_path, '-O', str(out_r2)])
    
    if quick:
        cmd.extend(['--disable_quality_filtering', '--disable_length_filtering'])
    else:
        # Standard quality filtering
        cmd.extend([
            '--qualified_quality_phred', '20',
            '--unqualified_percent_limit', '30',
            '--n_base_limit', '5'
        ])
    
    # Adapter trimming (auto-detection)
    cmd.extend(['--detect_adapter_for_pe' if r2_path else '--disable_adapter_trimming'])
    
    # Deduplication
    if not quick:
        cmd.append('--dedup')
    
    # JSON and HTML reports
    json_report = output_path / f"{sample_name}_fastp.json"
    html_report = output_path / f"{sample_name}_fastp.html"
    cmd.extend(['--json', str(json_report), '--html', str(html_report)])
    
    # Run fastp
    try:
        logging.debug(f"Running fastp: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if result.stderr:
            logging.debug(f"fastp stderr: {result.stderr}")
        
        logging.info(f"✅ Quality trimming completed for {sample_name}")
        return str(out_r1), str(out_r2) if out_r2 else None
        
    except subprocess.CalledProcessError as e:
        logging.error(f"fastp failed for {sample_name}: {e}")
        logging.error(f"fastp stderr: {e.stderr}")
        raise


def run_host_removal(r1_path: str, sample_name: str, output_dir: str,
                    host_index: str, r2_path: Optional[str] = None,
                    threads: int = 1) -> Tuple[str, Optional[str]]:
    """
    Remove host contamination using bowtie2.
    
    Args:
        host_index: Path to bowtie2 index (without .bt2 extension)
    
    Returns:
        Tuple of (clean_R1_path, clean_R2_path)
    """
    if not shutil.which('bowtie2'):
        logging.warning("bowtie2 not available, skipping host removal")
        return r1_path, r2_path
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Output file paths
    if r2_path:
        clean_r1 = output_path / f"{sample_name}_R1_clean.fastq.gz"
        clean_r2 = output_path / f"{sample_name}_R2_clean.fastq.gz"
        unpaired = output_path / f"{sample_name}_unpaired.fastq.gz"
    else:
        clean_r1 = output_path / f"{sample_name}_clean.fastq.gz"
        clean_r2 = None
        unpaired = None
    
    # SAM output (temporary)
    sam_file = output_path / f"{sample_name}_aligned.sam"
    
    try:
        # Build bowtie2 command
        cmd = [
            'bowtie2',
            '-x', host_index,
            '--threads', str(threads),
            '--very-sensitive',
            '--no-unal',  # Don't output unaligned reads to SAM
            '-S', str(sam_file)
        ]
        
        if r2_path:
            cmd.extend(['-1', r1_path, '-2', r2_path])
        else:
            cmd.extend(['-U', r1_path])
        
        # Run bowtie2
        logging.debug(f"Running bowtie2: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logging.warning(f"bowtie2 had issues for {sample_name}: {result.stderr}")
        
        # Extract unaligned reads (these are the clean reads)
        if r2_path:
            # Paired-end: extract unaligned pairs
            extract_cmd = [
                'samtools', 'fastq',
                '-f', '12',  # Both reads unmapped
                '-F', '256',  # Not secondary alignment
                '-1', str(clean_r1),
                '-2', str(clean_r2),
                '-0', '/dev/null',  # Discard singletons
                str(sam_file)
            ]
        else:
            # Single-end: extract unaligned reads
            extract_cmd = [
                'samtools', 'fastq',
                '-f', '4',  # Unmapped
                '-F', '256',  # Not secondary alignment
                str(sam_file)
            ]
            
            # Compress output
            with open(clean_r1, 'wb') as f:
                extract_proc = subprocess.Popen(extract_cmd, stdout=subprocess.PIPE)
                gzip_proc = subprocess.Popen(['gzip'], stdin=extract_proc.stdout, stdout=f)
                extract_proc.stdout.close()
                gzip_proc.communicate()
        
        if shutil.which('samtools'):
            subprocess.run(extract_cmd, check=True)
        else:
            logging.warning("samtools not available, keeping all reads (no host removal)")
            return r1_path, r2_path
        
        # Clean up temporary SAM file
        sam_file.unlink(missing_ok=True)
        
        logging.info(f"✅ Host removal completed for {sample_name}")
        return str(clean_r1), str(clean_r2) if clean_r2 else None
        
    except Exception as e:
        logging.error(f"Host removal failed for {sample_name}: {e}")
        # Return original files if host removal fails
        return r1_path, r2_path


def preprocess_sample(r1_path: str, sample_name: str, output_dir: str,
                     r2_path: Optional[str] = None, min_length: int = 50,
                     threads: int = 1, quick: bool = False,
                     host_index: Optional[str] = None) -> Tuple[str, Optional[str]]:
    """
    Complete preprocessing pipeline for a single sample.
    
    Returns:
        Tuple of (final_R1_path, final_R2_path)
    """
    logging.info(f"📥 Processing {sample_name}")
    
    # Step 1: Quality trimming and adapter removal
    trimmed_r1, trimmed_r2 = run_fastp(
        r1_path, sample_name, output_dir, r2_path,
        min_length=min_length, threads=threads, quick=quick
    )
    
    # Step 2: Host removal (optional)
    if host_index:
        clean_r1, clean_r2 = run_host_removal(
            trimmed_r1, sample_name, output_dir, host_index,
            trimmed_r2, threads=threads
        )
    else:
        clean_r1, clean_r2 = trimmed_r1, trimmed_r2
    
    logging.info(f"✅ Completed preprocessing for {sample_name}")
    return clean_r1, clean_r2


def create_parser():
    """Create command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="MetaGrouper preprocessing pipeline for raw FASTQ files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic preprocessing
  python preprocess.py raw_data/ -o clean_data/
  
  # Fast mode for testing
  python preprocess.py raw_data/ -o clean_data/ --quick
  
  # With host removal
  python preprocess.py raw_data/ -o clean_data/ --host-index human_genome
  
  # Custom parameters
  python preprocess.py raw_data/ -o clean_data/ --min-length 100 --threads 8
        """
    )
    
    parser.add_argument("input_dir", help="Directory containing raw FASTQ files")
    parser.add_argument("-o", "--output", default="preprocessed_data",
                       help="Output directory (default: preprocessed_data)")
    
    # Processing parameters
    parser.add_argument("--min-length", type=int, default=50,
                       help="Minimum read length after trimming (default: 50)")
    parser.add_argument("--threads", type=int, default=1,
                       help="Number of CPU threads to use (default: 1)")
    parser.add_argument("--quick", action="store_true",
                       help="Fast mode (minimal quality filtering)")
    
    # Host removal
    parser.add_argument("--host-index", 
                       help="Bowtie2 index for host genome (for contamination removal)")
    
    # Other options
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Verbose logging")
    
    return parser


def main():
    """Main preprocessing workflow."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    print("🧬 MetaGrouper Preprocessing Pipeline")
    print("=" * 50)
    
    # Check dependencies
    deps_ok, missing = check_dependencies()
    if not deps_ok:
        print(f"❌ Missing required dependencies: {', '.join(missing)}")
        print("Install with: conda install -c bioconda fastp")
        return 1
    
    # Validate host index if provided
    if args.host_index:
        if not Path(f"{args.host_index}.1.bt2").exists():
            logging.error(f"Host index not found: {args.host_index}")
            logging.error("Build with: bowtie2-build genome.fasta index_name")
            return 1
        logging.info(f"Using host index: {args.host_index}")
    
    try:
        # Find input files
        samples = find_fastq_files(args.input_dir)
        
        print(f"📁 Found {len(samples)} samples")
        print(f"🔧 Processing parameters:")
        print(f"   • Minimum length: {args.min_length}")
        print(f"   • Threads: {args.threads}")
        print(f"   • Quick mode: {'Yes' if args.quick else 'No'}")
        print(f"   • Host removal: {'Yes' if args.host_index else 'No'}")
        print()
        
        # Create output directory
        output_path = Path(args.output)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Process each sample
        processed_samples = []
        failed_samples = []
        
        for i, (r1_path, sample_name, r2_path) in enumerate(samples):
            print(f"📥 Processing sample {i+1}/{len(samples)}: {sample_name}")
            
            try:
                final_r1, final_r2 = preprocess_sample(
                    r1_path, sample_name, args.output, r2_path,
                    min_length=args.min_length, threads=args.threads,
                    quick=args.quick, host_index=args.host_index
                )
                processed_samples.append((final_r1, sample_name, final_r2))
                
            except Exception as e:
                logging.error(f"Failed to process {sample_name}: {e}")
                failed_samples.append(sample_name)
        
        # Summary
        print(f"\n🎯 Preprocessing Complete!")
        print(f"✅ Successfully processed: {len(processed_samples)}/{len(samples)} samples")
        
        if failed_samples:
            print(f"❌ Failed samples: {', '.join(failed_samples)}")
        
        print(f"📁 Clean data saved to: {args.output}")
        print(f"🚀 Ready for MetaGrouper analysis!")
        
        # Save sample list for MetaGrouper
        sample_list_file = output_path / "sample_list.txt"
        with open(sample_list_file, 'w') as f:
            f.write("sample_name\tR1_path\tR2_path\n")
            for r1, name, r2 in processed_samples:
                f.write(f"{name}\t{r1}\t{r2 or 'NA'}\n")
        
        print(f"📋 Sample list saved to: {sample_list_file}")
        
        return 0
        
    except Exception as e:
        logging.error(f"Preprocessing failed: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())