"""
K-mer profiling functionality for metagenomic samples.

This module contains the KmerProfiler class and related functions for extracting
and analyzing k-mer profiles from FASTQ files with support for both single-end
and paired-end reads.
"""

import logging
import gzip
import multiprocessing as mp
from multiprocessing import Pool, Manager
import time
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Union
from pathlib import Path

try:
    from .fast_kmer import FastKmerExtractor
    FAST_KMER_AVAILABLE = True
except ImportError:
    FAST_KMER_AVAILABLE = False
    logging.warning("Fast k-mer extraction not available")


def process_sample_worker(args):
    """
    Worker function for multiprocessing sample processing.
    
    Args:
        args: Tuple containing (filepath, sample_name, k, max_reads, min_kmer_freq, memory_efficient, use_fast_extraction, progress_dict, lock)
    
    Returns:
        Tuple containing (sample_name, profile, error_message)
    """
    filepath, sample_name, k, max_reads, min_kmer_freq, memory_efficient, use_fast_extraction, progress_dict, lock = args
    
    try:
        # Create a temporary profiler for this worker
        profiler = KmerProfiler(k=k, max_reads=max_reads, min_kmer_freq=min_kmer_freq, 
                               use_fast_extraction=use_fast_extraction)
        profile = profiler.profile_sample(filepath, sample_name, memory_efficient=memory_efficient)
        
        # Update progress safely
        with lock:
            progress_dict['completed'] += 1
            progress_dict['current_sample'] = sample_name
        
        return (sample_name, profile, None)
    
    except Exception as e:
        # Update progress and record error
        with lock:
            progress_dict['completed'] += 1
            progress_dict['failed'] += 1
            progress_dict['errors'][sample_name] = str(e)
        
        return (sample_name, None, str(e))


class KmerProfiler:
    """Handle k-mer profiling of metagenomic samples."""

    def __init__(self, k: int = 21, max_reads: Optional[int] = None, min_kmer_freq: int = 1, 
                 use_fast_extraction: bool = True):
        self.k = k
        self.max_reads = max_reads
        self.min_kmer_freq = min_kmer_freq  # Filter out low-frequency k-mers
        self.profiles = {}
        self.sample_names = []
        
        # Initialize fast extraction if available and requested
        self.use_fast_extraction = use_fast_extraction and FAST_KMER_AVAILABLE and k <= 32
        if self.use_fast_extraction:
            self.fast_extractor = FastKmerExtractor(k)
            logging.debug(f"Using fast k-mer extraction for k={k}")
        else:
            if use_fast_extraction and not FAST_KMER_AVAILABLE:
                logging.warning("Fast k-mer extraction requested but not available, using legacy method")
            elif use_fast_extraction and k > 32:
                logging.warning(f"Fast k-mer extraction not supported for k={k} > 32, using legacy method")
            logging.debug(f"Using legacy k-mer extraction for k={k}")

    def _open_file(self, filepath: str):
        """Open file, handling gzip compression."""
        if filepath.endswith(".gz"):
            return gzip.open(filepath, "rt")
        return open(filepath, "r")

    def _extract_kmers(self, sequence: str) -> Dict[str, int]:
        """Extract k-mers from a sequence."""
        if self.use_fast_extraction:
            return self.fast_extractor.extract_kmers_with_strings(sequence)
        else:
            return self._extract_kmers_legacy(sequence)
    
    def _extract_kmers_legacy(self, sequence: str) -> Dict[str, int]:
        """Legacy k-mer extraction using string operations."""
        kmers = defaultdict(int)
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i : i + self.k]
            if "N" not in kmer:  # Skip k-mers with ambiguous bases
                # Use canonical k-mer (lexicographically smaller of forward/reverse)
                rev_comp = self._reverse_complement(kmer)
                canonical = min(kmer, rev_comp)
                kmers[canonical] += 1
        return kmers

    def _reverse_complement(self, seq: str) -> str:
        """Return reverse complement of a DNA sequence."""
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return "".join(complement.get(base, base) for base in reversed(seq))

    def _validate_fastq_record(self, header: str, sequence: str, plus: str, quality: str, line_num: int, filepath: str) -> bool:
        """Validate a single FASTQ record."""
        if not header.startswith('@'):
            logging.error(f"Invalid FASTQ header at line {line_num} in {filepath}: {header[:50]}...")
            return False
        
        if not plus.startswith('+'):
            logging.error(f"Invalid FASTQ plus line at line {line_num + 2} in {filepath}: {plus[:50]}...")
            return False
        
        if len(sequence) != len(quality):
            logging.error(f"Sequence and quality length mismatch at line {line_num} in {filepath}")
            return False
        
        # Check for valid DNA bases
        valid_bases = set('ATCGN')
        if not all(base in valid_bases for base in sequence.upper()):
            logging.warning(f"Non-standard bases found in sequence at line {line_num} in {filepath}")
        
        return True

    def _stream_fastq_records(self, filepath: Union[str, List[str]]):
        """
        Generator that yields FASTQ sequences one at a time for memory efficiency.
        
        Args:
            filepath: Single file path (str) or list of paired-end file paths
            
        Yields:
            DNA sequences as strings
        """
        read_count = 0
        file_paths = [filepath] if isinstance(filepath, str) else filepath
        
        for file_path in file_paths:
            line_num = 0
            try:
                with self._open_file(file_path) as f:
                    while True:
                        if self.max_reads and read_count >= self.max_reads:
                            break

                        line_num += 1
                        header = f.readline().strip()
                        if not header:
                            break

                        line_num += 1
                        sequence = f.readline().strip()
                        line_num += 1
                        plus = f.readline().strip()
                        line_num += 1
                        quality = f.readline().strip()

                        if header and sequence and plus and quality:
                            if self._validate_fastq_record(header, sequence, plus, quality, line_num - 3, file_path):
                                yield sequence.upper()
                                read_count += 1
                        elif header:  # Incomplete record at end of file
                            logging.warning(f"Incomplete FASTQ record at end of {file_path}")
                            break
                            
            except Exception as e:
                raise ValueError(f"Error parsing FASTQ file {file_path} at line {line_num}: {str(e)}")

    def _parse_fastq(self, filepath: Union[str, List[str]]) -> List[str]:
        """
        Parse FASTQ file(s) and return sequences with validation.
        
        Args:
            filepath: Single file path (str) or list of paired-end file paths
            
        Returns:
            List of sequences from all input files
        """
        sequences = []
        read_count = 0
        
        # Handle both single files and paired-end file lists
        file_paths = [filepath] if isinstance(filepath, str) else filepath
        
        for file_path in file_paths:
            line_num = 0
            try:
                with self._open_file(file_path) as f:
                    while True:
                        if self.max_reads and read_count >= self.max_reads:
                            break

                        line_num += 1
                        header = f.readline().strip()
                        if not header:
                            break

                        line_num += 1
                        sequence = f.readline().strip()
                        line_num += 1
                        plus = f.readline().strip()
                        line_num += 1
                        quality = f.readline().strip()

                        if header and sequence and plus and quality:
                            if self._validate_fastq_record(header, sequence, plus, quality, line_num - 3, file_path):
                                sequences.append(sequence.upper())
                                read_count += 1
                        elif header:  # Incomplete record at end of file
                            logging.warning(f"Incomplete FASTQ record at end of {file_path}")
                            break
                            
            except Exception as e:
                raise ValueError(f"Error parsing FASTQ file {file_path} at line {line_num}: {str(e)}")

        return sequences

    def profile_sample(self, filepath: Union[str, List[str]], sample_name: str, 
                      memory_efficient: bool = True) -> Dict[str, int]:
        """
        Generate k-mer profile for a single sample (single-end or paired-end).
        
        Args:
            filepath: File path(s) to process
            sample_name: Name of the sample
            memory_efficient: Use streaming parser to reduce memory usage
        """
        if isinstance(filepath, list):
            logging.info(f"Processing paired-end sample: {sample_name} ({len(filepath)} files)")
        else:
            logging.info(f"Processing single-end sample: {sample_name}")

        total_kmers = defaultdict(int)
        sequence_count = 0
        
        if memory_efficient:
            # Use streaming parser for memory efficiency
            for sequence in self._stream_fastq_records(filepath):
                seq_kmers = self._extract_kmers(sequence)
                for kmer, count in seq_kmers.items():
                    total_kmers[kmer] += count
                sequence_count += 1
        else:
            # Load all sequences into memory (legacy behavior)
            sequences = self._parse_fastq(filepath)
            sequence_count = len(sequences)
            for seq in sequences:
                seq_kmers = self._extract_kmers(seq)
                for kmer, count in seq_kmers.items():
                    total_kmers[kmer] += count

        logging.info(f"Processed {sequence_count} sequences from {sample_name}")

        # Filter low-frequency k-mers to reduce memory usage
        if self.min_kmer_freq > 1:
            filtered_kmers = {kmer: count for kmer, count in total_kmers.items() 
                            if count >= self.min_kmer_freq}
            logging.info(f"Filtered {len(total_kmers) - len(filtered_kmers)} low-frequency k-mers")
            total_kmers = filtered_kmers

        # Normalize by total k-mer count
        total_count = sum(total_kmers.values())
        if total_count == 0:
            logging.warning(f"No k-mers found for sample {sample_name}")
            normalized_profile = {}
        else:
            normalized_profile = {
                kmer: count / total_count for kmer, count in total_kmers.items()
            }

        self.profiles[sample_name] = normalized_profile
        if sample_name not in self.sample_names:
            self.sample_names.append(sample_name)

        logging.info(f"Generated profile with {len(normalized_profile)} unique k-mers")
        return normalized_profile

    def process_samples_parallel(self, fastq_files: List[Tuple[Union[str, List[str]], str]], 
                                n_processes: Optional[int] = None, 
                                show_progress: bool = True,
                                memory_efficient: bool = True) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
        """
        Process multiple samples in parallel using multiprocessing.
        
        Args:
            fastq_files: List of (filepath, sample_name) tuples
            n_processes: Number of processes to use (default: CPU count)
            show_progress: Whether to show progress updates
            
        Returns:
            Tuple containing (profiles_dict, failed_samples_list)
        """
        if n_processes is None:
            n_processes = min(mp.cpu_count(), len(fastq_files))
        
        logging.info(f"Processing {len(fastq_files)} samples using {n_processes} processes")
        
        # Create shared progress tracking
        manager = Manager()
        progress_dict = manager.dict({
            'completed': 0,
            'failed': 0,
            'current_sample': '',
            'errors': manager.dict()
        })
        lock = manager.Lock()
        
        # Prepare arguments for worker processes
        worker_args = [
            (filepath, sample_name, self.k, self.max_reads, self.min_kmer_freq, memory_efficient, self.use_fast_extraction, progress_dict, lock)
            for filepath, sample_name in fastq_files
        ]
        
        start_time = time.time()
        failed_samples = []
        
        try:
            with Pool(processes=n_processes) as pool:
                if show_progress:
                    # Start progress monitoring
                    import threading
                    progress_thread = threading.Thread(
                        target=self._monitor_progress,
                        args=(progress_dict, len(fastq_files), start_time)
                    )
                    progress_thread.daemon = True
                    progress_thread.start()
                
                # Process samples in parallel
                results = pool.map(process_sample_worker, worker_args)
                
                # Collect results
                for sample_name, profile, error in results:
                    if error is None:
                        self.profiles[sample_name] = profile
                        if sample_name not in self.sample_names:
                            self.sample_names.append(sample_name)
                    else:
                        failed_samples.append(sample_name)
                        logging.error(f"Failed to process {sample_name}: {error}")
        
        except KeyboardInterrupt:
            logging.warning("Processing interrupted by user")
            raise
        except Exception as e:
            logging.error(f"Error in parallel processing: {e}")
            raise
        
        elapsed_time = time.time() - start_time
        success_count = len(fastq_files) - len(failed_samples)
        
        logging.info(f"Parallel processing completed in {elapsed_time:.1f}s")
        logging.info(f"Successfully processed: {success_count}/{len(fastq_files)} samples")
        
        if failed_samples:
            logging.warning(f"Failed samples: {', '.join(failed_samples)}")
        
        return self.profiles, failed_samples

    def _monitor_progress(self, progress_dict, total_samples: int, start_time: float):
        """Monitor and display progress updates."""
        while progress_dict['completed'] < total_samples:
            completed = progress_dict['completed']
            failed = progress_dict['failed']
            current = progress_dict.get('current_sample', '')
            
            if completed > 0:
                elapsed = time.time() - start_time
                rate = completed / elapsed
                eta = (total_samples - completed) / rate if rate > 0 else 0
                
                progress_msg = f"Progress: {completed}/{total_samples} ({completed/total_samples*100:.1f}%)"
                if failed > 0:
                    progress_msg += f" - {failed} failed"
                if current:
                    progress_msg += f" - Current: {current}"
                progress_msg += f" - ETA: {eta:.1f}s"
                
                logging.info(progress_msg)
            
            time.sleep(2)  # Update every 2 seconds