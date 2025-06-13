#!/usr/bin/env python3
"""
Compression-based similarity profiler for MetaGrouper.

This module implements Normalized Compression Distance (NCD) and related
compression-based similarity measures for metagenomic sequence analysis.

Based on:
- Li et al. (2004) - The similarity metric
- Chen et al. (2000) - An information-based sequence distance
- Rudi Cilibrasi & Paul Vitanyi (2005) - Clustering by compression
"""

import gzip
import lzma
import zlib
import bz2
import zstandard as zstd
import logging
import tempfile
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
from collections import defaultdict
import numpy as np
import multiprocessing as mp
from multiprocessing import Pool, Manager
import time


class CompressionProfiler:
    """
    Compression-based similarity profiler using Normalized Compression Distance (NCD).
    
    This class replaces k-mer based profiling with compression-based similarity
    measures that are faster and parameter-free.
    """
    
    def __init__(self, 
                 compressor: str = 'zstd',
                 max_reads: Optional[int] = None,
                 sample_sequences: bool = True,
                 compression_level: int = 3):
        """
        Initialize compression profiler.
        
        Args:
            compressor: Compression algorithm ('gzip', 'zstd', 'lzma', 'bz2')
            max_reads: Maximum number of reads to process per sample
            sample_sequences: Whether to sample sequences for speed
            compression_level: Compression level (higher = better compression, slower)
        """
        self.compressor = compressor
        self.max_reads = max_reads
        self.sample_sequences = sample_sequences
        self.compression_level = compression_level
        
        # Storage for sample data and results
        self.sample_data = {}  # sample_name -> compressed sequence data
        self.sample_names = []
        self.compression_sizes = {}  # sample_name -> compressed size
        
        # Set up compressor function
        self._setup_compressor()
        
        logging.info(f"Initialized CompressionProfiler with {compressor} compressor")
    
    def _setup_compressor(self):
        """Set up the compression function based on selected algorithm."""
        if self.compressor == 'gzip':
            self.compress_func = lambda data: gzip.compress(data, compresslevel=self.compression_level)
        elif self.compressor == 'zstd':
            cctx = zstd.ZstdCompressor(level=self.compression_level)
            self.compress_func = lambda data: cctx.compress(data)
        elif self.compressor == 'lzma':
            self.compress_func = lambda data: lzma.compress(data, preset=self.compression_level)
        elif self.compressor == 'bz2':
            self.compress_func = lambda data: bz2.compress(data, compresslevel=self.compression_level)
        else:
            raise ValueError(f"Unsupported compressor: {self.compressor}")
    
    def _stream_fastq_records(self, filepath: Union[str, List[str]]) -> str:
        """
        Stream FASTQ records and extract sequences only.
        
        Args:
            filepath: Path to FASTQ file or list of paths for paired-end
            
        Yields:
            sequence: DNA sequence string
        """
        files_to_process = [filepath] if isinstance(filepath, str) else filepath
        
        sequence_count = 0
        
        for file_path in files_to_process:
            try:
                if file_path.endswith('.gz'):
                    file_handle = gzip.open(file_path, 'rt')
                else:
                    file_handle = open(file_path, 'r')
                
                with file_handle as f:
                    lines = []
                    for line in f:
                        lines.append(line.strip())
                        
                        # Process complete FASTQ record (4 lines)
                        if len(lines) == 4:
                            if lines[0].startswith('@') and lines[2].startswith('+'):
                                sequence = lines[1]  # Extract just the sequence
                                yield sequence
                                sequence_count += 1
                                
                                # Respect max_reads limit
                                if self.max_reads and sequence_count >= self.max_reads:
                                    return
                            
                            lines = []
                            
            except Exception as e:
                logging.warning(f"Error reading {file_path}: {e}")
                continue
    
    def _concatenate_sequences(self, filepath: Union[str, List[str]], sample_name: str) -> bytes:
        """
        Concatenate sequences from FASTQ file(s) into a single string for compression.
        
        Args:
            filepath: Path to FASTQ file or list of paths for paired-end
            sample_name: Name of the sample
            
        Returns:
            bytes: Concatenated sequence data as bytes
        """
        sequences = []
        sequence_count = 0
        
        for sequence in self._stream_fastq_records(filepath):
            sequences.append(sequence)
            sequence_count += 1
        
        # Join all sequences with a separator
        concatenated = '|'.join(sequences)
        
        logging.info(f"Concatenated {sequence_count} sequences from {sample_name}")
        
        return concatenated.encode('utf-8')
    
    def profile_sample(self, filepath: Union[str, List[str]], sample_name: str) -> int:
        """
        Profile a single sample by compressing its sequence data.
        
        Args:
            filepath: Path to FASTQ file or list of paths for paired-end
            sample_name: Name of the sample
            
        Returns:
            int: Size of compressed data
        """
        if isinstance(filepath, list):
            logging.info(f"Processing paired-end sample: {sample_name} ({len(filepath)} files)")
        else:
            logging.info(f"Processing single-end sample: {sample_name}")
        
        # Concatenate sequences
        sequence_data = self._concatenate_sequences(filepath, sample_name)
        
        # Compress the data
        compressed_data = self.compress_func(sequence_data)
        compressed_size = len(compressed_data)
        
        # Store results
        self.sample_data[sample_name] = compressed_data
        self.compression_sizes[sample_name] = compressed_size
        
        if sample_name not in self.sample_names:
            self.sample_names.append(sample_name)
        
        logging.info(f"Compressed {sample_name}: {len(sequence_data)} -> {compressed_size} bytes "
                    f"(ratio: {compressed_size/len(sequence_data):.3f})")
        
        return compressed_size
    
    def compute_ncd(self, sample1: str, sample2: str) -> float:
        """
        Compute Normalized Compression Distance between two samples.
        
        NCD(x,y) = (C(xy) - min(C(x),C(y))) / max(C(x),C(y))
        
        Args:
            sample1, sample2: Names of samples to compare
            
        Returns:
            float: NCD value between 0 and 1
        """
        if sample1 not in self.sample_data or sample2 not in self.sample_data:
            raise ValueError(f"Samples {sample1} or {sample2} not found in profile data")
        
        # Get individual compressed sizes
        c_x = self.compression_sizes[sample1]
        c_y = self.compression_sizes[sample2]
        
        # Concatenate and compress together
        data_xy = self.sample_data[sample1] + b'|' + self.sample_data[sample2]
        c_xy = len(self.compress_func(data_xy))
        
        # Compute NCD
        ncd = (c_xy - min(c_x, c_y)) / max(c_x, c_y)
        
        return max(0.0, min(1.0, ncd))  # Clamp to [0, 1]
    
    def compute_distance_matrix(self) -> np.ndarray:
        """
        Compute pairwise NCD distance matrix for all samples.
        
        Returns:
            np.ndarray: Distance matrix of shape (n_samples, n_samples)
        """
        n_samples = len(self.sample_names)
        distance_matrix = np.zeros((n_samples, n_samples))
        
        logging.info(f"Computing {self.compressor}-based NCD matrix for {n_samples} samples")
        
        for i, sample1 in enumerate(self.sample_names):
            for j, sample2 in enumerate(self.sample_names):
                if i == j:
                    distance_matrix[i, j] = 0.0
                elif i < j:  # Only compute upper triangle
                    ncd_value = self.compute_ncd(sample1, sample2)
                    distance_matrix[i, j] = ncd_value
                    distance_matrix[j, i] = ncd_value  # Symmetric
        
        return distance_matrix
    
    def process_samples_parallel(self, 
                                fastq_files: List[Tuple[Union[str, List[str]], str]], 
                                n_processes: Optional[int] = None, 
                                show_progress: bool = True) -> Tuple[Dict[str, int], List[str]]:
        """
        Process multiple samples in parallel using compression-based profiling.
        
        Args:
            fastq_files: List of (filepath, sample_name) tuples
            n_processes: Number of processes to use (default: CPU count)
            show_progress: Whether to show progress updates
            
        Returns:
            Tuple containing (compression_sizes_dict, failed_samples_list)
        """
        if n_processes is None:
            n_processes = min(mp.cpu_count(), len(fastq_files))
        
        logging.info(f"Processing {len(fastq_files)} samples using {n_processes} processes "
                    f"with {self.compressor} compression")
        
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
            (filepath, sample_name, self.compressor, self.max_reads, 
             self.compression_level, progress_dict, lock)
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
                results = pool.map(compression_worker, worker_args)
                
                # Collect results
                for sample_name, compressed_size, compressed_data, error in results:
                    if error is None:
                        self.compression_sizes[sample_name] = compressed_size
                        self.sample_data[sample_name] = compressed_data
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
        
        logging.info(f"Compression-based processing completed in {elapsed_time:.1f}s")
        logging.info(f"Successfully processed: {success_count}/{len(fastq_files)} samples")
        
        if failed_samples:
            logging.warning(f"Failed samples: {', '.join(failed_samples)}")
        
        return self.compression_sizes, failed_samples
    
    def _monitor_progress(self, progress_dict, total_samples: int, start_time: float):
        """Monitor and log progress of parallel processing."""
        self.last_progress_time = 0
        while progress_dict['completed'] < total_samples:
            completed = progress_dict['completed']
            failed = progress_dict['failed']
            current = progress_dict.get('current_sample', '')
            
            elapsed = time.time() - start_time
            if completed > 0:
                avg_time = elapsed / completed
                eta = avg_time * (total_samples - completed)
                progress_msg = (f"Progress: {completed}/{total_samples} "
                              f"({100*completed/total_samples:.1f}%) - "
                              f"Current: {current} - ETA: {eta:.1f}s")
            else:
                progress_msg = f"Starting compression-based processing of {total_samples} samples..."
            
            # Only log progress every 30 seconds or 10% completion
            if (elapsed > self.last_progress_time + 30 or 
                completed > 0 and completed % max(1, total_samples // 10) == 0):
                logging.info(progress_msg)
                self.last_progress_time = elapsed
            
            time.sleep(5)  # Check every 5 seconds but log less frequently


def compression_worker(args):
    """
    Worker function for multiprocessing compression-based sample processing.
    
    Args:
        args: Tuple containing (filepath, sample_name, compressor, max_reads, 
              compression_level, progress_dict, lock)
    
    Returns:
        Tuple containing (sample_name, compressed_size, compressed_data, error_message)
    """
    filepath, sample_name, compressor, max_reads, compression_level, progress_dict, lock = args
    
    try:
        # Create a temporary profiler for this worker
        profiler = CompressionProfiler(
            compressor=compressor, 
            max_reads=max_reads, 
            compression_level=compression_level
        )
        
        # Get sequence data
        sequence_data = profiler._concatenate_sequences(filepath, sample_name)
        
        # Compress the data
        compressed_data = profiler.compress_func(sequence_data)
        compressed_size = len(compressed_data)
        
        # Update progress safely
        with lock:
            progress_dict['completed'] += 1
            progress_dict['current_sample'] = sample_name
        
        return (sample_name, compressed_size, compressed_data, None)
    
    except Exception as e:
        # Update progress and record error
        with lock:
            progress_dict['completed'] += 1
            progress_dict['failed'] += 1
            progress_dict['errors'][sample_name] = str(e)
        
        return (sample_name, None, None, str(e))