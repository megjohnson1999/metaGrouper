#!/usr/bin/env python3
"""
Sourmash-based k-mer profiler for MetaGrouper.

This module provides a fast alternative to the built-in k-mer profiling
using sourmash's MinHash sketches for improved performance on large datasets.
"""

import csv
import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import sourmash
from sourmash import SourmashSignature, MinHash
from sourmash.index import LinearIndex
from sourmash.search import SearchResult
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import json


class SourmashProfiler:
    """K-mer profiler using sourmash MinHash sketches."""
    
    def __init__(self, 
                 k: int = 21,
                 scaled: int = 100,
                 num_hashes: int = 0,
                 seed: int = 42,
                 processes: int = 1,
                 track_abundance: bool = False,
                 dna: bool = True,
                 dayhoff: bool = False,
                 hp: bool = False,
                 protein: bool = False,
                 additional_k_sizes: Optional[List[int]] = None):
        """
        Initialize sourmash profiler.
        
        Args:
            k: K-mer size (default: 21)
            scaled: Scaled sketch (1 in N hashes kept, default: 100 for high sensitivity)
            num_hashes: Number of hashes to keep (0 for scaled)
            seed: Random seed for MinHash
            processes: Number of parallel processes
            track_abundance: Track k-mer abundances (default: False, more robust to PCR bias)
            dna: DNA alphabet (default)
            dayhoff: Dayhoff alphabet
            hp: Hydrophobic-polar alphabet
            protein: Protein alphabet
            additional_k_sizes: Additional k-mer sizes to compute (e.g., [31, 51])
        """
        self.k = k
        self.scaled = scaled if num_hashes == 0 else 0
        self.num_hashes = num_hashes if num_hashes > 0 else 0
        self.seed = seed
        self.processes = processes
        self.track_abundance = track_abundance
        self.additional_k_sizes = additional_k_sizes or []
        
        # Set default additional k-sizes for multi-scale analysis
        if not self.additional_k_sizes and scaled <= 100:
            self.additional_k_sizes = [31, 51]  # Multi-scale like your previous analysis
        
        # Set molecule type
        if dna:
            self.moltype = 'DNA'
        elif dayhoff:
            self.moltype = 'dayhoff'
        elif hp:
            self.moltype = 'hp'
        elif protein:
            self.moltype = 'protein'
        else:
            self.moltype = 'DNA'
            
        logging.info(f"Initialized SourmashProfiler: k={k}, scaled={scaled}, "
                    f"num_hashes={num_hashes}, moltype={self.moltype}, "
                    f"track_abundance={track_abundance}, additional_k_sizes={self.additional_k_sizes}")
    
    def create_minhash(self, ksize: Optional[int] = None) -> MinHash:
        """Create a new MinHash object with current parameters."""
        k = ksize or self.k
        if self.scaled > 0:
            return MinHash(n=0, ksize=k, scaled=self.scaled, 
                          seed=self.seed, track_abundance=self.track_abundance,
                          is_protein=(self.moltype == 'protein'),
                          dayhoff=(self.moltype == 'dayhoff'),
                          hp=(self.moltype == 'hp'))
        else:
            return MinHash(n=self.num_hashes, ksize=k, 
                          seed=self.seed, track_abundance=self.track_abundance,
                          is_protein=(self.moltype == 'protein'),
                          dayhoff=(self.moltype == 'dayhoff'),
                          hp=(self.moltype == 'hp'))
    
    def sketch_sample(self, filepath: Union[str, List[str]], 
                     sample_name: Optional[str] = None) -> List[SourmashSignature]:
        """
        Create sourmash signatures for a sample (multi-scale analysis).
        
        Args:
            filepath: Path to FASTQ file(s)
            sample_name: Name for the signature
            
        Returns:
            List of SourmashSignature objects (one per k-mer size)
        """
        # Create multiple MinHash objects for different k-mer sizes
        k_sizes = [self.k] + self.additional_k_sizes
        minhashes = {k: self.create_minhash(k) for k in k_sizes}
        
        # Handle both single files and paired-end file lists
        file_paths = [filepath] if isinstance(filepath, str) else filepath
        
        # For paired-end files, concatenate them like in manual script
        if isinstance(filepath, list) and len(filepath) == 2:
            # This is paired-end data - concatenate R1 and R2 like manual script
            logging.debug(f"Processing paired-end sample: {filepath[0]} + {filepath[1]}")
            
            # Create temporary combined file
            import tempfile
            import gzip
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as temp_combined:
                temp_combined_path = temp_combined.name
            
            try:
                # Concatenate R1 and R2 files like in manual script: cat ${r1_file} ${r2_file} > combined.fastq
                with open(temp_combined_path, 'w') as outfile:
                    for file_path in file_paths:
                        logging.debug(f"Adding {file_path} to combined file")
                        
                        # Handle both gzipped and plain text files
                        if file_path.endswith('.gz'):
                            with gzip.open(file_path, 'rt', encoding='utf-8') as infile:
                                outfile.write(infile.read())
                        else:
                            with open(file_path, 'r', encoding='utf-8') as infile:
                                outfile.write(infile.read())
                
                # Process the combined file
                import screed
                for record in screed.open(temp_combined_path):
                    # Add sequence to all k-mer sizes
                    for mh in minhashes.values():
                        mh.add_sequence(record.sequence, force=True)
                
            finally:
                # Clean up temporary file
                try:
                    os.unlink(temp_combined_path)
                except:
                    pass
        else:
            # Single-end or single file processing
            for file_path in file_paths:
                logging.debug(f"Processing {file_path}")
                
                # Use screed to parse FASTQ files (screed handles gzipped files automatically)
                import screed
                for record in screed.open(file_path):
                    # Add sequence to all k-mer sizes
                    for mh in minhashes.values():
                        mh.add_sequence(record.sequence, force=True)
        
        # Create signatures for all k-mer sizes
        if sample_name is None:
            sample_name = Path(file_paths[0]).stem
            
        signatures = []
        for k, mh in minhashes.items():
            sig_name = f"{sample_name}_k{k}" if len(minhashes) > 1 else sample_name
            sig = SourmashSignature(mh, name=sig_name)
            signatures.append(sig)
        
        return signatures
    
    def process_samples_parallel(self, samples: Union[Dict[str, Union[str, List[str]]], List[Tuple[Union[str, List[str]], str]]]) -> Dict[str, List[SourmashSignature]]:
        """
        Process multiple samples in parallel.
        
        Args:
            samples: Dictionary mapping sample names to file paths, or MetaGrouper format list
                    MetaGrouper format: List[Tuple[Union[str, List[str]], str]]
                    Where each tuple is (filepath_or_list, sample_name)
            
        Returns:
            Dictionary mapping sample names to lists of signatures (one per k-mer size)
        """
        signatures = {}
        
        # Handle MetaGrouper format vs dict format
        if isinstance(samples, dict):
            sample_dict = samples
        else:
            # Convert MetaGrouper list format to dict
            sample_dict = {}
            for item in samples:
                # MetaGrouper format: each item is (filepath_or_list, sample_name)
                if isinstance(item, (tuple, list)) and len(item) == 2:
                    filepath, sample_name = item
                    sample_dict[sample_name] = filepath
                else:
                    # Fallback for unexpected format
                    logging.warning(f"Unexpected sample format: {item}")
                    if isinstance(item, str):
                        # Single file path without sample name
                        filepath = item
                        sample_name = Path(item).stem
                        sample_dict[sample_name] = filepath
                    else:
                        logging.error(f"Cannot process sample: {item}")
                        continue
        
        if self.processes == 1:
            # Single process
            for sample_name, filepath in sample_dict.items():
                logging.info(f"Processing {sample_name}")
                try:
                    sigs = self.sketch_sample(filepath, sample_name)
                    signatures[sample_name] = sigs
                except Exception as e:
                    logging.error(f"Error processing {sample_name}: {e}")
        else:
            # Parallel processing
            with ProcessPoolExecutor(max_workers=self.processes) as executor:
                future_to_sample = {
                    executor.submit(self.sketch_sample, filepath, sample_name): sample_name
                    for sample_name, filepath in sample_dict.items()
                }
                
                for future in as_completed(future_to_sample):
                    sample_name = future_to_sample[future]
                    try:
                        sigs = future.result()
                        signatures[sample_name] = sigs
                        logging.info(f"Completed {sample_name}")
                    except Exception as e:
                        logging.error(f"Error processing {sample_name}: {e}")
        
        return signatures
    
    def compute_similarity_matrix(self, signatures: Dict[str, List[SourmashSignature]], 
                                 use_k_size: Optional[int] = None) -> np.ndarray:
        """
        Compute pairwise Jaccard similarity matrix using sourmash compare.
        
        Args:
            signatures: Dictionary of sourmash signature lists
            use_k_size: Specific k-mer size to use (default: primary k-mer size)
            
        Returns:
            Similarity matrix as numpy array
        """
        sample_names = list(signatures.keys())
        n_samples = len(sample_names)
        
        # Select signatures for the specified k-mer size
        if use_k_size is None:
            use_k_size = self.k
            
        selected_signatures = []
        for sample_name in sample_names:
            sig_list = signatures[sample_name]
            # Find signature with matching k-mer size
            matching_sig = None
            for sig in sig_list:
                if sig.minhash.ksize == use_k_size:
                    matching_sig = sig
                    break
            if matching_sig is None:
                # Fall back to first signature if no exact match
                matching_sig = sig_list[0]
            selected_signatures.append(matching_sig)
        
        # Create temporary file for signatures
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sig', delete=False) as temp_sig_file:
            temp_sig_path = temp_sig_file.name
            sourmash.save_signatures(selected_signatures, temp_sig_file)
        
        # Create temporary file for output matrix
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as temp_matrix_file:
            temp_matrix_path = temp_matrix_file.name
        
        try:
            # Run sourmash compare with CSV output
            cmd = ['sourmash', 'compare', temp_sig_path, '--csv', temp_matrix_path]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logging.debug(f"sourmash compare completed: {result.stderr}")
            
            # Read the similarity matrix from CSV
            similarity_matrix = np.zeros((n_samples, n_samples))
            with open(temp_matrix_path, 'r') as f:
                reader = csv.reader(f)
                # Skip header row (sample names)
                next(reader, None)
                for i, row in enumerate(reader):
                    # Skip the first column (sample name) and read the similarity values
                    for j, value in enumerate(row[1:]):
                        similarity_matrix[i, j] = float(value)
            
            return similarity_matrix
            
        except subprocess.CalledProcessError as e:
            logging.error(f"sourmash compare failed: {e.stderr}")
            # Fall back to manual computation
            logging.warning("Falling back to manual similarity computation")
            return self._compute_similarity_matrix_manual(dict(zip(sample_names, selected_signatures)))
        
        finally:
            # Clean up temporary files
            try:
                os.unlink(temp_sig_path)
                os.unlink(temp_matrix_path)
            except:
                pass
    
    def _compute_similarity_matrix_manual(self, signatures: Dict[str, SourmashSignature]) -> np.ndarray:
        """
        Manual computation of similarity matrix as fallback.
        
        Args:
            signatures: Dictionary of sourmash signatures
            
        Returns:
            Similarity matrix as numpy array
        """
        sample_names = list(signatures.keys())
        n_samples = len(sample_names)
        similarity_matrix = np.zeros((n_samples, n_samples))
        
        # Fill diagonal with 1s
        np.fill_diagonal(similarity_matrix, 1.0)
        
        # Compute pairwise similarities
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                sig1 = signatures[sample_names[i]]
                sig2 = signatures[sample_names[j]]
                
                # Compute Jaccard similarity
                similarity = sig1.jaccard(sig2)
                similarity_matrix[i, j] = similarity
                similarity_matrix[j, i] = similarity
                
        return similarity_matrix
    
    def save_signatures(self, signatures: Dict[str, List[SourmashSignature]], 
                       output_path: str) -> None:
        """
        Save signatures to a file.
        
        Args:
            signatures: Dictionary of signature lists
            output_path: Output file path (.sig or .sig.gz)
        """
        # Flatten all signatures into one list
        all_signatures = []
        for sig_list in signatures.values():
            all_signatures.extend(sig_list)
        
        with open(output_path, 'w') as f:
            sourmash.save_signatures(all_signatures, f)
        logging.info(f"Saved {len(all_signatures)} signatures ({len(signatures)} samples) to {output_path}")
    
    def load_signatures(self, signature_path: str) -> Dict[str, SourmashSignature]:
        """
        Load signatures from a file.
        
        Args:
            signature_path: Path to signature file
            
        Returns:
            Dictionary mapping sample names to signatures
        """
        signatures = {}
        
        for sig in sourmash.load_file_as_signatures(signature_path):
            signatures[sig.name] = sig
                
        logging.info(f"Loaded {len(signatures)} signatures from {signature_path}")
        return signatures
    
    def export_to_metagrouper_format(self, signatures: Dict[str, List[SourmashSignature]],
                                    similarity_matrix: np.ndarray, 
                                    use_k_size: Optional[int] = None,
                                    prevalence_threshold: float = 0.1) -> Tuple[Dict, List[str]]:
        """
        Convert sourmash results to MetaGrouper's expected format.
        
        Args:
            signatures: Dictionary of signature lists
            similarity_matrix: Similarity matrix
            use_k_size: Specific k-mer size to use (default: primary k-mer size)
            
        Returns:
            Tuple of (profiles dict, sample names list)
        """
        sample_names = list(signatures.keys())
        
        # Select signatures for the specified k-mer size
        if use_k_size is None:
            use_k_size = self.k
        
        # Convert signatures to a format compatible with MetaGrouper
        # We'll use the hash values as "k-mers" for compatibility
        profiles = {}
        
        for sample_name, sig_list in signatures.items():
            # Find signature with matching k-mer size
            matching_sig = None
            for sig in sig_list:
                if sig.minhash.ksize == use_k_size:
                    matching_sig = sig
                    break
            if matching_sig is None:
                # Fall back to first signature if no exact match
                matching_sig = sig_list[0]
            
            # Get the MinHash object
            mh = matching_sig.minhash
            
            # Get hashes as a proxy for k-mers
            if self.track_abundance:
                # Use abundance information if available
                hashes = mh.hashes
                profiles[sample_name] = {str(h): count for h, count in hashes.items()}
            else:
                # Just use presence/absence
                hashes = mh.hashes
                profiles[sample_name] = {str(h): 1 for h in hashes}
        
        # Apply prevalence filtering to reduce memory usage
        if prevalence_threshold > 0 and len(sample_names) > 10:
            from collections import defaultdict
            
            # Count how many samples each hash appears in
            hash_counts = defaultdict(int)
            for profile in profiles.values():
                for hash_val in profile.keys():
                    hash_counts[hash_val] += 1
            
            # Early warning for very large feature spaces
            total_unique_hashes = len(hash_counts)
            if total_unique_hashes > 50000000:  # 50M threshold
                logging.warning(f"Very large feature space detected: {total_unique_hashes:,} unique hashes")
                logging.warning("Consider using --aggressive-mode or higher --scaled value for better memory efficiency")
            
            # Filter hashes by prevalence
            n_samples = len(sample_names)
            
            # Auto-adjust prevalence threshold for very large datasets
            adjusted_threshold = prevalence_threshold
            if n_samples > 500 and prevalence_threshold < 0.2:
                adjusted_threshold = 0.2  # 20% for large datasets
                logging.info(f"Auto-adjusted prevalence threshold to {adjusted_threshold:.1%} for large dataset ({n_samples} samples)")
            elif n_samples > 1000 and prevalence_threshold < 0.3:
                adjusted_threshold = 0.3  # 30% for very large datasets
                logging.info(f"Auto-adjusted prevalence threshold to {adjusted_threshold:.1%} for very large dataset ({n_samples} samples)")
            
            min_samples = max(2, int(adjusted_threshold * n_samples))
            common_hashes = {h for h, count in hash_counts.items() if count >= min_samples}
            
            # Apply filtering to all profiles
            filtered_profiles = {}
            for sample_name, profile in profiles.items():
                filtered_profiles[sample_name] = {h: count for h, count in profile.items() if h in common_hashes}
            
            # Log filtering results
            original_count = len(hash_counts)
            filtered_count = len(common_hashes)
            reduction_percent = 100 * (1 - filtered_count / original_count) if original_count > 0 else 0
            
            logging.info(f"Sourmash hash filtering: {original_count:,} total → {filtered_count:,} retained "
                        f"({reduction_percent:.1f}% reduction, min_samples={min_samples})")
            
            profiles = filtered_profiles
        
        return profiles, sample_names
    
    def create_analysis_summary(self, signatures: Dict[str, List[SourmashSignature]]) -> Dict:
        """
        Create analysis summary statistics.
        
        Args:
            signatures: Dictionary of signature lists
            
        Returns:
            Summary statistics dictionary
        """
        summary = {
            'num_samples': len(signatures),
            'k': self.k,
            'scaled': self.scaled,
            'num_hashes': self.num_hashes,
            'moltype': self.moltype,
            'track_abundance': self.track_abundance,
            'additional_k_sizes': self.additional_k_sizes,
            'samples': {}
        }
        
        for sample_name, sig_list in signatures.items():
            sample_info = {
                'num_signatures': len(sig_list),
                'k_sizes': [sig.minhash.ksize for sig in sig_list],
                'signatures': []
            }
            
            for sig in sig_list:
                mh = sig.minhash
                sample_info['signatures'].append({
                    'k_size': mh.ksize,
                    'num_hashes': len(mh),
                    'md5sum': sig.md5sum(),
                    'name': sig.name
                })
            
            summary['samples'][sample_name] = sample_info
            
        return summary