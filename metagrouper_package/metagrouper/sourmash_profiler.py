#!/usr/bin/env python3
"""
Sourmash-based k-mer profiler for MetaGrouper.

This module provides a fast alternative to the built-in k-mer profiling
using sourmash's MinHash sketches for improved performance on large datasets.
"""

import logging
import os
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
                 scaled: int = 1000,
                 num_hashes: int = 0,
                 seed: int = 42,
                 processes: int = 1,
                 track_abundance: bool = False,
                 dna: bool = True,
                 dayhoff: bool = False,
                 hp: bool = False,
                 protein: bool = False):
        """
        Initialize sourmash profiler.
        
        Args:
            k: K-mer size
            scaled: Scaled sketch (1 in N hashes kept)
            num_hashes: Number of hashes to keep (0 for scaled)
            seed: Random seed for MinHash
            processes: Number of parallel processes
            track_abundance: Track k-mer abundances
            dna: DNA alphabet (default)
            dayhoff: Dayhoff alphabet
            hp: Hydrophobic-polar alphabet
            protein: Protein alphabet
        """
        self.k = k
        self.scaled = scaled if num_hashes == 0 else 0
        self.num_hashes = num_hashes if num_hashes > 0 else 0
        self.seed = seed
        self.processes = processes
        self.track_abundance = track_abundance
        
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
                    f"num_hashes={num_hashes}, moltype={self.moltype}")
    
    def create_minhash(self) -> MinHash:
        """Create a new MinHash object with current parameters."""
        if self.scaled > 0:
            return MinHash(n=0, ksize=self.k, scaled=self.scaled, 
                          seed=self.seed, track_abundance=self.track_abundance,
                          is_protein=(self.moltype == 'protein'),
                          dayhoff=(self.moltype == 'dayhoff'),
                          hp=(self.moltype == 'hp'))
        else:
            return MinHash(n=self.num_hashes, ksize=self.k, 
                          seed=self.seed, track_abundance=self.track_abundance,
                          is_protein=(self.moltype == 'protein'),
                          dayhoff=(self.moltype == 'dayhoff'),
                          hp=(self.moltype == 'hp'))
    
    def sketch_sample(self, filepath: Union[str, List[str]], 
                     sample_name: Optional[str] = None) -> SourmashSignature:
        """
        Create sourmash signature for a sample.
        
        Args:
            filepath: Path to FASTQ file(s)
            sample_name: Name for the signature
            
        Returns:
            SourmashSignature object
        """
        mh = self.create_minhash()
        
        # Handle both single files and paired-end file lists
        file_paths = [filepath] if isinstance(filepath, str) else filepath
        
        for file_path in file_paths:
            logging.debug(f"Processing {file_path}")
            
            # Use screed to parse FASTQ files
            import screed
            for record in screed.open(file_path):
                mh.add_sequence(record.sequence, force=True)
        
        # Create signature
        if sample_name is None:
            sample_name = Path(file_paths[0]).stem
            
        sig = SourmashSignature(mh, name=sample_name)
        return sig
    
    def process_samples_parallel(self, samples: Union[Dict[str, Union[str, List[str]]], List[Union[str, Tuple[str, ...]]]]) -> Dict[str, SourmashSignature]:
        """
        Process multiple samples in parallel.
        
        Args:
            samples: Dictionary mapping sample names to file paths, or MetaGrouper format dict
            
        Returns:
            Dictionary mapping sample names to signatures
        """
        signatures = {}
        
        # Handle MetaGrouper format (dict with sample names as keys, file paths as values)
        # vs list format that needs to be converted to dict
        if isinstance(samples, dict):
            sample_dict = samples
        else:
            # Convert list to dict using basename as sample name
            sample_dict = {}
            for item in samples:
                # Handle both single files (str) and paired-end files (tuple/list)
                if isinstance(item, (tuple, list)):
                    # Paired-end files - use first file for sample name
                    filepath = item
                    sample_name = Path(item[0]).stem
                    # Remove common suffixes like _1, _R1, etc.
                    sample_name = sample_name.replace('_1', '').replace('_R1', '').replace('_r1', '')
                else:
                    # Single file
                    filepath = item
                    sample_name = Path(item).stem
                
                sample_dict[sample_name] = filepath
        
        if self.processes == 1:
            # Single process
            for sample_name, filepath in sample_dict.items():
                logging.info(f"Processing {sample_name}")
                try:
                    sig = self.sketch_sample(filepath, sample_name)
                    signatures[sample_name] = sig
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
                        sig = future.result()
                        signatures[sample_name] = sig
                        logging.info(f"Completed {sample_name}")
                    except Exception as e:
                        logging.error(f"Error processing {sample_name}: {e}")
        
        return signatures
    
    def compute_similarity_matrix(self, signatures: Dict[str, SourmashSignature]) -> np.ndarray:
        """
        Compute pairwise Jaccard similarity matrix.
        
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
    
    def save_signatures(self, signatures: Dict[str, SourmashSignature], 
                       output_path: str) -> None:
        """
        Save signatures to a file.
        
        Args:
            signatures: Dictionary of signatures
            output_path: Output file path (.sig or .sig.gz)
        """
        with sourmash.save_signatures(signatures.values(), output_path) as f:
            pass
        logging.info(f"Saved {len(signatures)} signatures to {output_path}")
    
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
    
    def export_to_metagrouper_format(self, signatures: Dict[str, SourmashSignature],
                                    similarity_matrix: np.ndarray) -> Tuple[Dict, List[str]]:
        """
        Convert sourmash results to MetaGrouper's expected format.
        
        Args:
            signatures: Dictionary of signatures
            similarity_matrix: Similarity matrix
            
        Returns:
            Tuple of (profiles dict, sample names list)
        """
        sample_names = list(signatures.keys())
        
        # Convert signatures to a format compatible with MetaGrouper
        # We'll use the hash values as "k-mers" for compatibility
        profiles = {}
        
        for sample_name, sig in signatures.items():
            # Get the MinHash object
            mh = sig.minhash
            
            # Get hashes as a proxy for k-mers
            if self.track_abundance:
                # Use abundance information if available
                hashes = mh.hashes
                profiles[sample_name] = {str(h): count for h, count in hashes.items()}
            else:
                # Just use presence/absence
                hashes = mh.hashes
                profiles[sample_name] = {str(h): 1 for h in hashes}
        
        return profiles, sample_names
    
    def create_analysis_summary(self, signatures: Dict[str, SourmashSignature]) -> Dict:
        """
        Create analysis summary statistics.
        
        Args:
            signatures: Dictionary of signatures
            
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
            'samples': {}
        }
        
        for sample_name, sig in signatures.items():
            mh = sig.minhash
            summary['samples'][sample_name] = {
                'num_hashes': len(mh),
                'md5sum': sig.md5sum(),
                'name': sig.name
            }
            
        return summary