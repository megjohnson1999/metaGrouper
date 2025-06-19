#!/usr/bin/env python3
"""
Updated test suite for the modular MetaGrouper architecture.

This test suite validates the functionality of the refactored modular components.
"""

import unittest
import tempfile
import shutil
import sys
from pathlib import Path

# Add the package to the path for testing
sys.path.insert(0, str(Path(__file__).parent))

from metagrouper import (
    KmerProfiler,
    SimilarityAnalyzer,
    Visualizer,
    find_fastq_files,
    setup_logging,
    save_results,
    MetaGrouperConfig
)
from metagrouper.utils import (
    validate_input_directory,
    create_output_directory,
    summarize_fastq_files,
    check_dependencies,
    get_system_info
)


class TestModularKmerProfiler(unittest.TestCase):
    """Test the modular KmerProfiler implementation."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.profiler = KmerProfiler(k=17, min_kmer_freq=1)
        
        # Create test FASTQ files
        self.create_test_fastq_files()
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def create_test_fastq_files(self):
        """Create test FASTQ files for testing."""
        # Single-end sample
        single_file = self.temp_dir / "sample_001.fastq"
        with open(single_file, "w") as f:
            f.write("@read1\nATCGATCGATCGATCGATCG\n+\n~~~~~~~~~~~~~~~~~~~~\n")
            f.write("@read2\nGCTAGCTAGCTAGCTAGCTA\n+\n~~~~~~~~~~~~~~~~~~~~\n")
        
        # Paired-end sample
        r1_file = self.temp_dir / "sample_002_R1.fastq"
        r2_file = self.temp_dir / "sample_002_R2.fastq"
        
        with open(r1_file, "w") as f:
            f.write("@read1\nATCGATCGATCGATCGATCG\n+\n~~~~~~~~~~~~~~~~~~~~\n")
            f.write("@read2\nTTCCGGAATTCCGGAATTCC\n+\n~~~~~~~~~~~~~~~~~~~~\n")
        
        with open(r2_file, "w") as f:
            f.write("@read1\nCGATCGATCGATCGATCGAT\n+\n~~~~~~~~~~~~~~~~~~~~\n")
            f.write("@read2\nGGAATTCCGGAATTCCGGAA\n+\n~~~~~~~~~~~~~~~~~~~~\n")
    
    def test_kmer_extraction(self):
        """Test k-mer extraction functionality."""
        sequence = "ATCGATCGATCGATCGATCG"
        kmers = self.profiler._extract_kmers(sequence)
        
        self.assertIsInstance(kmers, dict)
        self.assertTrue(len(kmers) > 0)
        
        # Check canonical k-mer representation
        for kmer in kmers.keys():
            self.assertEqual(len(kmer), 17)
            rev_comp = self.profiler._reverse_complement(kmer)
            self.assertEqual(kmer, min(kmer, rev_comp))
    
    def test_reverse_complement(self):
        """Test reverse complement functionality."""
        seq = "ATCG"
        rev_comp = self.profiler._reverse_complement(seq)
        self.assertEqual(rev_comp, "CGAT")
    
    def test_single_end_profiling(self):
        """Test single-end sample profiling."""
        single_file = str(self.temp_dir / "sample_001.fastq")
        profile = self.profiler.profile_sample(single_file, "sample_001")
        
        self.assertIsInstance(profile, dict)
        self.assertTrue(len(profile) > 0)
        self.assertAlmostEqual(sum(profile.values()), 1.0, places=5)
    
    def test_paired_end_profiling(self):
        """Test paired-end sample profiling."""
        paired_files = [
            str(self.temp_dir / "sample_002_R1.fastq"),
            str(self.temp_dir / "sample_002_R2.fastq")
        ]
        profile = self.profiler.profile_sample(paired_files, "sample_002")
        
        self.assertIsInstance(profile, dict)
        self.assertTrue(len(profile) > 0)
        self.assertAlmostEqual(sum(profile.values()), 1.0, places=5)
    
    def test_memory_efficient_processing(self):
        """Test memory-efficient processing mode."""
        single_file = str(self.temp_dir / "sample_001.fastq")
        
        # Test with memory_efficient=True
        profile_efficient = self.profiler.profile_sample(
            single_file, "sample_001", memory_efficient=True
        )
        
        # Test with memory_efficient=False
        profile_standard = self.profiler.profile_sample(
            single_file, "sample_001_std", memory_efficient=False
        )
        
        # Both should produce similar results
        self.assertEqual(len(profile_efficient), len(profile_standard))


class TestModularSimilarityAnalyzer(unittest.TestCase):
    """Test the modular SimilarityAnalyzer implementation."""
    
    def setUp(self):
        """Set up test environment."""
        # Create mock profiles
        self.profiles = {
            "sample1": {"ATCGATCGATCGATCGA": 0.5, "GCTAGCTAGCTAGCTAG": 0.5},
            "sample2": {"ATCGATCGATCGATCGA": 0.3, "TTCCGGAATTCCGGAAT": 0.7},
            "sample3": {"GCTAGCTAGCTAGCTAG": 0.4, "TTCCGGAATTCCGGAAT": 0.6}
        }
        self.analyzer = SimilarityAnalyzer(self.profiles)
    
    def test_distance_matrix_computation(self):
        """Test distance matrix computation."""
        distance_matrix = self.analyzer.compute_distance_matrix()
        
        self.assertEqual(distance_matrix.shape, (3, 3))
        
        # Diagonal should be 0 (distance from sample to itself)
        for i in range(3):
            self.assertAlmostEqual(distance_matrix[i, i], 0.0, places=5)
        
        # Matrix should be symmetric
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(distance_matrix[i, j], distance_matrix[j, i], places=5)
    
    def test_pca_analysis(self):
        """Test PCA analysis."""
        pca_result, pca = self.analyzer.perform_pca()
        
        self.assertEqual(pca_result.shape, (3, 2))
        self.assertTrue(len(pca.explained_variance_ratio_) == 2)
        self.assertTrue(sum(pca.explained_variance_ratio_) <= 1.0)
    
    def test_mds_analysis(self):
        """Test MDS analysis."""
        self.analyzer.compute_distance_matrix()
        mds_result = self.analyzer.perform_mds()
        
        self.assertEqual(mds_result.shape, (3, 2))


class TestModularVisualizer(unittest.TestCase):
    """Test the modular Visualizer implementation."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.sample_names = ["sample1", "sample2", "sample3"]
        self.visualizer = Visualizer(self.sample_names)
        
        # Create mock data
        import numpy as np
        self.distance_matrix = np.array([
            [0.0, 0.5, 0.8],
            [0.5, 0.0, 0.6],
            [0.8, 0.6, 0.0]
        ])
        self.pca_result = np.array([[1, 2], [3, 4], [5, 6]])
        self.mds_result = np.array([[1, 1], [2, 2], [3, 3]])
        
        # Create mock PCA object
        class MockPCA:
            explained_variance_ratio_ = [0.6, 0.3]
        self.pca = MockPCA()
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_distance_heatmap(self):
        """Test distance heatmap generation."""
        output_path = self.temp_dir / "test_heatmap.png"
        self.visualizer.plot_distance_heatmap(self.distance_matrix, str(output_path))
        self.assertTrue(output_path.exists())
    
    def test_pca_plot(self):
        """Test PCA plot generation."""
        output_path = self.temp_dir / "test_pca.png"
        self.visualizer.plot_pca(self.pca_result, self.pca, str(output_path))
        self.assertTrue(output_path.exists())
    
    def test_mds_plot(self):
        """Test MDS plot generation."""
        output_path = self.temp_dir / "test_mds.png"
        self.visualizer.plot_mds(self.mds_result, str(output_path))
        self.assertTrue(output_path.exists())


class TestModularUtils(unittest.TestCase):
    """Test the modular utility functions."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create test FASTQ files
        single_file = self.temp_dir / "sample_001.fastq"
        with open(single_file, "w") as f:
            f.write("@read1\nATCG\n+\n~~~~\n")
        
        r1_file = self.temp_dir / "sample_002_R1.fastq"
        r2_file = self.temp_dir / "sample_002_R2.fastq"
        
        with open(r1_file, "w") as f:
            f.write("@read1\nATCG\n+\n~~~~\n")
        with open(r2_file, "w") as f:
            f.write("@read1\nCGAT\n+\n~~~~\n")
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_find_fastq_files(self):
        """Test FASTQ file discovery."""
        files = find_fastq_files(str(self.temp_dir))
        
        self.assertEqual(len(files), 2)  # 1 single-end + 1 paired-end
        
        # Check file types
        single_end_found = False
        paired_end_found = False
        
        for filepath, sample_name in files:
            if isinstance(filepath, str):
                single_end_found = True
                self.assertEqual(sample_name, "sample_001")
            elif isinstance(filepath, list):
                paired_end_found = True
                self.assertEqual(sample_name, "sample_002")
                self.assertEqual(len(filepath), 2)
        
        self.assertTrue(single_end_found)
        self.assertTrue(paired_end_found)
    
    def test_validate_input_directory(self):
        """Test input directory validation."""
        # Valid directory with FASTQ files
        self.assertTrue(validate_input_directory(str(self.temp_dir)))
        
        # Non-existent directory
        self.assertFalse(validate_input_directory("/non/existent/path"))
    
    def test_summarize_fastq_files(self):
        """Test FASTQ file summarization."""
        files = find_fastq_files(str(self.temp_dir))
        summary = summarize_fastq_files(files)
        
        self.assertEqual(summary["total_samples"], 2)
        self.assertEqual(summary["paired_end_samples"], 1)
        self.assertEqual(summary["single_end_samples"], 1)
        self.assertEqual(summary["total_files"], 3)
    
    def test_check_dependencies(self):
        """Test dependency checking."""
        deps = check_dependencies()
        
        self.assertIsInstance(deps, dict)
        self.assertIn("numpy", deps)
        self.assertIn("pandas", deps)
        self.assertTrue(deps["numpy"])  # Should be available
    
    def test_get_system_info(self):
        """Test system information gathering."""
        info = get_system_info()
        
        self.assertIsInstance(info, dict)
        self.assertIn("platform", info)
        self.assertIn("python_version", info)
        self.assertIn("cpu_count", info)
        self.assertGreater(info["cpu_count"], 0)


class TestModularConfig(unittest.TestCase):
    """Test the modular configuration system."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = Path(tempfile.mkdtemp())
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_default_config(self):
        """Test default configuration creation."""
        config = MetaGrouperConfig()
        
        self.assertEqual(config.profiling.k_size, 21)
        self.assertEqual(config.profiling.min_kmer_freq, 1)
        self.assertTrue(config.profiling.memory_efficient)
        self.assertEqual(config.analysis.distance_metric, "braycurtis")
    
    def test_config_validation(self):
        """Test configuration validation."""
        config = MetaGrouperConfig()
        
        # Valid configuration should pass
        config.validate()
        
        # Invalid k-mer size should fail
        config.profiling.k_size = -1
        with self.assertRaises(ValueError):
            config.validate()
    
    def test_config_file_operations(self):
        """Test configuration file save/load."""
        config = MetaGrouperConfig()
        config.profiling.k_size = 25
        config.analysis.distance_metric = "jaccard"
        
        config_file = self.temp_dir / "test_config.json"
        
        # Save configuration
        config.save_to_file(str(config_file))
        self.assertTrue(config_file.exists())
        
        # Load configuration
        new_config = MetaGrouperConfig(str(config_file))
        self.assertEqual(new_config.profiling.k_size, 25)
        self.assertEqual(new_config.analysis.distance_metric, "jaccard")
    
    def test_config_summary(self):
        """Test configuration summary generation."""
        config = MetaGrouperConfig()
        summary = config.get_summary()
        
        self.assertIsInstance(summary, str)
        self.assertIn("MetaGrouper Configuration", summary)
        self.assertIn("K-mer size", summary)


class TestModularIntegration(unittest.TestCase):
    """Integration tests for the modular architecture."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.output_dir = self.temp_dir / "output"
        
        # Create test dataset
        self.create_test_dataset()
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def create_test_dataset(self):
        """Create a realistic test dataset."""
        sequences = [
            "ATCGATCGATCGATCGATCGATCGATCGATCG",
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
            "TTCCGGAATTCCGGAATTCCGGAATTCCGGAA"
        ]
        
        # Create 3 samples (1 single-end, 2 paired-end)
        single_file = self.temp_dir / "sample_001.fastq"
        with open(single_file, "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f"@read_{i}\n{seq}\n+\n{'~' * len(seq)}\n")
        
        for sample_num in [2, 3]:
            for direction in ["R1", "R2"]:
                filename = f"sample_00{sample_num}_{direction}.fastq"
                filepath = self.temp_dir / filename
                with open(filepath, "w") as f:
                    for i, seq in enumerate(sequences):
                        f.write(f"@read_{i}\n{seq}\n+\n{'~' * len(seq)}\n")
    
    def test_complete_workflow(self):
        """Test complete analysis workflow."""
        # Initialize configuration
        config = MetaGrouperConfig()
        config.profiling.k_size = 17
        config.profiling.max_reads = 10
        config.output.output_dir = str(self.output_dir)
        
        # Find FASTQ files
        fastq_files = find_fastq_files(str(self.temp_dir))
        self.assertEqual(len(fastq_files), 3)
        
        # Profile samples
        profiler = KmerProfiler(
            k=config.profiling.k_size,
            max_reads=config.profiling.max_reads,
            min_kmer_freq=config.profiling.min_kmer_freq
        )
        
        for filepath, sample_name in fastq_files:
            profiler.profile_sample(filepath, sample_name)
        
        self.assertEqual(len(profiler.profiles), 3)
        
        # Analyze similarities
        analyzer = SimilarityAnalyzer(profiler.profiles)
        distance_matrix = analyzer.compute_distance_matrix()
        pca_result, pca = analyzer.perform_pca()
        mds_result = analyzer.perform_mds()
        
        # Generate visualizations
        self.output_dir.mkdir(parents=True, exist_ok=True)
        visualizer = Visualizer(profiler.sample_names)
        
        visualizer.plot_distance_heatmap(
            distance_matrix, str(self.output_dir / "heatmap.png")
        )
        visualizer.plot_pca(
            pca_result, pca, str(self.output_dir / "pca.png")
        )
        visualizer.plot_mds(
            mds_result, str(self.output_dir / "mds.png")
        )
        
        # Save results
        save_results(
            profiler.profiles, distance_matrix, 
            profiler.sample_names, str(self.output_dir)
        )
        
        # Verify outputs
        expected_files = [
            "heatmap.png", "pca.png", "mds.png",
            "kmer_profiles.pkl", "distance_matrix.npy", 
            "distance_matrix.csv", "sample_names.json"
        ]
        
        for filename in expected_files:
            filepath = self.output_dir / filename
            self.assertTrue(filepath.exists(), f"Missing output file: {filename}")


def run_modular_tests():
    """Run all modular architecture tests."""
    # Set up logging for tests
    setup_logging(verbose=False)
    
    # Create test suite
    suite = unittest.TestSuite()
    
    # Add test classes
    test_classes = [
        TestModularKmerProfiler,
        TestModularSimilarityAnalyzer,
        TestModularVisualizer,
        TestModularUtils,
        TestModularConfig,
        TestModularIntegration
    ]
    
    for test_class in test_classes:
        tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()


if __name__ == "__main__":
    import sys
    
    print("🧪 Running Modular MetaGrouper Test Suite")
    print("=" * 50)
    
    success = run_modular_tests()
    
    if success:
        print("\n✅ All tests passed!")
        sys.exit(0)
    else:
        print("\n❌ Some tests failed!")
        sys.exit(1)