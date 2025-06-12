#!/usr/bin/env python3
"""
Test suite for MetaGrouper

Comprehensive tests for all three phases of MetaGrouper functionality.
"""

import unittest
import tempfile
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
import os
import sys

# Add the current directory to the path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from metagrouper import KmerProfiler, SimilarityAnalyzer, Visualizer
from metadata_analyzer import MetadataAnalyzer, PermanovaAnalyzer
from assembly_recommender import AssemblyRecommender, AssemblyStrategyEngine


class TestKmerProfiler(unittest.TestCase):
    """Test the KmerProfiler class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.profiler = KmerProfiler(k=15, max_reads=100)

        # Create test FASTQ file
        self.test_fastq = Path(self.temp_dir) / "test_sample.fastq"
        with open(self.test_fastq, "w") as f:
            for i in range(10):
                f.write(f"@read_{i}\n")
                f.write("ATCGATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("~~~~~~~~~~~~~~~~~~~~\n")

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)

    def test_reverse_complement(self):
        """Test reverse complement function."""
        result = self.profiler._reverse_complement("ATCG")
        self.assertEqual(result, "CGAT")

        result = self.profiler._reverse_complement("AAAA")
        self.assertEqual(result, "TTTT")

    def test_extract_kmers(self):
        """Test k-mer extraction."""
        sequence = "ATCGATCGATCG"
        kmers = self.profiler._extract_kmers(sequence)

        # Should extract k-mers of length 15, but sequence is only 12 bp
        # So no k-mers should be extracted
        self.assertEqual(len(kmers), 0)

        # Test with longer sequence
        long_sequence = "ATCGATCGATCGATCGATCGATCG"
        kmers = self.profiler._extract_kmers(long_sequence)
        self.assertGreater(len(kmers), 0)

    def test_profile_sample(self):
        """Test sample profiling."""
        profile = self.profiler.profile_sample(str(self.test_fastq), "test_sample")

        self.assertIsInstance(profile, dict)
        self.assertIn("test_sample", self.profiler.sample_names)
        self.assertIn("test_sample", self.profiler.profiles)

        # Check that profile values sum to 1 (normalized)
        total = sum(profile.values())
        self.assertAlmostEqual(total, 1.0, places=5)


class TestSimilarityAnalyzer(unittest.TestCase):
    """Test the SimilarityAnalyzer class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock k-mer profiles
        self.profiles = {
            "sample_A": {"ATG": 0.3, "TGC": 0.3, "GCA": 0.4},
            "sample_B": {"ATG": 0.2, "TGC": 0.4, "GCA": 0.4},
            "sample_C": {"ATG": 0.1, "TGC": 0.1, "XXX": 0.8},  # Very different
        }
        self.analyzer = SimilarityAnalyzer(self.profiles)

    def test_compute_distance_matrix(self):
        """Test distance matrix computation."""
        distance_matrix = self.analyzer.compute_distance_matrix()

        # Check shape
        self.assertEqual(distance_matrix.shape, (3, 3))

        # Check diagonal is zero
        np.testing.assert_array_almost_equal(np.diag(distance_matrix), [0, 0, 0])

        # Check symmetry
        np.testing.assert_array_almost_equal(distance_matrix, distance_matrix.T)

        # Check that sample_C is most distant from others
        self.assertGreater(distance_matrix[0, 2], distance_matrix[0, 1])
        self.assertGreater(distance_matrix[1, 2], distance_matrix[0, 1])

    def test_perform_pca(self):
        """Test PCA analysis."""
        pca_result, pca = self.analyzer.perform_pca()

        # Check output shape
        self.assertEqual(pca_result.shape, (3, 2))

        # Check that PCA object is returned
        self.assertIsNotNone(pca)
        self.assertEqual(len(pca.explained_variance_ratio_), 2)

    def test_perform_mds(self):
        """Test MDS analysis."""
        self.analyzer.compute_distance_matrix()  # Required for MDS
        mds_result = self.analyzer.perform_mds()

        # Check output shape
        self.assertEqual(mds_result.shape, (3, 2))


class TestPermanovaAnalyzer(unittest.TestCase):
    """Test the PermanovaAnalyzer class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock distance matrix
        self.distance_matrix = np.array(
            [
                [0.0, 0.5, 0.8, 0.9],
                [0.5, 0.0, 0.7, 0.8],
                [0.8, 0.7, 0.0, 0.3],
                [0.9, 0.8, 0.3, 0.0],
            ]
        )
        self.sample_names = ["A", "B", "C", "D"]
        self.analyzer = PermanovaAnalyzer(self.distance_matrix, self.sample_names)

    def test_permanova_test(self):
        """Test PERMANOVA statistical test."""
        # Test with simple grouping variable
        groups = np.array([0, 0, 1, 1])  # Two groups

        result = self.analyzer.permanova_test(groups, n_permutations=99)

        # Check that all expected keys are present
        expected_keys = ["f_statistic", "p_value", "r_squared", "n_samples", "n_groups"]
        for key in expected_keys:
            self.assertIn(key, result)

        # Check that values are reasonable
        self.assertGreaterEqual(result["r_squared"], 0)
        self.assertLessEqual(result["r_squared"], 1)
        self.assertGreaterEqual(result["p_value"], 0)
        self.assertLessEqual(result["p_value"], 1)
        self.assertEqual(result["n_samples"], 4)
        self.assertEqual(result["n_groups"], 2)


class TestMetadataAnalyzer(unittest.TestCase):
    """Test the MetadataAnalyzer class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

        # Create mock distance matrix and sample names
        self.distance_matrix = np.array(
            [
                [0.0, 0.3, 0.7, 0.8],
                [0.3, 0.0, 0.6, 0.7],
                [0.7, 0.6, 0.0, 0.2],
                [0.8, 0.7, 0.2, 0.0],
            ]
        )
        self.sample_names = ["sample_A", "sample_B", "sample_C", "sample_D"]

        # Create test metadata file
        self.metadata_file = Path(self.temp_dir) / "metadata.csv"
        metadata_data = {
            "sample_id": ["sample_A", "sample_B", "sample_C", "sample_D"],
            "group": ["control", "control", "treatment", "treatment"],
            "site": ["A", "A", "B", "B"],
            "age": [25, 30, 35, 40],
        }
        pd.DataFrame(metadata_data).to_csv(self.metadata_file, index=False)

        self.analyzer = MetadataAnalyzer(self.distance_matrix, self.sample_names)

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)

    def test_load_metadata(self):
        """Test metadata loading."""
        self.analyzer.load_metadata(str(self.metadata_file))

        self.assertIsNotNone(self.analyzer.metadata)
        self.assertEqual(len(self.analyzer.metadata), 4)
        self.assertIn("group", self.analyzer.metadata.columns)
        self.assertIn("site", self.analyzer.metadata.columns)
        self.assertIn("age", self.analyzer.metadata.columns)

    def test_analyze_variables(self):
        """Test variable analysis with PERMANOVA."""
        self.analyzer.load_metadata(str(self.metadata_file))
        results = self.analyzer.analyze_variables(n_permutations=99)

        self.assertIsInstance(results, pd.DataFrame)
        self.assertGreater(len(results), 0)

        # Check expected columns
        expected_columns = ["variable", "r_squared", "p_value", "variable_type"]
        for col in expected_columns:
            self.assertIn(col, results.columns)

    def test_identify_clusters(self):
        """Test cluster identification."""
        cluster_results = self.analyzer.identify_clusters(n_clusters_range=(2, 4))

        self.assertIn("kmeans", cluster_results)
        self.assertIn("hierarchical", cluster_results)

        # Check that optimal clusters are identified
        for method in ["kmeans", "hierarchical"]:
            if "optimal" in cluster_results[method]:
                optimal = cluster_results[method]["optimal"]
                self.assertIn("n_clusters", optimal)
                self.assertIn("silhouette_score", optimal)


class TestAssemblyStrategyEngine(unittest.TestCase):
    """Test the AssemblyStrategyEngine class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock distance matrix with clear groups
        self.distance_matrix = np.array(
            [
                [0.0, 0.1, 0.8, 0.9],  # A and B are similar
                [0.1, 0.0, 0.9, 0.8],
                [0.8, 0.9, 0.0, 0.1],  # C and D are similar
                [0.9, 0.8, 0.1, 0.0],
            ]
        )
        self.sample_names = ["sample_A", "sample_B", "sample_C", "sample_D"]
        self.engine = AssemblyStrategyEngine(self.distance_matrix, self.sample_names)

    def test_assess_group_quality(self):
        """Test group quality assessment."""
        # Test high-quality group (similar samples)
        high_quality_group = [0, 1]  # A and B
        quality_score = self.engine._assess_group_quality(high_quality_group)
        self.assertGreater(quality_score, 0.5)

        # Test low-quality group (dissimilar samples)
        low_quality_group = [0, 2]  # A and C
        quality_score = self.engine._assess_group_quality(low_quality_group)
        self.assertLess(quality_score, 0.5)

    def test_recommend_by_similarity(self):
        """Test similarity-based recommendations."""
        groups = self.engine.recommend_by_similarity()

        self.assertIsInstance(groups, list)

        # Should find groups for similar samples
        if groups:
            for group in groups:
                self.assertGreater(len(group.sample_names), 1)
                self.assertLessEqual(
                    group.avg_distance, self.engine.similarity_threshold_medium
                )


class TestAssemblyRecommender(unittest.TestCase):
    """Test the complete AssemblyRecommender class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

        # Create mock data
        self.distance_matrix = np.array(
            [
                [0.0, 0.2, 0.8, 0.7],
                [0.2, 0.0, 0.7, 0.8],
                [0.8, 0.7, 0.0, 0.3],
                [0.7, 0.8, 0.3, 0.0],
            ]
        )
        self.sample_names = ["sample_A", "sample_B", "sample_C", "sample_D"]

        # Create metadata
        metadata_data = {
            "sample_id": self.sample_names,
            "group": ["X", "X", "Y", "Y"],
            "site": ["1", "1", "2", "2"],
        }
        self.metadata = pd.DataFrame(metadata_data).set_index("sample_id")

        # Create mock PERMANOVA results
        self.metadata_results = pd.DataFrame(
            {
                "variable": ["group", "site"],
                "r_squared": [0.6, 0.4],
                "p_value": [0.01, 0.05],
                "variable_type": ["categorical", "categorical"],
            }
        )

        self.recommender = AssemblyRecommender(self.distance_matrix, self.sample_names)

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)

    def test_generate_recommendations(self):
        """Test comprehensive recommendation generation."""
        recommendation = self.recommender.generate_recommendations(
            metadata_results=self.metadata_results, metadata=self.metadata
        )

        # Check that recommendation object has expected attributes
        self.assertIsNotNone(recommendation.strategy)
        self.assertIn(recommendation.strategy, ["individual", "grouped", "global"])

        self.assertIsInstance(recommendation.groups, list)
        self.assertIsInstance(recommendation.overall_confidence, (int, float))
        self.assertIsInstance(recommendation.assembly_commands, dict)

        # Check assembly commands
        self.assertIn("megahit", recommendation.assembly_commands)
        self.assertIn("spades", recommendation.assembly_commands)

        # Check performance predictions
        self.assertIn("strategy_summary", recommendation.performance_predictions)


class TestIntegration(unittest.TestCase):
    """Integration tests for the complete workflow."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

        # Create test FASTQ files
        self.fastq_dir = Path(self.temp_dir) / "fastq_files"
        self.fastq_dir.mkdir()

        for sample in ["sample_A", "sample_B", "sample_C"]:
            fastq_file = self.fastq_dir / f"{sample}.fastq"
            with open(fastq_file, "w") as f:
                for i in range(5):  # Small files for testing
                    f.write(f"@{sample}_read_{i}\n")
                    f.write("ATCGATCGATCGATCGATCGATCG\n")
                    f.write("+\n")
                    f.write("~~~~~~~~~~~~~~~~~~~~~~~~\n")

        # Create metadata file
        self.metadata_file = Path(self.temp_dir) / "metadata.csv"
        metadata_data = {
            "sample_id": ["sample_A", "sample_B", "sample_C"],
            "group": ["control", "control", "treatment"],
            "site": ["A", "A", "B"],
        }
        pd.DataFrame(metadata_data).to_csv(self.metadata_file, index=False)

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)

    def test_end_to_end_workflow(self):
        """Test complete end-to-end workflow."""
        # This would test the complete pipeline
        # For now, just test that components can be integrated

        # Phase 1: K-mer profiling
        profiler = KmerProfiler(k=15, max_reads=10)

        fastq_files = list(self.fastq_dir.glob("*.fastq"))
        for fastq_file in fastq_files:
            sample_name = fastq_file.stem
            profiler.profile_sample(str(fastq_file), sample_name)

        self.assertEqual(len(profiler.profiles), 3)

        # Phase 2: Similarity analysis
        analyzer = SimilarityAnalyzer(profiler.profiles)
        distance_matrix = analyzer.compute_distance_matrix()

        self.assertEqual(distance_matrix.shape, (3, 3))

        # Phase 3: Metadata analysis
        meta_analyzer = MetadataAnalyzer(distance_matrix, profiler.sample_names)
        meta_analyzer.load_metadata(str(self.metadata_file))
        results = meta_analyzer.analyze_variables(n_permutations=99)

        self.assertGreater(len(results), 0)

        # Phase 4: Assembly recommendations
        recommender = AssemblyRecommender(distance_matrix, profiler.sample_names)
        recommendation = recommender.generate_recommendations(
            metadata_results=results, metadata=meta_analyzer.metadata
        )

        self.assertIsNotNone(recommendation.strategy)


def run_tests():
    """Run all tests."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add test classes
    test_classes = [
        TestKmerProfiler,
        TestSimilarityAnalyzer,
        TestPermanovaAnalyzer,
        TestMetadataAnalyzer,
        TestAssemblyStrategyEngine,
        TestAssemblyRecommender,
        TestIntegration,
    ]

    for test_class in test_classes:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Return success status
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
