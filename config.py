#!/usr/bin/env python3
"""
MetaGrouper Configuration Management

Handles configuration files, parameter presets, and performance optimization.
"""

import json
import yaml
from pathlib import Path
from typing import Dict, Any, Optional, List
import logging


class MetaGrouperConfig:
    """Configuration manager for MetaGrouper."""

    def __init__(self, config_file: Optional[str] = None):
        """Initialize configuration with default values."""

        # Default configuration
        self.config = {
            "kmer_analysis": {
                "default_k": 21,
                "min_k": 11,
                "max_k": 31,
                "default_distance_metric": "braycurtis",
                "supported_metrics": ["braycurtis", "jaccard", "cosine", "euclidean"],
            },
            "metadata_analysis": {
                "default_permutations": 999,
                "min_permutations": 99,
                "max_permutations": 99999,
                "significance_threshold": 0.05,
                "cluster_range": [2, 8],
            },
            "assembly_recommendations": {
                "similarity_thresholds": {
                    "stringent": 0.15,
                    "medium": 0.25,
                    "permissive": 0.40,
                },
                "default_similarity_threshold": 0.30,
                "min_group_size": 2,
                "max_group_size": 10,
                "supported_tools": ["megahit", "spades", "flye"],
            },
            "performance": {
                "memory_limit_gb": 8,
                "max_samples_in_memory": 100,
                "default_max_reads": None,
                "chunk_size": 1000,
            },
            "output": {
                "save_intermediate": True,
                "compression_level": 6,
                "plot_dpi": 300,
                "plot_format": "png",
            },
            "presets": {
                "quick": {
                    "kmer_size": 15,
                    "max_reads": 1000,
                    "permutations": 99,
                    "similarity_threshold": 0.35,
                },
                "default": {
                    "kmer_size": 21,
                    "max_reads": None,
                    "permutations": 999,
                    "similarity_threshold": 0.30,
                },
                "high_quality": {
                    "kmer_size": 25,
                    "max_reads": None,
                    "permutations": 9999,
                    "similarity_threshold": 0.20,
                },
                "large_dataset": {
                    "kmer_size": 17,
                    "max_reads": 5000,
                    "permutations": 499,
                    "similarity_threshold": 0.25,
                },
            },
        }

        # Load config file if provided
        if config_file:
            self.load_config(config_file)

    def load_config(self, config_file: str):
        """Load configuration from file."""
        config_path = Path(config_file)

        if not config_path.exists():
            logging.warning(f"Config file not found: {config_file}")
            return

        try:
            if config_file.endswith(".json"):
                with open(config_path, "r") as f:
                    user_config = json.load(f)
            elif config_file.endswith((".yml", ".yaml")):
                with open(config_path, "r") as f:
                    user_config = yaml.safe_load(f)
            else:
                logging.error(f"Unsupported config format: {config_file}")
                return

            # Merge user config with defaults
            self._merge_config(self.config, user_config)
            logging.info(f"Loaded configuration from {config_file}")

        except Exception as e:
            logging.error(f"Failed to load config file {config_file}: {e}")

    def _merge_config(self, default: Dict, user: Dict):
        """Recursively merge user config with defaults."""
        for key, value in user.items():
            if (
                key in default
                and isinstance(default[key], dict)
                and isinstance(value, dict)
            ):
                self._merge_config(default[key], value)
            else:
                default[key] = value

    def save_config(self, config_file: str):
        """Save current configuration to file."""
        config_path = Path(config_file)

        try:
            if config_file.endswith(".json"):
                with open(config_path, "w") as f:
                    json.dump(self.config, f, indent=2)
            elif config_file.endswith((".yml", ".yaml")):
                with open(config_path, "w") as f:
                    yaml.dump(self.config, f, indent=2)
            else:
                logging.error(f"Unsupported config format: {config_file}")
                return

            logging.info(f"Configuration saved to {config_file}")

        except Exception as e:
            logging.error(f"Failed to save config file {config_file}: {e}")

    def get_preset(self, preset_name: str) -> Dict[str, Any]:
        """Get configuration preset."""
        if preset_name not in self.config["presets"]:
            available = list(self.config["presets"].keys())
            raise ValueError(f"Unknown preset: {preset_name}. Available: {available}")

        return self.config["presets"][preset_name].copy()

    def get_performance_config(self, dataset_size: str = "medium") -> Dict[str, Any]:
        """Get performance configuration based on dataset size."""

        base_config = self.config["performance"].copy()

        if dataset_size == "small":
            # < 10 samples
            base_config.update(
                {"max_reads": None, "permutations": 999, "kmer_size": 21}
            )
        elif dataset_size == "medium":
            # 10-50 samples
            base_config.update(
                {"max_reads": 10000, "permutations": 999, "kmer_size": 19}
            )
        elif dataset_size == "large":
            # 50-200 samples
            base_config.update(
                {"max_reads": 5000, "permutations": 499, "kmer_size": 17}
            )
        elif dataset_size == "xlarge":
            # > 200 samples
            base_config.update(
                {"max_reads": 2000, "permutations": 199, "kmer_size": 15}
            )

        return base_config

    def validate_parameters(self, params: Dict[str, Any]) -> List[str]:
        """Validate parameters against configuration constraints."""
        errors = []

        # Validate k-mer size
        if "kmer_size" in params:
            k = params["kmer_size"]
            min_k = self.config["kmer_analysis"]["min_k"]
            max_k = self.config["kmer_analysis"]["max_k"]

            if k < min_k or k > max_k:
                errors.append(f"K-mer size must be between {min_k} and {max_k}")

            if k % 2 == 0:
                errors.append("K-mer size must be odd")

        # Validate distance metric
        if "distance_metric" in params:
            metric = params["distance_metric"]
            supported = self.config["kmer_analysis"]["supported_metrics"]
            if metric not in supported:
                errors.append(f"Distance metric must be one of: {supported}")

        # Validate permutations
        if "permutations" in params:
            perms = params["permutations"]
            min_perms = self.config["metadata_analysis"]["min_permutations"]
            max_perms = self.config["metadata_analysis"]["max_permutations"]

            if perms < min_perms or perms > max_perms:
                errors.append(
                    f"Permutations must be between {min_perms} and {max_perms}"
                )

        # Validate similarity threshold
        if "similarity_threshold" in params:
            threshold = params["similarity_threshold"]
            if not 0.0 <= threshold <= 1.0:
                errors.append("Similarity threshold must be between 0.0 and 1.0")

        # Validate group sizes
        if "min_group_size" in params and "max_group_size" in params:
            min_size = params["min_group_size"]
            max_size = params["max_group_size"]

            if min_size >= max_size:
                errors.append("Minimum group size must be less than maximum group size")

        return errors

    def optimize_for_system(self, memory_gb: float, cpu_cores: int) -> Dict[str, Any]:
        """Optimize configuration for system specifications."""

        optimized = {}

        # Memory-based optimizations
        if memory_gb < 4:
            optimized.update({"kmer_size": 15, "max_reads": 1000, "permutations": 99})
        elif memory_gb < 8:
            optimized.update({"kmer_size": 17, "max_reads": 5000, "permutations": 499})
        elif memory_gb < 16:
            optimized.update({"kmer_size": 19, "max_reads": 10000, "permutations": 999})
        else:
            optimized.update({"kmer_size": 21, "max_reads": None, "permutations": 999})

        # CPU-based optimizations
        if cpu_cores >= 8:
            optimized["permutations"] = min(
                optimized.get("permutations", 999) * 2, 9999
            )

        return optimized

    def get_assembly_tool_config(self, tool: str) -> Dict[str, Any]:
        """Get tool-specific assembly configuration."""

        configs = {
            "megahit": {
                "min_contig_len": 500,
                "k_list": [21, 29, 39, 59, 79, 99],
                "memory_flag": "--memory",
                "threads_flag": "--num-cpu-threads",
            },
            "spades": {
                "mode": "--meta",
                "memory_flag": "--memory",
                "threads_flag": "--threads",
                "careful": True,
            },
            "flye": {
                "mode": "--meta",
                "read_type": "--nano-raw",  # or --nano-corr, --pacbio-raw, etc.
                "threads_flag": "--threads",
            },
        }

        if tool not in configs:
            raise ValueError(f"Unknown assembly tool: {tool}")

        return configs[tool].copy()


def create_default_config_file(output_path: str):
    """Create a default configuration file."""

    config = MetaGrouperConfig()
    config.save_config(output_path)

    print(f"Default configuration saved to: {output_path}")
    print("Edit this file to customize MetaGrouper behavior.")


def estimate_resource_requirements(
    n_samples: int, avg_reads_per_sample: int, kmer_size: int
) -> Dict[str, Any]:
    """Estimate computational resource requirements."""

    # Rough estimates based on typical performance

    # Memory estimation
    unique_kmers_per_sample = min(
        avg_reads_per_sample * (150 - kmer_size + 1), 4**kmer_size
    )  # Theoretical max

    memory_per_sample_mb = (
        unique_kmers_per_sample * 8 / 1024 / 1024
    )  # 8 bytes per k-mer
    total_memory_mb = memory_per_sample_mb * n_samples * 1.5  # 1.5x overhead

    # Runtime estimation
    base_time_per_sample = 30  # seconds
    kmer_factor = (kmer_size / 21) ** 2
    reads_factor = avg_reads_per_sample / 10000

    estimated_runtime_minutes = (
        base_time_per_sample * n_samples * kmer_factor * reads_factor
    ) / 60

    return {
        "estimated_memory_gb": total_memory_mb / 1024,
        "estimated_runtime_minutes": estimated_runtime_minutes,
        "estimated_storage_gb": n_samples * 0.1,  # Conservative estimate
        "recommended_preset": (
            "quick"
            if n_samples > 100 or total_memory_mb > 8192
            else "high_quality" if n_samples < 20 else "default"
        ),
    }


def get_system_info() -> Dict[str, Any]:
    """Get system information for optimization."""
    import psutil
    import multiprocessing

    return {
        "memory_gb": psutil.virtual_memory().total / (1024**3),
        "cpu_cores": multiprocessing.cpu_count(),
        "available_memory_gb": psutil.virtual_memory().available / (1024**3),
        "disk_space_gb": psutil.disk_usage(".").free / (1024**3),
    }


def main():
    """Main function for config command-line interface."""
    import argparse

    parser = argparse.ArgumentParser(description="MetaGrouper Configuration Tools")
    parser.add_argument("--create-config", help="Create default config file")
    parser.add_argument(
        "--system-info", action="store_true", help="Show system information"
    )
    parser.add_argument(
        "--estimate-resources",
        nargs=3,
        type=int,
        metavar=("samples", "reads_per_sample", "kmer_size"),
        help="Estimate resource requirements",
    )

    args = parser.parse_args()

    if args.create_config:
        create_default_config_file(args.create_config)

    if args.system_info:
        info = get_system_info()
        print("System Information:")
        for key, value in info.items():
            print(
                f"  {key}: {value:.1f}"
                if isinstance(value, float)
                else f"  {key}: {value}"
            )

    if args.estimate_resources:
        n_samples, reads_per_sample, kmer_size = args.estimate_resources
        estimates = estimate_resource_requirements(
            n_samples, reads_per_sample, kmer_size
        )

        print("Resource Estimates:")
        for key, value in estimates.items():
            if isinstance(value, float):
                print(f"  {key}: {value:.1f}")
            else:
                print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
