"""Configuration classes for PDU extraction, training, and clustering."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List


@dataclass
class Config:
    """Base configuration for PDU processing."""

    # Database
    db: str = "per_aa_sqlite/pdus_L.sqlite"

    # Directories
    out_dir: str = "analysis"
    data_dir: str = "data"

    # Feature extraction
    radius: float = 15.0
    bin_width: float = 1.0
    min_pdus: int = 25

    # Amino acid to process
    aa: Optional[str] = None

    # Logging
    log_level: str = "INFO"
    log_file: Optional[str] = None

    def __post_init__(self):
        """Validate and convert paths."""
        self.out_dir = str(Path(self.out_dir))
        self.data_dir = str(Path(self.data_dir))

    @property
    def n_distance_bins(self) -> int:
        """Calculate number of distance bins."""
        import numpy as np
        return len(np.arange(0.0, self.radius + self.bin_width, self.bin_width)) - 1


@dataclass
class TrainingConfig:
    """Configuration for autoencoder training."""

    # Model architecture
    latent_dim: int = 16
    hidden_dim_1: Optional[int] = None
    hidden_dim_2: Optional[int] = None

    # Training parameters
    epochs: int = 10
    batch_size: int = 32
    learning_rate: float = 1e-3

    # Optimization
    optimizer: str = "adam"  # adam, sgd, rmsprop
    loss_function: str = "mse"  # mse, bce

    # Hardware
    device: str = "auto"  # auto, cpu, cuda, mps

    # Early stopping (optional)
    early_stopping_patience: Optional[int] = None
    validation_split: float = 0.0

    # Checkpointing
    save_checkpoints: bool = False
    checkpoint_dir: Optional[str] = None

    def __post_init__(self):
        """Resolve device and validate parameters."""
        if self.device == "auto":
            self.device = self._resolve_device()

        # Validate parameters
        if self.latent_dim <= 0:
            raise ValueError(f"latent_dim must be > 0, got {self.latent_dim}")
        if self.epochs <= 0:
            raise ValueError(f"epochs must be > 0, got {self.epochs}")
        if self.batch_size <= 0:
            raise ValueError(f"batch_size must be > 0, got {self.batch_size}")
        if self.learning_rate <= 0:
            raise ValueError(f"learning_rate must be > 0, got {self.learning_rate}")
        if not 0 <= self.validation_split < 1:
            raise ValueError(f"validation_split must be in [0,1), got {self.validation_split}")

    @staticmethod
    def _resolve_device() -> str:
        """Resolve best available device."""
        try:
            import torch
            if torch.cuda.is_available():
                return "cuda"
            elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
                return "mps"
        except ImportError:
            pass
        return "cpu"

    def get_device(self):
        """Get torch device object."""
        import torch
        return torch.device(self.device)


@dataclass
class ClusteringConfig:
    """Configuration for HDBSCAN clustering."""

    # HDBSCAN parameters
    min_cluster_size: int = 200
    min_samples: int = 5
    cluster_selection_epsilon: float = 0.0
    cluster_selection_method: str = "eom"  # eom, leaf

    # Metrics
    metric: str = "euclidean"  # euclidean, manhattan, cosine, etc.

    # Post-processing
    allow_single_cluster: bool = False
    filter_orphans: bool = True

    def __post_init__(self):
        """Validate parameters."""
        if self.min_cluster_size <= 0:
            raise ValueError(f"min_cluster_size must be > 0, got {self.min_cluster_size}")
        if self.min_samples <= 0:
            raise ValueError(f"min_samples must be > 0, got {self.min_samples}")
        if self.min_samples > self.min_cluster_size:
            raise ValueError(
                f"min_samples ({self.min_samples}) cannot exceed "
                f"min_cluster_size ({self.min_cluster_size})"
            )


@dataclass
class EnrichmentConfig:
    """Configuration for enrichment analysis."""

    # PyCoM
    pycom_db: str = "pycom_db/pycom.db"

    # Statistical tests
    test_method: str = "chi2"  # chi2, fisher
    min_category_size: int = 5
    significance_level: float = 0.05

    # Caching
    cache_alignments: bool = True
    cache_dir: str = ".cache/pycom"


@dataclass
class ValidationConfig:
    """Configuration for validation pipelines."""

    # Grid search
    latent_dims: List[int] = field(default_factory=lambda: [8, 16, 32, 64, 128, 256])
    min_cluster_sizes: List[int] = field(default_factory=lambda: [50, 100, 200, 300, 500])

    # Bootstrap/perturbation
    n_bootstrap_iterations: int = 10
    perturbation_fraction: float = 0.1  # Remove 10% of points
    n_perturbation_iterations: int = 10

    # Metrics to track
    metrics: List[str] = field(
        default_factory=lambda: ["noise_percent", "avg_confidence", "n_clusters", "largest_cluster"]
    )


@dataclass
class FeatureConfig:
    """Configuration for feature extraction and encoding."""

    # Feature type
    encoding: str = "aa_20"  # aa_20, functional_5, aa_20_conservation

    # Conservation metrics (if using conservation encoding)
    compute_conservation: bool = False
    conservation_source: str = "pycom"  # pycom, uniprot, msa

    # Distance shells
    distance_bins: Optional[List[float]] = None

    # Feature normalization
    normalize: bool = False
    normalization_method: str = "zscore"  # zscore, minmax, none


class ConfigFactory:
    """Factory for creating configurations."""

    @staticmethod
    def default_extraction() -> Config:
        """Default configuration for feature extraction."""
        return Config()

    @staticmethod
    def default_training() -> TrainingConfig:
        """Default configuration for training."""
        return TrainingConfig()

    @staticmethod
    def default_clustering() -> ClusteringConfig:
        """Default configuration for clustering."""
        return ClusteringConfig()

    @staticmethod
    def default_enrichment() -> EnrichmentConfig:
        """Default configuration for enrichment."""
        return EnrichmentConfig()

    @staticmethod
    def default_validation() -> ValidationConfig:
        """Default configuration for validation."""
        return ValidationConfig()

    @staticmethod
    def compression_hypothesis() -> tuple:
        """Configuration for compression hypothesis test."""
        feature_config = FeatureConfig(encoding="aa_20")
        training_config = TrainingConfig(epochs=10, batch_size=32)
        clustering_config = ClusteringConfig(min_cluster_size=200, min_samples=5)
        return feature_config, training_config, clustering_config
