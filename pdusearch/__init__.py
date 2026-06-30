"""PDU Search: Structure-driven rational protein design framework.

A modular library for extracting, analyzing, and clustering protein dynamic units (PDUs)
from protein structures to enable rational protein design.

Example usage:
    from pdusearch.features import FeatureExtractor, Encoding
    from pdusearch.models import train_autoencoder, TrainingConfig
    from pdusearch.clustering import HDBSCANClusterer
    from pdusearch.config import Config

    # Extract features
    config = Config(db="per_aa_sqlite/pdus_L.sqlite", aa="L")
    extractor = FeatureExtractor(config, encoding=Encoding.AA_20)
    X, pdu_ids = extractor.extract_for_aa("L")

    # Train model
    training_cfg = TrainingConfig(latent_dim=16, epochs=10)
    model, Z, losses = train_autoencoder(X, training_cfg)

    # Cluster
    clusterer = HDBSCANClusterer()
    labels, confidences = clusterer.cluster(Z)
"""

__version__ = "0.1.0"
__author__ = "Sarath Dantu"
__email__ = "sarathchandradantu@gmail.com"

# Core imports
from .config import (
    Config,
    TrainingConfig,
    ClusteringConfig,
    EnrichmentConfig,
    ValidationConfig,
    FeatureConfig,
    ConfigFactory,
)

from .logging_utils import (
    configure_logging,
    get_logger,
    add_logging_args,
    log_timestamp,
    log_separator,
)

from .utils import (
    AA_ORDER,
    AA_TO_INDEX,
    SS_ORDER,
    SS_TO_INDEX,
    FunctionalClass,
    AA_TO_FUNCTIONAL,
)

from .io import (
    PDUDatabase,
    BatchPDUFetcher,
)

__all__ = [
    # Config
    "Config",
    "TrainingConfig",
    "ClusteringConfig",
    "EnrichmentConfig",
    "ValidationConfig",
    "FeatureConfig",
    "ConfigFactory",
    # Logging
    "configure_logging",
    "get_logger",
    "add_logging_args",
    "log_timestamp",
    "log_separator",
    # Utils
    "AA_ORDER",
    "AA_TO_INDEX",
    "SS_ORDER",
    "SS_TO_INDEX",
    "FunctionalClass",
    "AA_TO_FUNCTIONAL",
    # I/O
    "PDUDatabase",
    "BatchPDUFetcher",
]
