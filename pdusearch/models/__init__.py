"""ML model implementations for PDU feature encoding.

This module provides autoencoder architectures for learning compact representations
of PDU features, with comprehensive training utilities.

Models:
    - Autoencoder: Standard encoder-decoder for dimensionality reduction
    - VariationalAutoencoder: Probabilistic autoencoder with KL regularization

Training:
    - train_autoencoder(): High-level training function
    - train_autoencoder_ddp(): Distributed training with DDP
    - train_autoencoder_with_validation(): Training with validation set
    - save_model(), load_model(): Model persistence

Example:
    >>> import numpy as np
    >>> from pdusearch.models import Autoencoder, train_autoencoder, TrainingConfig
    >>>
    >>> # Create and train
    >>> X = np.random.randn(10000, 900).astype(np.float32)
    >>> config = TrainingConfig(latent_dim=16, epochs=10, batch_size=32)
    >>> model, Z, losses = train_autoencoder(X, config)
    >>>
    >>> # Use model
    >>> print(f"Learned {Z.shape[1]}D embeddings")
"""

from .autoencoder import Autoencoder, VariationalAutoencoder
from .training import (
    train_autoencoder,
    train_autoencoder_ddp,
    train_autoencoder_with_validation,
    save_model,
    load_model,
)
from pdusearch.config import TrainingConfig

__all__ = [
    "Autoencoder",
    "VariationalAutoencoder",
    "TrainingConfig",
    "train_autoencoder",
    "train_autoencoder_ddp",
    "train_autoencoder_with_validation",
    "save_model",
    "load_model",
]
