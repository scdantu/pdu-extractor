"""Autoencoder architecture for PDU feature compression.

This module defines neural network autoencoders for compressing high-dimensional PDU
features into lower-dimensional latent representations. Autoencoders learn a compact
embedding of neighborhood structure suitable for downstream clustering.

Architecture:
    Symmetric encoder-decoder with configurable hidden layer dimensions.

    Encoder:  input_dim → hidden_1 → hidden_2 → latent_dim
    Decoder:  latent_dim → hidden_2 → hidden_1 → input_dim

Example:
    >>> import torch
    >>> from pdusearch.models import Autoencoder
    >>>
    >>> model = Autoencoder(input_dim=900, latent_dim=16)
    >>> x = torch.randn(32, 900)  # Batch of 32 PDUs
    >>> x_reconstructed = model(x)  # Encode + decode
    >>> latent = model.encode(x)  # Just embeddings
    >>> print(latent.shape)  # (32, 16)
"""

import torch
import torch.nn as nn
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class Autoencoder(nn.Module):
    """Variational autoencoder for PDU feature compression.

    Learns to encode high-dimensional PDU feature vectors into a compact latent
    representation and reconstruct the original features from this encoding.
    Used for dimensionality reduction before clustering.

    Architecture:
        - Encoder: input_dim → hidden_1 → hidden_2 → latent_dim
        - Decoder: latent_dim → hidden_2 → hidden_1 → output_dim
        - Activation: ReLU in hidden layers, Sigmoid in output (bounded to [0,1])

    The hidden layer dimensions are automatically calculated based on input size
    to create a smooth bottleneck:
        - hidden_1 = min(512, max(64, input_dim // 2))
        - hidden_2 = min(256, max(32, input_dim // 4))

    Attributes:
        encoder: Encoder network
        decoder: Decoder network
        input_dim: Input feature dimension
        latent_dim: Latent embedding dimension
    """

    def __init__(
        self,
        input_dim: int,
        latent_dim: int,
        hidden_dim_1: Optional[int] = None,
        hidden_dim_2: Optional[int] = None,
    ):
        """Initialize autoencoder architecture.

        Args:
            input_dim: Input feature dimension (e.g., 900 for AA20, 225 for FUNCTIONAL_5)
            latent_dim: Target latent embedding dimension (e.g., 16, 32, 64)
            hidden_dim_1: First hidden layer size. If None, auto-calculated as
                         min(512, max(64, input_dim // 2))
            hidden_dim_2: Second hidden layer size. If None, auto-calculated as
                         min(256, max(32, input_dim // 4))

        Raises:
            ValueError: If dimensions are invalid (<=0 or latent >= input)
        """
        super().__init__()

        # Validate dimensions
        if input_dim <= 0:
            raise ValueError(f"input_dim must be > 0, got {input_dim}")
        if latent_dim <= 0:
            raise ValueError(f"latent_dim must be > 0, got {latent_dim}")
        if latent_dim >= input_dim:
            raise ValueError(
                f"latent_dim ({latent_dim}) must be < input_dim ({input_dim}) "
                "for meaningful compression"
            )

        self.input_dim = input_dim
        self.latent_dim = latent_dim

        # Auto-calculate hidden dimensions if not provided
        if hidden_dim_1 is None:
            hidden_dim_1 = min(512, max(64, input_dim // 2))
        if hidden_dim_2 is None:
            hidden_dim_2 = min(256, max(32, input_dim // 4))

        self.hidden_dim_1 = hidden_dim_1
        self.hidden_dim_2 = hidden_dim_2

        logger.debug(
            f"Autoencoder: {input_dim}D → {hidden_dim_1}D → {hidden_dim_2}D → "
            f"{latent_dim}D → {hidden_dim_2}D → {hidden_dim_1}D → {input_dim}D"
        )

        # Encoder: progressively compress features
        # Each layer reduces dimensionality by ~50% with ReLU activation
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim_1),
            nn.ReLU(),
            nn.Linear(hidden_dim_1, hidden_dim_2),
            nn.ReLU(),
            nn.Linear(hidden_dim_2, latent_dim),
            # No activation on latent layer - allows unbounded embeddings
        )

        # Decoder: progressively expand back to original dimension
        # Symmetric to encoder but with Sigmoid output to bound to [0,1]
        # (PDU features are counts, naturally in [0, ∞) but typically small)
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim_2),
            nn.ReLU(),
            nn.Linear(hidden_dim_2, hidden_dim_1),
            nn.ReLU(),
            nn.Linear(hidden_dim_1, input_dim),
            nn.Sigmoid(),  # Bound output to [0, 1]
        )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Encode features to latent representation.

        Args:
            x: Input tensor of shape (batch_size, input_dim)

        Returns:
            Latent embeddings of shape (batch_size, latent_dim)

        Example:
            >>> model = Autoencoder(900, 16)
            >>> x = torch.randn(32, 900)
            >>> z = model.encode(x)  # (32, 16)
        """
        return self.encoder(x)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """Decode latent representation back to feature space.

        Args:
            z: Latent embeddings of shape (batch_size, latent_dim)

        Returns:
            Reconstructed features of shape (batch_size, input_dim)

        Example:
            >>> z = torch.randn(32, 16)
            >>> x_recon = model.decode(z)  # (32, 900)
        """
        return self.decoder(z)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Forward pass: encode and decode features.

        Args:
            x: Input tensor of shape (batch_size, input_dim)

        Returns:
            Reconstructed features of shape (batch_size, input_dim)

        Example:
            >>> model = Autoencoder(900, 16)
            >>> x = torch.randn(32, 900)
            >>> x_recon = model(x)  # (32, 900)
        """
        z = self.encode(x)
        return self.decode(z)

    def get_config(self) -> dict:
        """Get architecture configuration as dictionary.

        Returns:
            Configuration dict with dimensions

        Example:
            >>> model = Autoencoder(900, 16)
            >>> config = model.get_config()
            >>> print(config["input_dim"])  # 900
        """
        return {
            "input_dim": self.input_dim,
            "latent_dim": self.latent_dim,
            "hidden_dim_1": self.hidden_dim_1,
            "hidden_dim_2": self.hidden_dim_2,
        }

    def count_parameters(self) -> int:
        """Count total trainable parameters.

        Returns:
            Total number of parameters

        Example:
            >>> model = Autoencoder(900, 16)
            >>> params = model.count_parameters()
            >>> print(f"{params:,} parameters")
        """
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"Autoencoder("
            f"{self.input_dim}→{self.hidden_dim_1}→{self.hidden_dim_2}→"
            f"{self.latent_dim}→{self.hidden_dim_2}→{self.hidden_dim_1}→{self.input_dim}, "
            f"params={self.count_parameters():,}"
            f")"
        )


class VariationalAutoencoder(nn.Module):
    """Variational autoencoder with KL divergence regularization.

    Extends standard autoencoder with probabilistic interpretation: learns to
    encode features as mean and variance of latent distribution, enforcing
    smooth latent space via KL divergence regularization.

    This enables better generalization and meaningful interpolation in latent
    space, but at cost of slight reconstruction fidelity.

    Architecture:
        - Encoder: input_dim → hidden → latent_dim×2 (mean + log_var)
        - Decoder: latent_dim → hidden → input_dim
        - Loss: MSE + β×KL(N(μ,σ) || N(0,I))

    Attributes:
        encoder: Maps features to latent distribution parameters
        decoder: Maps latent samples back to feature space
        beta: Weight of KL divergence term (default 0.1)
    """

    def __init__(
        self,
        input_dim: int,
        latent_dim: int,
        hidden_dim: Optional[int] = None,
        beta: float = 0.1,
    ):
        """Initialize variational autoencoder.

        Args:
            input_dim: Input feature dimension
            latent_dim: Latent embedding dimension
            hidden_dim: Hidden layer size (default auto-calculated)
            beta: Weight of KL divergence term in loss (0 = standard AE, 1 = full VAE)

        Note:
            VAE is more computationally expensive than standard AE but provides
            better generalization properties.
        """
        super().__init__()

        if input_dim <= 0:
            raise ValueError(f"input_dim must be > 0, got {input_dim}")
        if latent_dim <= 0:
            raise ValueError(f"latent_dim must be > 0, got {latent_dim}")
        if not 0 <= beta <= 1:
            raise ValueError(f"beta must be in [0,1], got {beta}")

        self.input_dim = input_dim
        self.latent_dim = latent_dim
        self.beta = beta

        if hidden_dim is None:
            hidden_dim = min(256, max(64, input_dim // 2))

        self.hidden_dim = hidden_dim

        # Encoder: project to latent mean and log variance
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
        )

        # Latent distribution parameters
        self.fc_mean = nn.Linear(hidden_dim, latent_dim)
        self.fc_log_var = nn.Linear(hidden_dim, latent_dim)

        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, input_dim),
            nn.Sigmoid(),
        )

    def encode(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """Encode to latent distribution parameters.

        Args:
            x: Input tensor (batch_size, input_dim)

        Returns:
            Tuple of (mean, log_variance)
        """
        h = self.encoder(x)
        mean = self.fc_mean(h)
        log_var = self.fc_log_var(h)
        return mean, log_var

    def reparameterize(self, mean: torch.Tensor, log_var: torch.Tensor) -> torch.Tensor:
        """Sample from latent distribution using reparameterization trick.

        Args:
            mean: Mean of latent distribution
            log_var: Log variance of latent distribution

        Returns:
            Latent samples z = μ + σ × ε where ε ~ N(0,I)
        """
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        z = mean + eps * std
        return z

    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass with KL regularization.

        Args:
            x: Input tensor (batch_size, input_dim)

        Returns:
            Tuple of (reconstructed, mean, log_variance)
        """
        mean, log_var = self.encode(x)
        z = self.reparameterize(mean, log_var)
        recon = self.decoder(z)
        return recon, mean, log_var

    def get_config(self) -> dict:
        """Get configuration."""
        return {
            "input_dim": self.input_dim,
            "latent_dim": self.latent_dim,
            "hidden_dim": self.hidden_dim,
            "beta": self.beta,
        }

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"VariationalAutoencoder("
            f"{self.input_dim}→{self.hidden_dim}→{self.latent_dim}, "
            f"beta={self.beta}"
            f")"
        )
