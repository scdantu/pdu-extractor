"""HDBSCAN clustering for PDU embeddings.

This module provides a wrapper around HDBSCAN for clustering learned PDU embeddings.
HDBSCAN is a hierarchical density-based clustering algorithm that:
- Finds natural clusters without pre-specifying count
- Assigns confidence scores (0=noise, 1=core point)
- Handles variable-density clusters well

Example:
    >>> from pdusearch.clustering import HDBSCANClusterer
    >>> import numpy as np
    >>>
    >>> Z = np.random.randn(1000, 16).astype(np.float32)  # Embeddings
    >>> clusterer = HDBSCANClusterer(min_cluster_size=200)
    >>> labels, confidences = clusterer.cluster(Z)
    >>> noise_pct = (labels == -1).mean() * 100
    >>> print(f"Clustering complete: {noise_pct:.1f}% noise")
"""

import numpy as np
import logging
from typing import Tuple, Optional

logger = logging.getLogger(__name__)


class HDBSCANClusterer:
    """HDBSCAN clustering wrapper for PDU embeddings.

    Provides a unified interface to HDBSCAN clustering with logging and
    parameter validation. Cluster labels range from 0+ for true clusters
    and -1 for noise points.

    Attributes:
        min_cluster_size: Minimum cluster size
        min_samples: Minimum samples for core point
        metric: Distance metric (euclidean, manhattan, cosine, etc.)
        cluster_selection_method: Method for selecting clusters (eom, leaf)
    """

    def __init__(
        self,
        min_cluster_size: int = 200,
        min_samples: int = 5,
        metric: str = "euclidean",
        cluster_selection_method: str = "eom",
        cluster_selection_epsilon: float = 0.0,
        allow_single_cluster: bool = False,
    ):
        """Initialize HDBSCAN clusterer.

        Args:
            min_cluster_size: Minimum cluster size for cluster membership
            min_samples: Minimum samples in neighborhood for core point
            metric: Distance metric for computing distances
            cluster_selection_method: "eom" (default) or "leaf"
            cluster_selection_epsilon: Epsilon parameter for cluster selection
            allow_single_cluster: If True, single cluster is valid result

        Raises:
            ValueError: If parameters invalid
        """
        if min_cluster_size <= 0:
            raise ValueError(f"min_cluster_size must be > 0, got {min_cluster_size}")
        if min_samples <= 0:
            raise ValueError(f"min_samples must be > 0, got {min_samples}")
        if min_samples > min_cluster_size:
            raise ValueError(
                f"min_samples ({min_samples}) cannot exceed "
                f"min_cluster_size ({min_cluster_size})"
            )

        self.min_cluster_size = min_cluster_size
        self.min_samples = min_samples
        self.metric = metric
        self.cluster_selection_method = cluster_selection_method
        self.cluster_selection_epsilon = cluster_selection_epsilon
        self.allow_single_cluster = allow_single_cluster

        logger.debug(
            f"HDBSCANClusterer initialized: "
            f"min_cluster_size={min_cluster_size}, "
            f"min_samples={min_samples}, metric={metric}"
        )

    def cluster(
        self,
        Z: np.ndarray,
        verbose: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Cluster embeddings using HDBSCAN.

        Args:
            Z: Embedding matrix of shape (n_samples, n_features)
            verbose: Whether to log clustering results

        Returns:
            Tuple of (labels, confidences) where:
                - labels: Cluster assignment per sample (-1 for noise)
                - confidences: Cluster probability per sample (0-1)

        Raises:
            ImportError: If hdbscan package not installed
            ValueError: If input data invalid

        Example:
            >>> Z = np.random.randn(1000, 16).astype(np.float32)
            >>> labels, conf = clusterer.cluster(Z)
            >>> print(f"Clusters: {labels.max() + 1}, Noise: {(labels==-1).sum()}")
        """
        # Import here to avoid hard dependency
        try:
            from hdbscan import HDBSCAN
        except ImportError:
            raise ImportError(
                "hdbscan package required. Install with: pip install hdbscan"
            )

        # Validate input
        if Z.ndim != 2:
            raise ValueError(f"Z must be 2D, got shape {Z.shape}")
        if Z.dtype != np.float32:
            Z = Z.astype(np.float32)

        n_samples = Z.shape[0]

        if verbose:
            logger.info(
                f"Clustering {n_samples:,} samples with HDBSCAN "
                f"(min_cluster_size={self.min_cluster_size})"
            )

        # Run clustering
        clusterer = HDBSCAN(
            min_cluster_size=self.min_cluster_size,
            min_samples=self.min_samples,
            metric=self.metric,
            cluster_selection_method=self.cluster_selection_method,
            cluster_selection_epsilon=self.cluster_selection_epsilon,
            allow_single_cluster=self.allow_single_cluster,
        )

        labels = clusterer.fit_predict(Z)
        confidences = clusterer.probabilities_

        # Log results
        n_clusters = labels.max() + 1
        n_noise = (labels == -1).sum()
        noise_pct = 100 * n_noise / n_samples

        if verbose:
            logger.info(
                f"✓ Clustering complete: "
                f"{n_clusters} clusters, {n_noise:,} noise ({noise_pct:.1f}%)"
            )

            # Show cluster size distribution
            unique, counts = np.unique(labels[labels != -1], return_counts=True)
            if len(counts) > 0:
                logger.debug(
                    f"  Cluster sizes: min={counts.min()}, "
                    f"max={counts.max()}, mean={counts.mean():.0f}"
                )

        return labels, confidences

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"HDBSCANClusterer("
            f"min_cluster_size={self.min_cluster_size}, "
            f"min_samples={self.min_samples}, "
            f"metric={self.metric}"
            f")"
        )
