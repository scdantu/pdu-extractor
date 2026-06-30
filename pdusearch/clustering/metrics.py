"""Cluster quality metrics for PDU clustering evaluation.

This module provides metrics to assess clustering quality including noise percentage,
confidence scores, cluster cohesion, and enrichment statistics.

Metrics:
    - noise_percent: Percentage of unassigned points
    - avg_confidence: Average cluster membership probability
    - n_clusters: Number of discovered clusters
    - cluster_size_stats: Min/max/mean cluster sizes
    - silhouette_score: Cluster cohesion metric

Example:
    >>> from pdusearch.clustering import ClusterMetrics
    >>> labels = np.array([0, 0, 1, 1, -1])  # -1 = noise
    >>> confidences = np.array([0.9, 0.8, 0.95, 0.85, 0.0])
    >>> metrics = ClusterMetrics(labels, confidences)
    >>> print(f"Noise: {metrics.noise_percent:.1f}%")
    >>> print(f"Quality: {metrics.avg_confidence:.3f}")
"""

import numpy as np
from typing import Dict, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


class ClusterMetrics:
    """Compute quality metrics for clustering results.

    Provides comprehensive assessment of cluster quality including noise
    levels, confidence distributions, and cluster sizes.

    Attributes:
        labels: Cluster assignment array (-1 for noise)
        confidences: Cluster probability array (0-1)
    """

    def __init__(
        self,
        labels: np.ndarray,
        confidences: np.ndarray,
        embeddings: Optional[np.ndarray] = None,
    ):
        """Initialize cluster metrics.

        Args:
            labels: Cluster assignments from HDBSCAN (-1 for noise)
            confidences: Cluster probabilities from HDBSCAN (0-1)
            embeddings: Optional embedding matrix for silhouette score

        Raises:
            ValueError: If input shapes mismatch
        """
        if labels.shape[0] != confidences.shape[0]:
            raise ValueError(
                f"labels and confidences must have same length, "
                f"got {labels.shape[0]} vs {confidences.shape[0]}"
            )

        self.labels = np.asarray(labels, dtype=np.int32)
        self.confidences = np.asarray(confidences, dtype=np.float32)
        self.embeddings = embeddings

    @property
    def n_samples(self) -> int:
        """Total number of samples."""
        return len(self.labels)

    @property
    def n_clusters(self) -> int:
        """Number of clusters (excluding noise)."""
        return int(self.labels.max() + 1)

    @property
    def n_noise(self) -> int:
        """Number of noise points."""
        return int((self.labels == -1).sum())

    @property
    def noise_percent(self) -> float:
        """Percentage of samples classified as noise."""
        return 100 * self.n_noise / self.n_samples

    @property
    def assigned_points(self) -> np.ndarray:
        """Boolean array of assigned (non-noise) points."""
        return self.labels != -1

    @property
    def n_assigned(self) -> int:
        """Number of assigned points."""
        return int(self.assigned_points.sum())

    @property
    def avg_confidence(self) -> float:
        """Average confidence of assigned points.

        Noise points typically have confidence ≈ 0, so this metric
        reflects quality of point assignments (ignoring noise).
        """
        if self.n_assigned == 0:
            return 0.0
        return float(self.confidences[self.assigned_points].mean())

    @property
    def min_confidence(self) -> float:
        """Minimum confidence among assigned points."""
        if self.n_assigned == 0:
            return 0.0
        return float(self.confidences[self.assigned_points].min())

    @property
    def max_confidence(self) -> float:
        """Maximum confidence among assigned points."""
        if self.n_assigned == 0:
            return 0.0
        return float(self.confidences[self.assigned_points].max())

    def cluster_sizes(self) -> Dict[int, int]:
        """Get size of each cluster.

        Returns:
            Dictionary mapping cluster_id -> size
        """
        sizes = {}
        for cluster_id in range(self.n_clusters):
            size = (self.labels == cluster_id).sum()
            sizes[cluster_id] = int(size)
        return sizes

    def cluster_size_stats(self) -> Dict[str, float]:
        """Get statistics about cluster sizes.

        Returns:
            Dictionary with min_size, max_size, mean_size, median_size, std_size
        """
        sizes = list(self.cluster_sizes().values())
        if not sizes:
            return {
                "min_size": 0,
                "max_size": 0,
                "mean_size": 0.0,
                "median_size": 0.0,
                "std_size": 0.0,
            }

        return {
            "min_size": float(np.min(sizes)),
            "max_size": float(np.max(sizes)),
            "mean_size": float(np.mean(sizes)),
            "median_size": float(np.median(sizes)),
            "std_size": float(np.std(sizes)),
        }

    def confidence_stats(self) -> Dict[str, float]:
        """Get statistics about cluster confidences.

        Returns:
            Dictionary with min, max, mean, median, std of assigned confidences
        """
        if self.n_assigned == 0:
            return {
                "min_conf": 0.0,
                "max_conf": 0.0,
                "mean_conf": 0.0,
                "median_conf": 0.0,
                "std_conf": 0.0,
            }

        conf = self.confidences[self.assigned_points]
        return {
            "min_conf": float(np.min(conf)),
            "max_conf": float(np.max(conf)),
            "mean_conf": float(np.mean(conf)),
            "median_conf": float(np.median(conf)),
            "std_conf": float(np.std(conf)),
        }

    def get_summary(self) -> Dict:
        """Get comprehensive summary of clustering metrics.

        Returns:
            Dictionary with all metrics

        Example:
            >>> summary = metrics.get_summary()
            >>> for key, value in summary.items():
            ...     print(f"{key}: {value}")
        """
        cluster_stats = self.cluster_size_stats()
        conf_stats = self.confidence_stats()

        return {
            "n_samples": self.n_samples,
            "n_clusters": self.n_clusters,
            "n_noise": self.n_noise,
            "noise_percent": self.noise_percent,
            "n_assigned": self.n_assigned,
            "avg_confidence": self.avg_confidence,
            **cluster_stats,
            **conf_stats,
        }

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"ClusterMetrics("
            f"n_clusters={self.n_clusters}, "
            f"noise={self.noise_percent:.1f}%, "
            f"avg_conf={self.avg_confidence:.3f}"
            f")"
        )

    def __str__(self) -> str:
        """Formatted string representation."""
        lines = [
            f"Cluster Metrics Summary:",
            f"  Samples:        {self.n_samples:,}",
            f"  Clusters:       {self.n_clusters}",
            f"  Assigned:       {self.n_assigned:,}",
            f"  Noise:          {self.n_noise:,} ({self.noise_percent:.1f}%)",
            f"  Avg Confidence: {self.avg_confidence:.4f}",
        ]

        stats = self.cluster_size_stats()
        if stats["mean_size"] > 0:
            lines.extend([
                f"  Cluster Sizes:",
                f"    Min:    {stats['min_size']:.0f}",
                f"    Max:    {stats['max_size']:.0f}",
                f"    Mean:   {stats['mean_size']:.0f}",
                f"    Median: {stats['median_size']:.0f}",
            ])

        return "\n".join(lines)


def compute_silhouette_score(
    Z: np.ndarray,
    labels: np.ndarray,
    metric: str = "euclidean",
    sample_size: Optional[int] = None,
) -> float:
    """Compute silhouette score for clustering.

    Silhouette score measures how well-separated clusters are:
        - Score near +1: Well-separated clusters
        - Score near 0: Overlapping clusters
        - Score near -1: Misclassified points

    Can be slow for large datasets; use sample_size to subsample.

    Args:
        Z: Embeddings (n_samples, n_features)
        labels: Cluster assignments (-1 for noise)
        metric: Distance metric
        sample_size: If specified, subsample to this size (for speed)

    Returns:
        Silhouette score (-1 to 1)

    Note:
        Requires sklearn. Ignores noise points in scoring.
    """
    try:
        from sklearn.metrics import silhouette_score
    except ImportError:
        logger.warning("sklearn not available for silhouette score")
        return 0.0

    # Filter noise points
    mask = labels != -1
    if mask.sum() < 2:
        logger.warning("Not enough assigned points for silhouette score")
        return 0.0

    Z_assigned = Z[mask]
    labels_assigned = labels[mask]

    # Optionally subsample for speed
    if sample_size is not None and len(Z_assigned) > sample_size:
        idx = np.random.choice(len(Z_assigned), sample_size, replace=False)
        Z_assigned = Z_assigned[idx]
        labels_assigned = labels_assigned[idx]

    score = silhouette_score(Z_assigned, labels_assigned, metric=metric)
    return float(score)
