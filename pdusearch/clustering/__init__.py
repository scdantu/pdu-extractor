"""Clustering algorithms and quality metrics for PDU analysis.

This module provides HDBSCAN clustering and comprehensive quality metrics
for assessing clustering results.

Main Classes:
    - HDBSCANClusterer: Density-based clustering wrapper
    - ClusterMetrics: Cluster quality assessment

Metrics:
    - noise_percent: Fraction of noise points
    - avg_confidence: Average cluster membership probability
    - cluster_size_stats: Statistics about cluster sizes
    - compute_silhouette_score: External cluster cohesion metric

Example:
    >>> from pdusearch.clustering import HDBSCANClusterer, ClusterMetrics
    >>> import numpy as np
    >>>
    >>> Z = np.random.randn(1000, 16).astype(np.float32)
    >>> clusterer = HDBSCANClusterer(min_cluster_size=200)
    >>> labels, conf = clusterer.cluster(Z)
    >>>
    >>> metrics = ClusterMetrics(labels, conf)
    >>> print(f"Noise: {metrics.noise_percent:.1f}%")
    >>> print(f"Clusters: {metrics.n_clusters}")
"""

from .hdbscan import HDBSCANClusterer
from .metrics import ClusterMetrics, compute_silhouette_score

__all__ = [
    "HDBSCANClusterer",
    "ClusterMetrics",
    "compute_silhouette_score",
]
