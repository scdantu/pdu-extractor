#!/usr/bin/env python3
"""
Cluster stability validation: Bootstrap + Perturbation + Cross-validation.

Tests if clusters are robust or artifacts of the pipeline.

Metrics:
  - Adjusted Rand Index (ARI) >0.85 = stable
  - Perturbation: add noise to embeddings → re-cluster → stability
  - Bootstrap: resample 90% → re-cluster → ARI
  - Silhouette: cluster cohesion

Usage:
  python3 validate_clusters_stability.py \
    --embeddings analysis/embeddings/pdu_embedding_L.npz \
    --aa L \
    --n-bootstrap 10 \
    --out-dir analysis/validation
"""

import argparse
from pathlib import Path
import logging

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, silhouette_score
import hdbscan

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def cluster_embeddings(embeddings, min_cluster_size=200, min_samples=10):
    """Cluster embeddings with HDBSCAN."""
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        core_dist_n_jobs=-1
    )
    return clusterer.fit_predict(embeddings)


def bootstrap_stability(embeddings, labels, n_bootstrap=10, min_cluster_size=200):
    """Test cluster stability via bootstrap resampling."""
    logger.info(f"\nBootstrap validation ({n_bootstrap} iterations)...")

    ari_scores = []
    rng = np.random.RandomState(42)

    for i in range(n_bootstrap):
        # Resample 90%
        indices = rng.choice(len(embeddings), size=int(0.9 * len(embeddings)), replace=False)
        boot_embeddings = embeddings[indices]
        boot_labels_true = labels[indices]

        # Re-cluster
        boot_labels_pred = cluster_embeddings(boot_embeddings, min_cluster_size)

        # Compare to original assignment at same indices
        ari = adjusted_rand_score(boot_labels_true, boot_labels_pred)
        ari_scores.append(ari)

        logger.info(f"  Bootstrap {i+1}/{n_bootstrap}: ARI = {ari:.3f}")

    mean_ari = np.mean(ari_scores)
    std_ari = np.std(ari_scores)

    logger.info(f"\nBootstrap Results:")
    logger.info(f"  Mean ARI: {mean_ari:.3f} ± {std_ari:.3f}")
    logger.info(f"  Min ARI: {np.min(ari_scores):.3f}")
    logger.info(f"  Max ARI: {np.max(ari_scores):.3f}")

    if mean_ari > 0.85:
        logger.info(f"  ✓ STABLE (ARI >0.85)")
    else:
        logger.info(f"  ⚠ UNSTABLE (ARI <0.85)")

    return ari_scores


def perturbation_stability(embeddings, labels, n_perturbations=5, noise_std=0.01):
    """Test cluster stability under perturbation."""
    logger.info(f"\nPerturbation validation ({n_perturbations} iterations, noise={noise_std})...")

    ari_scores = []
    rng = np.random.RandomState(42)

    for i in range(n_perturbations):
        # Add Gaussian noise
        noise = rng.normal(0, noise_std, embeddings.shape)
        pert_embeddings = embeddings + noise

        # Re-cluster
        pert_labels = cluster_embeddings(pert_embeddings)

        # Compare to original
        ari = adjusted_rand_score(labels, pert_labels)
        ari_scores.append(ari)

        logger.info(f"  Perturbation {i+1}/{n_perturbations}: ARI = {ari:.3f}")

    mean_ari = np.mean(ari_scores)
    std_ari = np.std(ari_scores)

    logger.info(f"\nPerturbation Results:")
    logger.info(f"  Mean ARI: {mean_ari:.3f} ± {std_ari:.3f}")

    if mean_ari > 0.85:
        logger.info(f"  ✓ ROBUST (ARI >0.85)")
    else:
        logger.info(f"  ⚠ SENSITIVE (ARI <0.85)")

    return ari_scores


def silhouette_analysis(embeddings, labels):
    """Calculate silhouette score (cluster cohesion)."""
    logger.info(f"\nSilhouette analysis...")

    # Remove noise points (-1)
    valid_mask = labels != -1
    valid_embeddings = embeddings[valid_mask]
    valid_labels = labels[valid_mask]

    if len(set(valid_labels)) < 2:
        logger.warning("  Not enough clusters for silhouette score")
        return None

    score = silhouette_score(valid_embeddings, valid_labels, sample_size=10000)
    logger.info(f"  Silhouette Score: {score:.3f}")
    logger.info(f"  (Range: -1 to 1; >0.5 = good separation)")

    return score


def main():
    parser = argparse.ArgumentParser(description="Cluster stability validation")
    parser.add_argument("--embeddings", required=True, help="Embeddings NPZ file")
    parser.add_argument("--clusters-csv", default=None, help="Cluster assignments CSV (optional)")
    parser.add_argument("--aa", required=True, help="Amino acid")
    parser.add_argument("--n-bootstrap", type=int, default=10, help="Bootstrap iterations")
    parser.add_argument("--n-perturbations", type=int, default=5, help="Perturbation iterations")
    parser.add_argument("--min-cluster-size", type=int, default=200, help="HDBSCAN min_cluster_size")
    parser.add_argument("--out-dir", default="analysis/validation", help="Output directory")

    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info("Cluster Stability Validation")
    logger.info("=" * 80)

    # Load embeddings
    logger.info(f"\nLoading embeddings from {args.embeddings}...")
    data = np.load(args.embeddings)
    embeddings = data['embeddings'].astype(np.float32)
    logger.info(f"Shape: {embeddings.shape}")

    # Get or generate cluster assignments
    if args.clusters_csv:
        logger.info(f"Loading cluster assignments from {args.clusters_csv}...")
        df = pd.read_csv(args.clusters_csv)
        labels = df['cluster'].values
    else:
        logger.info(f"Generating initial clusters...")
        labels = cluster_embeddings(embeddings, args.min_cluster_size)

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)
    logger.info(f"Clusters: {n_clusters}, Noise: {n_noise} ({100*n_noise/len(labels):.1f}%)")

    # Validations
    boot_scores = bootstrap_stability(embeddings, labels, args.n_bootstrap, args.min_cluster_size)
    pert_scores = perturbation_stability(embeddings, labels, args.n_perturbations)
    silhouette = silhouette_analysis(embeddings, labels)

    # Summary
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results = {
        'aa': args.aa,
        'n_clusters': n_clusters,
        'n_pdus': len(embeddings),
        'n_noise': n_noise,
        'bootstrap_mean_ari': np.mean(boot_scores),
        'bootstrap_std_ari': np.std(boot_scores),
        'perturbation_mean_ari': np.mean(pert_scores),
        'perturbation_std_ari': np.std(pert_scores),
        'silhouette_score': silhouette if silhouette else None,
    }

    results_df = pd.DataFrame([results])
    results_file = out_dir / f"stability_results_{args.aa}.csv"
    results_df.to_csv(results_file, index=False)

    logger.info(f"\n✓ Results saved: {results_file}")
    logger.info("\n" + "=" * 80)
    logger.info("SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Bootstrap ARI: {results['bootstrap_mean_ari']:.3f} ± {results['bootstrap_std_ari']:.3f}")
    logger.info(f"Perturbation ARI: {results['perturbation_mean_ari']:.3f} ± {results['perturbation_std_ari']:.3f}")
    if silhouette:
        logger.info(f"Silhouette Score: {silhouette:.3f}")

    if results['bootstrap_mean_ari'] > 0.85 and results['perturbation_mean_ari'] > 0.85:
        logger.info(f"\n✓ CLUSTERS ARE STABLE AND ROBUST")
    else:
        logger.info(f"\n⚠ CLUSTERS MAY BE UNSTABLE")

    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
