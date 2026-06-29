#!/usr/bin/env python3
"""
Parallel stability validation: Bootstrap + Perturbation via array jobs.

Each array task runs one bootstrap or perturbation iteration independently.
Results aggregated at the end.

Usage:
  # Array job: task 0-19 (10 bootstrap + 10 perturbation)
  sbatch -N 1 --array=0-19 submit_stability_parallel.slurm
"""

import argparse
from pathlib import Path
import logging
import os

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
import hdbscan

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def cluster_embeddings(embeddings, min_cluster_size=200):
    """Cluster with multi-core support."""
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=10,
        core_dist_n_jobs=-1
    )
    return clusterer.fit_predict(embeddings)


def bootstrap_iteration(embeddings, labels, iteration, min_cluster_size=200):
    """Single bootstrap iteration."""
    logger.info(f"Bootstrap iteration {iteration}...")

    rng = np.random.RandomState(iteration)
    indices = rng.choice(len(embeddings), size=int(0.9 * len(embeddings)), replace=False)

    boot_embeddings = embeddings[indices]
    boot_labels_true = labels[indices]

    boot_labels_pred = cluster_embeddings(boot_embeddings, min_cluster_size)

    ari = adjusted_rand_score(boot_labels_true, boot_labels_pred)
    logger.info(f"  ARI: {ari:.3f}")

    return {
        'iteration': iteration,
        'type': 'bootstrap',
        'ari': ari,
    }


def perturbation_iteration(embeddings, labels, iteration, noise_std=0.01, min_cluster_size=200):
    """Single perturbation iteration."""
    logger.info(f"Perturbation iteration {iteration}...")

    rng = np.random.RandomState(iteration + 100)
    noise = rng.normal(0, noise_std, embeddings.shape)
    pert_embeddings = embeddings + noise

    pert_labels = cluster_embeddings(pert_embeddings, min_cluster_size)

    ari = adjusted_rand_score(labels, pert_labels)
    logger.info(f"  ARI: {ari:.3f}")

    return {
        'iteration': iteration,
        'type': 'perturbation',
        'ari': ari,
    }


def main():
    parser = argparse.ArgumentParser(description="Parallel stability validation")
    parser.add_argument("--task-id", type=int, required=True, help="Task ID (0-19)")
    parser.add_argument("--embeddings", required=True, help="Embeddings NPZ file")
    parser.add_argument("--aa", required=True, help="Amino acid")
    parser.add_argument("--out-dir", default="analysis/stability_parallel", help="Output directory")
    parser.add_argument("--min-cluster-size", type=int, default=200)

    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info(f"ARCHER2 Parallel Stability Validation - Task {args.task_id}")
    logger.info("=" * 80)
    logger.info(f"CPU cores: {os.cpu_count()}")

    # Load embeddings
    logger.info(f"\nLoading embeddings...")
    data = np.load(args.embeddings)
    embeddings = data['embeddings'].astype(np.float32)
    logger.info(f"Shape: {embeddings.shape}")

    # Generate or load initial labels
    logger.info(f"Generating initial cluster labels...")
    labels = cluster_embeddings(embeddings, args.min_cluster_size)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    logger.info(f"Clusters: {n_clusters}")

    # Task mapping: 0-9 = bootstrap, 10-19 = perturbation
    if args.task_id < 10:
        iteration = args.task_id
        result = bootstrap_iteration(embeddings, labels, iteration, args.min_cluster_size)
    else:
        iteration = args.task_id - 10
        result = perturbation_iteration(embeddings, labels, iteration, min_cluster_size=args.min_cluster_size)

    # Save result
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    result_file = out_dir / f"task_{args.task_id:02d}_{result['type']}.csv"
    result_df = pd.DataFrame([result])
    result_df.to_csv(result_file, index=False)

    logger.info(f"\n✓ Task {args.task_id} complete")
    logger.info(f"  Type: {result['type']}")
    logger.info(f"  Iteration: {iteration}")
    logger.info(f"  ARI: {result['ari']:.3f}")
    logger.info(f"  Result saved: {result_file}")
    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
