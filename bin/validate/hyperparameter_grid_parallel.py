#!/usr/bin/env python3
"""
Parallel hyperparameter grid search for ARCHER2 array jobs.

Each array task processes one (latent_dim, min_cluster_size) combination.
Runs 4+ tasks in parallel across multiple nodes.

Usage:
  # Single task (for testing):
  python3 hyperparameter_grid_parallel.py --task-id 0 --aa L --features analysis/features/pdu_features_L.npz

  # Via Slurm array:
  sbatch -N 2 --array=0-19 submit_grid_search_parallel.slurm
"""

import argparse
from pathlib import Path
import logging
import os

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from torch.optim import Adam
from sklearn.metrics import silhouette_score
import hdbscan

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class PDUAutoencoder(nn.Module):
    """Autoencoder with configurable latent dimension."""

    def __init__(self, input_dim=900, latent_dim=16):
        super().__init__()

        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 512),
            nn.ReLU(),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Linear(128, latent_dim)
        )

        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 128),
            nn.ReLU(),
            nn.Linear(128, 256),
            nn.ReLU(),
            nn.Linear(256, 512),
            nn.ReLU(),
            nn.Linear(512, input_dim)
        )

    def encode(self, x):
        return self.encoder(x)

    def forward(self, x):
        z = self.encode(x)
        return self.decoder(z)


def get_hyperparameters(task_id):
    """Map task ID to (latent_dim, min_cluster_size) combination."""

    latent_dims = [8, 16, 32, 64]
    min_cluster_sizes = [50, 100, 200, 300, 500]

    # Flatten grid: task_id → (latent_dim, min_cluster_size)
    num_latent = len(latent_dims)
    num_cluster = len(min_cluster_sizes)

    if task_id >= num_latent * num_cluster:
        raise ValueError(f"Task ID {task_id} out of range (max {num_latent * num_cluster - 1})")

    latent_idx = task_id // num_cluster
    cluster_idx = task_id % num_cluster

    latent_dim = latent_dims[latent_idx]
    min_cluster_size = min_cluster_sizes[cluster_idx]

    return latent_dim, min_cluster_size


def train_autoencoder(features, latent_dim, epochs=5, batch_size=1024):
    """Quick autoencoder training."""
    device = torch.device('cpu')
    dataset = TensorDataset(torch.from_numpy(features))
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    model = PDUAutoencoder(input_dim=900, latent_dim=latent_dim).to(device)
    optimizer = Adam(model.parameters(), lr=0.001)
    loss_fn = nn.MSELoss()

    for epoch in range(epochs):
        for batch, in loader:
            batch = batch.to(device)
            optimizer.zero_grad()
            recon = model(batch)
            loss = loss_fn(recon, batch)
            loss.backward()
            optimizer.step()

    # Extract embeddings
    model.eval()
    embeddings = []
    with torch.no_grad():
        for batch, in loader:
            batch = batch.to(device)
            z = model.encode(batch)
            embeddings.append(z.cpu().numpy())

    return np.concatenate(embeddings, axis=0)


def evaluate_clustering(embeddings, min_cluster_size, min_samples=10):
    """Evaluate clustering quality with multi-core support."""
    logger.info(f"HDBSCAN clustering with {os.cpu_count()} cores available...")

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        core_dist_n_jobs=-1  # Use all available cores
    )
    labels = clusterer.fit_predict(embeddings)

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)
    noise_pct = 100 * n_noise / len(labels)

    # Silhouette (on non-noise points)
    valid_mask = labels != -1
    if len(set(labels[valid_mask])) > 1:
        silhouette = silhouette_score(embeddings[valid_mask], labels[valid_mask], sample_size=5000)
    else:
        silhouette = -1

    return {
        'n_clusters': n_clusters,
        'n_noise': n_noise,
        'noise_pct': noise_pct,
        'silhouette': silhouette,
    }


def main():
    parser = argparse.ArgumentParser(description="Parallel hyperparameter grid search")
    parser.add_argument("--task-id", type=int, required=True, help="Task ID (0-19 for 4×5 grid)")
    parser.add_argument("--features", required=True, help="Features NPZ file")
    parser.add_argument("--aa", required=True, help="Amino acid")
    parser.add_argument("--out-dir", default="analysis/grid_search_parallel", help="Output directory")
    parser.add_argument("--epochs", type=int, default=5, help="Training epochs")

    args = parser.parse_args()

    # Get hyperparameters for this task
    latent_dim, min_cluster_size = get_hyperparameters(args.task_id)

    logger.info("=" * 80)
    logger.info(f"ARCHER2 Parallel Grid Search - Task {args.task_id}")
    logger.info("=" * 80)
    logger.info(f"AA: {args.aa}")
    logger.info(f"latent_dim: {latent_dim}")
    logger.info(f"min_cluster_size: {min_cluster_size}")
    logger.info(f"CPU cores: {os.cpu_count()}")

    # Load features
    logger.info(f"\nLoading features...")
    data = np.load(args.features)
    features = data['X'].astype(np.float32)
    logger.info(f"Shape: {features.shape}")

    # Train autoencoder
    logger.info(f"\n[1/2] Training autoencoder (latent_dim={latent_dim})...")
    embeddings = train_autoencoder(features, latent_dim, epochs=args.epochs)
    logger.info(f"Embeddings shape: {embeddings.shape}")

    # Evaluate clustering
    logger.info(f"\n[2/2] HDBSCAN clustering (min_cluster_size={min_cluster_size})...")
    eval_results = evaluate_clustering(embeddings, min_cluster_size)

    result = {
        'task_id': args.task_id,
        'aa': args.aa,
        'latent_dim': latent_dim,
        'min_cluster_size': min_cluster_size,
        'n_clusters': eval_results['n_clusters'],
        'n_noise': eval_results['n_noise'],
        'noise_pct': eval_results['noise_pct'],
        'silhouette': eval_results['silhouette'],
    }

    # Save result
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    result_file = out_dir / f"task_{args.task_id:02d}_L{latent_dim}_C{min_cluster_size}.csv"
    result_df = pd.DataFrame([result])
    result_df.to_csv(result_file, index=False)

    logger.info(f"\n✓ Task {args.task_id} complete")
    logger.info(f"  Clusters: {eval_results['n_clusters']}")
    logger.info(f"  Noise: {eval_results['noise_pct']:.1f}%")
    logger.info(f"  Silhouette: {eval_results['silhouette']:.3f}")
    logger.info(f"  Result saved: {result_file}")
    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
