#!/usr/bin/env python3
"""
Hyperparameter grid search for autoencoder + clustering.

Tests combinations of:
  - latent_dim: [8, 16, 32, 64]
  - min_cluster_size: [50, 100, 200, 300, 500]

Evaluates on: silhouette score, cluster count, noise percentage.

Usage:
  python3 hyperparameter_grid_search.py \
    --features analysis/features/pdu_features_L.npz \
    --aa L \
    --out-dir analysis/grid_search
"""

import argparse
from pathlib import Path
import logging
from itertools import product

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
    """Evaluate clustering quality."""
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        core_dist_n_jobs=-1
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
        'labels': labels,
    }


def main():
    parser = argparse.ArgumentParser(description="Hyperparameter grid search")
    parser.add_argument("--features", required=True, help="Features NPZ file")
    parser.add_argument("--aa", required=True, help="Amino acid")
    parser.add_argument("--latent-dims", nargs='+', type=int, default=[8, 16, 32, 64])
    parser.add_argument("--min-cluster-sizes", nargs='+', type=int, default=[50, 100, 200, 300, 500])
    parser.add_argument("--out-dir", default="analysis/grid_search", help="Output directory")
    parser.add_argument("--epochs", type=int, default=5, help="Training epochs (quick search)")

    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info("Hyperparameter Grid Search")
    logger.info("=" * 80)

    # Load features
    logger.info(f"Loading features...")
    data = np.load(args.features)
    features = data['X'].astype(np.float32)
    logger.info(f"Shape: {features.shape}")

    # Grid search
    results = []
    total = len(args.latent_dims) * len(args.min_cluster_sizes)
    iteration = 0

    for latent_dim, min_cluster_size in product(args.latent_dims, args.min_cluster_sizes):
        iteration += 1
        logger.info(f"\n[{iteration}/{total}] latent_dim={latent_dim}, min_cluster_size={min_cluster_size}")

        # Train autoencoder
        embeddings = train_autoencoder(features, latent_dim, epochs=args.epochs)
        logger.info(f"  Embeddings shape: {embeddings.shape}")

        # Evaluate clustering
        eval_results = evaluate_clustering(embeddings, min_cluster_size)

        result = {
            'aa': args.aa,
            'latent_dim': latent_dim,
            'min_cluster_size': min_cluster_size,
            'n_clusters': eval_results['n_clusters'],
            'n_noise': eval_results['n_noise'],
            'noise_pct': eval_results['noise_pct'],
            'silhouette': eval_results['silhouette'],
        }
        results.append(result)

        logger.info(f"  → Clusters: {eval_results['n_clusters']}, Noise: {eval_results['noise_pct']:.1f}%, "
                   f"Silhouette: {eval_results['silhouette']:.3f}")

    # Save results
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results_df = pd.DataFrame(results)
    results_file = out_dir / f"grid_search_results_{args.aa}.csv"
    results_df.to_csv(results_file, index=False)

    logger.info(f"\n✓ Results saved: {results_file}")

    # Find optimal hyperparameters (maximize silhouette, reasonable cluster count)
    logger.info("\n" + "=" * 80)
    logger.info("TOP CONFIGURATIONS (by silhouette score)")
    logger.info("=" * 80)

    top_results = results_df.nlargest(5, 'silhouette')
    for idx, row in top_results.iterrows():
        logger.info(f"latent_dim={row['latent_dim']}, min_cluster_size={row['min_cluster_size']}: "
                   f"Silhouette={row['silhouette']:.3f}, Clusters={row['n_clusters']}, "
                   f"Noise={row['noise_pct']:.1f}%")

    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
