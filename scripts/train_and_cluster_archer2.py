#!/usr/bin/env python3
"""
Combined autoencoder training + HDBSCAN clustering on ARCHER2.

Assumes features already exist on ARCHER2.
Runs: features → autoencoder → embeddings → HDBSCAN → clusters

Much faster than separate jobs (no intermediate data transfer).

Usage:
  python3 train_and_cluster_archer2.py --aa L --features-dir analysis/features --out-dir analysis
"""

import os
import argparse
from pathlib import Path
import logging

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from torch.optim import Adam

try:
    import hdbscan
    HAS_HDBSCAN = True
except ImportError:
    HAS_HDBSCAN = False
    print("Warning: hdbscan not available, will skip clustering")

# Set up MKL for CPU optimization
os.environ['MKL_NUM_THREADS'] = '4'
os.environ['OPENBLAS_NUM_THREADS'] = '4'
os.environ['OMP_NUM_THREADS'] = '4'


class PDUAutoencoder(nn.Module):
    """Autoencoder: 900-dim → 16-dim latent → 900-dim"""

    def __init__(self, input_dim=900, latent_dim=16):
        super().__init__()

        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 512),
            nn.ReLU(),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Linear(128, latent_dim)
        )

        # Decoder
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


def train_autoencoder(features, args, logger, device):
    """Train autoencoder on features."""

    logger.info(f"Training autoencoder on {len(features):,} features...")

    # Split train/val
    n_val = max(1, int(len(features) * args.validation_fraction))
    n_train = len(features) - n_val

    indices = np.random.permutation(len(features))
    train_indices = indices[:n_train]
    val_indices = indices[n_train:]

    train_features = features[train_indices]
    val_features = features[val_indices]

    # Create datasets and loaders
    train_dataset = TensorDataset(torch.from_numpy(train_features))
    val_dataset = TensorDataset(torch.from_numpy(val_features))

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_dataset, batch_size=args.batch_size, shuffle=False, num_workers=0)

    # Model
    model = PDUAutoencoder(input_dim=900, latent_dim=args.latent_dim)
    model = model.to(device)

    # Optimizer
    optimizer = Adam(model.parameters(), lr=args.learning_rate)
    loss_fn = nn.MSELoss()

    # Training loop
    best_val_loss = float('inf')

    for epoch in range(args.epochs):
        model.train()
        train_loss = 0.0

        for batch_idx, (batch_features,) in enumerate(train_loader):
            batch_features = batch_features.to(device)

            optimizer.zero_grad()
            reconstructed = model(batch_features)
            loss = loss_fn(reconstructed, batch_features)
            loss.backward()
            optimizer.step()

            train_loss += loss.item()

            if (batch_idx + 1) % 500 == 0:
                logger.info(
                    f"Epoch {epoch+1}/{args.epochs}, Batch {batch_idx+1}/{len(train_loader)}, "
                    f"Loss: {loss.item():.6f}"
                )

        # Validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for batch_features, in val_loader:
                batch_features = batch_features.to(device)
                reconstructed = model(batch_features)
                loss = loss_fn(reconstructed, batch_features)
                val_loss += loss.item()

        avg_train_loss = train_loss / len(train_loader)
        avg_val_loss = val_loss / len(val_loader)

        logger.info(
            f"Epoch {epoch+1}/{args.epochs}: "
            f"train_loss={avg_train_loss:.6f}, val_loss={avg_val_loss:.6f}"
        )

        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss

    # Extract embeddings
    logger.info("Extracting embeddings...")
    model.eval()
    embeddings = []

    with torch.no_grad():
        for batch_features, in DataLoader(TensorDataset(torch.from_numpy(features)), batch_size=args.batch_size):
            batch_features = batch_features.to(device)
            z = model.encode(batch_features)
            embeddings.append(z.cpu().numpy())

    embeddings = np.concatenate(embeddings, axis=0)
    logger.info(f"Embeddings shape: {embeddings.shape}")

    return embeddings, best_val_loss


def cluster_embeddings(embeddings, pdu_ids, args, logger):
    """Cluster embeddings using HDBSCAN."""

    if not HAS_HDBSCAN:
        logger.warning("hdbscan not available, skipping clustering")
        return None

    logger.info(f"Clustering {len(embeddings):,} embeddings with HDBSCAN...")
    logger.info(f"min_cluster_size={args.min_cluster_size}, min_samples={args.min_samples}")

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        core_dist_n_jobs=-1  # Use all cores
    )

    labels = clusterer.fit_predict(embeddings)
    probabilities = clusterer.probabilities_

    # Statistics
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)

    logger.info(f"Found {n_clusters} clusters")
    logger.info(f"Noise points: {n_noise} ({100*n_noise/len(labels):.1f}%)")

    # Save results
    import pandas as pd
    results_df = pd.DataFrame({
        'pdu_id': pdu_ids,
        'cluster': labels,
        'cluster_probability': probabilities
    })

    return results_df


def main():
    parser = argparse.ArgumentParser(
        description="Train autoencoder + cluster on ARCHER2 (combined job)"
    )
    parser.add_argument("--aa", required=True, help="Amino acid (e.g., L)")
    parser.add_argument("--features-dir", required=True, help="Directory with NPZ feature files")
    parser.add_argument("--out-dir", required=True, help="Output directory for embeddings/clusters")

    # Autoencoder params
    parser.add_argument("--latent-dim", type=int, default=16)
    parser.add_argument("--epochs", type=int, default=15)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--learning-rate", type=float, default=0.001)
    parser.add_argument("--validation-fraction", type=float, default=0.1)

    # Clustering params
    parser.add_argument("--min-cluster-size", type=int, default=200)
    parser.add_argument("--min-samples", type=int, default=10)

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)

    device = torch.device('cpu')

    logger.info("=" * 80)
    logger.info("ARCHER2: Train Autoencoder + Cluster")
    logger.info("=" * 80)
    logger.info(f"AA: {args.aa}")
    logger.info(f"Latent dim: {args.latent_dim}, Epochs: {args.epochs}")
    logger.info(f"Batch size: {args.batch_size}")
    logger.info(f"Clustering: min_cluster_size={args.min_cluster_size}")

    # Load features
    features_file = Path(args.features_dir) / f"pdu_features_{args.aa}.npz"
    if not features_file.exists():
        logger.error(f"Features not found: {features_file}")
        exit(1)

    logger.info(f"Loading features from {features_file}...")
    data = np.load(features_file)

    # Handle multiple NPZ formats
    logger.info(f"NPZ keys available: {list(data.files)}")

    # Try different key names
    if 'features' in data.files:
        features = data['features'].astype(np.float32)
    elif 'arr_0' in data.files:
        features = data['arr_0'].astype(np.float32)
    elif len(data.files) == 1:
        # Single array in NPZ, use it
        key = data.files[0]
        features = data[key].astype(np.float32)
        logger.info(f"Using key: {key}")
    else:
        logger.error(f"Cannot find features array. Available keys: {list(data.files)}")
        exit(1)

    pdu_ids = data.get('pdu_ids', None)
    if pdu_ids is None:
        pdu_ids = np.arange(len(features))

    logger.info(f"Loaded {len(features):,} PDU features (shape: {features.shape})")

    # Step 1: Train autoencoder
    logger.info("\n" + "=" * 80)
    logger.info("STEP 1: Train Autoencoder")
    logger.info("=" * 80)

    embeddings, best_val_loss = train_autoencoder(features, args, logger, device)

    logger.info(f"Best validation loss: {best_val_loss:.6f}")

    # Step 2: Cluster
    logger.info("\n" + "=" * 80)
    logger.info("STEP 2: HDBSCAN Clustering")
    logger.info("=" * 80)

    clusters_df = cluster_embeddings(embeddings, pdu_ids, args, logger)

    # Save outputs
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    embeddings_file = out_dir / f"pdu_embedding_{args.aa}.npz"
    np.savez_compressed(embeddings_file, embeddings=embeddings)
    logger.info(f"Saved embeddings: {embeddings_file}")

    if clusters_df is not None:
        clusters_file = out_dir / f"pdu_clusters_{args.aa}_umap.csv"
        clusters_df.to_csv(clusters_file, index=False)
        logger.info(f"Saved clusters: {clusters_file}")

    logger.info("\n" + "=" * 80)
    logger.info("✓ Complete: Autoencoder trained + Clustering done")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
