#!/usr/bin/env python3
"""Train autoencoders and cluster for compression hypothesis test.

For each (distance, latent_dim) combination:
1. Load features
2. Train autoencoder
3. Cluster embeddings
4. Record stats
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from hdbscan import HDBSCAN
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


class Autoencoder(nn.Module):
    def __init__(self, input_dim, latent_dim):
        super().__init__()
        h1 = min(512, max(64, input_dim // 2))
        h2 = min(256, max(32, input_dim // 4))
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, h1), nn.ReLU(),
            nn.Linear(h1, h2), nn.ReLU(),
            nn.Linear(h2, latent_dim)
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, h2), nn.ReLU(),
            nn.Linear(h2, h1), nn.ReLU(),
            nn.Linear(h1, input_dim), nn.Sigmoid()
        )

    def forward(self, x):
        return self.decoder(self.encoder(x))


def train_autoencoder(X, latent_dim, epochs, batch_size, device='cpu'):
    """Train autoencoder and return embeddings."""
    model = Autoencoder(X.shape[1], latent_dim).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    loss_fn = nn.MSELoss()

    X_tensor = torch.FloatTensor(X).to(device)
    loader = DataLoader(TensorDataset(X_tensor), batch_size=batch_size, shuffle=True)

    for epoch in range(epochs):
        losses = []
        for (batch,) in loader:
            opt.zero_grad()
            loss = loss_fn(model(batch), batch)
            loss.backward()
            opt.step()
            losses.append(loss.item())

        if (epoch + 1) % max(1, epochs // 5) == 0:
            logger.info(f"  Epoch {epoch+1}/{epochs}: loss={np.mean(losses):.4f}")

    # Get embeddings
    with torch.no_grad():
        Z = model.encoder(X_tensor).cpu().numpy()

    return Z


def main():
    if len(sys.argv) < 4:
        print("Usage: train_and_cluster_compression_test.py <aa> <distance> <latent_dim>")
        sys.exit(1)

    aa = sys.argv[1]
    distance = int(sys.argv[2])
    latent_dim = int(sys.argv[3])

    # Paths
    temp_dir = Path(f"analysis/compression_test/temp_{aa}_{distance}a")
    features_file = temp_dir / f"pdu_features_{aa}.npz"

    if not features_file.exists():
        logger.error(f"Features not found: {features_file}")
        sys.exit(1)

    # Load features
    logger.info(f"Loading features: {latent_dim}D from {features_file.name}")
    data = np.load(features_file, allow_pickle=True)
    X = data['X'].astype(np.float32)
    pdu_ids = data.get('pdu_ids', np.arange(len(X)))

    logger.info(f"  {X.shape[0]:,} PDUs, {X.shape[1]} features")

    # Train autoencoder
    logger.info(f"Training autoencoder ({X.shape[1]}→{latent_dim}D)...")
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    Z = train_autoencoder(X, latent_dim, epochs=10, batch_size=32, device=device)

    # Save embeddings
    emb_file = temp_dir / f"embedding_{aa}_{distance}a_{latent_dim}d.npz"
    np.savez(emb_file, Z=Z, pdu_ids=pdu_ids, coords=np.zeros((len(Z), 2)))
    logger.info(f"Saved embeddings: {emb_file.name}")

    # Cluster
    logger.info(f"Clustering...")
    clusterer = HDBSCAN(min_cluster_size=200, min_samples=5)
    labels = clusterer.fit_predict(Z)
    confidences = clusterer.probabilities_

    # Save clusters
    cluster_file = temp_dir / f"clusters_{aa}_{distance}a_{latent_dim}d.csv"
    df = pd.DataFrame({
        'pdu_id': pdu_ids,
        'cluster': labels,
        'cluster_probability': confidences
    })
    df.to_csv(cluster_file, index=False)

    # Stats
    assigned = df[df['cluster'] != -1]
    noise_pct = 100 * (df['cluster'] == -1).sum() / len(df)
    avg_conf = assigned['cluster_probability'].mean() if len(assigned) > 0 else 0

    logger.info(f"Results: {labels.max()+1} clusters, {noise_pct:.1f}% noise, {avg_conf:.3f} conf")
    print(f"{noise_pct:.1f},{avg_conf:.3f}")


if __name__ == "__main__":
    main()
