#!/usr/bin/env python3
"""
Train PDU autoencoder on ARCHER2 using distributed training.

Uses PyTorch DistributedDataParallel (DDP) for multi-node, multi-core parallelization.
Optimized for CPU nodes (128-core AMD EPYC with 256GB RAM).

Expected speedup: 20-30x on single node, 60-100x+ on 4 nodes vs Mac.

Run via Slurm:
  sbatch submit_archer2.slurm
"""

import os
import argparse
from pathlib import Path
import logging

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from torch.nn.parallel import DistributedDataParallel as DDP
import torch.distributed as dist
from torch.optim import Adam

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

    def decode(self, z):
        return self.decoder(z)

    def forward(self, x):
        z = self.encode(x)
        return self.decode(z)


def setup_distributed():
    """Initialize distributed training."""
    if 'RANK' in os.environ and 'WORLD_SIZE' in os.environ:
        rank = int(os.environ['RANK'])
        world_size = int(os.environ['WORLD_SIZE'])
        dist.init_process_group(backend='gloo')  # CPU-friendly backend
        return rank, world_size
    return 0, 1


def main():
    parser = argparse.ArgumentParser(description="Train PDU autoencoder on ARCHER2")
    parser.add_argument("--features-dir", required=True, help="Directory with NPZ feature files")
    parser.add_argument("--out-dir", default="embeddings", help="Output directory")
    parser.add_argument("--aa", required=True, help="Amino acid (e.g., L)")
    parser.add_argument("--latent-dim", type=int, default=16, help="Latent dimension")
    parser.add_argument("--epochs", type=int, default=15, help="Training epochs")
    parser.add_argument("--batch-size", type=int, default=512, help="Batch size (can be large on ARCHER2)")
    parser.add_argument("--learning-rate", type=float, default=0.001, help="Learning rate")
    parser.add_argument("--validation-fraction", type=float, default=0.1, help="Validation fraction")

    args = parser.parse_args()

    # Setup distributed
    rank, world_size = setup_distributed()
    device = torch.device('cpu')

    # Logging
    if rank == 0:
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        logger = logging.getLogger(__name__)
        logger.info(f"Training on ARCHER2: {world_size} processes, rank {rank}")
        logger.info(f"Batch size: {args.batch_size}, Epochs: {args.epochs}")
    else:
        logging.basicConfig(level=logging.WARNING)
        logger = logging.getLogger(__name__)

    # Load features
    features_file = Path(args.features_dir) / f"pdu_features_{args.aa}.npz"
    if rank == 0:
        logger.info(f"Loading features from {features_file}...")

    data = np.load(features_file)
    features = data['features'].astype(np.float32)

    if rank == 0:
        logger.info(f"Loaded {len(features):,} PDU features (shape: {features.shape})")

    # Split train/val
    n_val = max(1, int(len(features) * args.validation_fraction))
    n_train = len(features) - n_val

    indices = np.random.permutation(len(features))
    train_indices = indices[:n_train]
    val_indices = indices[n_train:]

    train_features = features[train_indices]
    val_features = features[val_indices]

    # Create datasets
    train_dataset = TensorDataset(torch.from_numpy(train_features))
    val_dataset = TensorDataset(torch.from_numpy(val_features))

    # Distributed sampler (each rank gets different subset)
    from torch.utils.data import DistributedSampler

    train_sampler = DistributedSampler(
        train_dataset,
        num_replicas=world_size,
        rank=rank,
        shuffle=True
    )

    val_sampler = DistributedSampler(
        val_dataset,
        num_replicas=world_size,
        rank=rank,
        shuffle=False
    )

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, sampler=train_sampler, num_workers=4)
    val_loader = DataLoader(val_dataset, batch_size=args.batch_size, sampler=val_sampler, num_workers=4)

    # Model
    model = PDUAutoencoder(input_dim=900, latent_dim=args.latent_dim)
    model = model.to(device)

    # Distributed wrapper
    if world_size > 1:
        model = DDP(model, device_ids=None)  # CPU training

    # Optimizer
    optimizer = Adam(model.parameters(), lr=args.learning_rate)
    loss_fn = nn.MSELoss()

    # Training loop
    best_val_loss = float('inf')

    for epoch in range(args.epochs):
        model.train()
        train_sampler.set_epoch(epoch)

        train_loss = 0.0
        for batch_idx, (features,) in enumerate(train_loader):
            features = features.to(device)

            optimizer.zero_grad()
            reconstructed = model(features)
            loss = loss_fn(reconstructed, features)
            loss.backward()
            optimizer.step()

            train_loss += loss.item()

            if rank == 0 and (batch_idx + 1) % 100 == 0:
                logger.info(
                    f"Epoch {epoch+1}/{args.epochs}, Batch {batch_idx+1}/{len(train_loader)}, "
                    f"Loss: {loss.item():.6f}"
                )

        # Validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for features, in val_loader:
                features = features.to(device)
                reconstructed = model(features)
                loss = loss_fn(reconstructed, features)
                val_loss += loss.item()

        # Average across ranks
        if world_size > 1:
            train_loss_tensor = torch.tensor([train_loss], device=device)
            val_loss_tensor = torch.tensor([val_loss], device=device)
            dist.all_reduce(train_loss_tensor)
            dist.all_reduce(val_loss_tensor)
            train_loss = train_loss_tensor.item() / world_size
            val_loss = val_loss_tensor.item() / world_size

        avg_train_loss = train_loss / len(train_loader)
        avg_val_loss = val_loss / len(val_loader)

        if rank == 0:
            logger.info(
                f"Epoch {epoch+1}/{args.epochs}: "
                f"train_loss={avg_train_loss:.6f}, val_loss={avg_val_loss:.6f}"
            )

            if avg_val_loss < best_val_loss:
                best_val_loss = avg_val_loss
                logger.info(f"  → New best validation loss!")

    # Save embeddings (only on rank 0)
    if rank == 0:
        logger.info("Extracting embeddings...")

        model.eval()
        embeddings = []

        with torch.no_grad():
            for features, in DataLoader(TensorDataset(torch.from_numpy(features)), batch_size=args.batch_size):
                features = features.to(device)
                z = model.module.encode(features) if isinstance(model, DDP) else model.encode(features)
                embeddings.append(z.cpu().numpy())

        embeddings = np.concatenate(embeddings, axis=0)

        # Save
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        embedding_file = out_dir / f"pdu_embedding_{args.aa}.npz"
        np.savez_compressed(embedding_file, embeddings=embeddings)

        logger.info(f"Saved embeddings: {embedding_file}")
        logger.info(f"Shape: {embeddings.shape}")
        logger.info(f"Best validation loss: {best_val_loss:.6f}")

    # Cleanup
    if world_size > 1:
        dist.destroy_process_group()


if __name__ == "__main__":
    main()
