#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


class Autoencoder(nn.Module):
    def __init__(self, input_dim, latent_dim):
        super().__init__()
        hidden_1 = min(512, max(64, input_dim // 2))
        hidden_2 = min(256, max(32, input_dim // 4))
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_1),
            nn.ReLU(),
            nn.Linear(hidden_1, hidden_2),
            nn.ReLU(),
            nn.Linear(hidden_2, latent_dim),
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_2),
            nn.ReLU(),
            nn.Linear(hidden_2, hidden_1),
            nn.ReLU(),
            nn.Linear(hidden_1, input_dim),
            nn.Sigmoid(),
        )

    def forward(self, x):
        z = self.encoder(x)
        return self.decoder(z)


def main():
    parser = argparse.ArgumentParser(description="Train AA-specific PDU autoencoders and export embeddings.")
    parser.add_argument("--features-dir", default="analysis/features", help="Directory containing pdu_features_*.npz.")
    parser.add_argument("--out-dir", default="analysis/embeddings", help="Output directory for embeddings and maps.")
    parser.add_argument("--aa", default=None, help="Optional single reference amino acid to process, e.g. A.")
    parser.add_argument("--latent-dim", type=int, default=16, help="Autoencoder embedding dimension.")
    parser.add_argument("--epochs", type=int, default=40, help="Training epochs per AA class.")
    parser.add_argument("--batch-size", type=int, default=256, help="Batch size.")
    parser.add_argument("--learning-rate", type=float, default=1e-3, help="Learning rate.")
    parser.add_argument("--seed", type=int, default=1, help="Random seed.")
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    features_dir = Path(args.features_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    paths = sorted(features_dir.glob("pdu_features_*.npz"))
    if args.aa:
        paths = [features_dir / f"pdu_features_{args.aa}.npz"]

    if not paths:
        raise SystemExit(f"No feature files found in {features_dir}")

    device = select_device()
    print(f"Using PyTorch device: {device}")
    for path in paths:
        if not path.exists():
            print(f"Skipping missing feature file: {path}")
            continue
        data = np.load(path, allow_pickle=True)
        reference_aa = str(data["reference_aa"][0])
        X = data["X"].astype(np.float32)
        pdu_ids = data["pdu_ids"].astype(np.int64)

        model = Autoencoder(X.shape[1], args.latent_dim).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate)
        loss_fn = nn.MSELoss()
        loader = DataLoader(TensorDataset(torch.from_numpy(X)), batch_size=args.batch_size, shuffle=True)

        model.train()
        for epoch in range(1, args.epochs + 1):
            total_loss = 0.0
            total_rows = 0
            for (batch,) in loader:
                batch = batch.to(device)
                optimizer.zero_grad()
                reconstructed = model(batch)
                loss = loss_fn(reconstructed, batch)
                loss.backward()
                optimizer.step()
                total_loss += float(loss.item()) * len(batch)
                total_rows += len(batch)
            if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
                print(f"{reference_aa} epoch {epoch}/{args.epochs} loss={total_loss / total_rows:.6f}")

        model.eval()
        with torch.no_grad():
            Z = model.encoder(torch.from_numpy(X).to(device)).cpu().numpy()

        coords, method = reduce_to_2d(Z, args.seed)
        np.savez_compressed(out_dir / f"pdu_embedding_{reference_aa}.npz", Z=Z, coords=coords, pdu_ids=pdu_ids)
        write_csv(out_dir / f"pdu_embedding_{reference_aa}.csv", pdu_ids, Z, coords, method)
        torch.save(model.state_dict(), out_dir / f"pdu_autoencoder_{reference_aa}.pt")
        print(f"{reference_aa}: wrote embeddings with {method} coordinates")


def reduce_to_2d(Z, seed):
    try:
        import umap

        coords = umap.UMAP(n_components=2, random_state=seed).fit_transform(Z)
        return coords.astype(np.float32), "umap"
    except ModuleNotFoundError:
        centered = Z - Z.mean(axis=0, keepdims=True)
        _, _, vt = np.linalg.svd(centered, full_matrices=False)
        coords = centered @ vt[:2].T
        return coords.astype(np.float32), "pca_fallback"


def select_device():
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def write_csv(path, pdu_ids, Z, coords, method):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pdu_id", "x", "y", "reduction_method"] + [f"z_{i}" for i in range(Z.shape[1])])
        for pdu_id, coord, embedding in zip(pdu_ids, coords, Z):
            writer.writerow([int(pdu_id), float(coord[0]), float(coord[1]), method] + [float(v) for v in embedding])


if __name__ == "__main__":
    main()
