#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import hdbscan
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Cluster PDU embeddings per reference amino-acid class.")
    parser.add_argument("--embeddings-dir", default="analysis/embeddings", help="Directory with pdu_embedding_*.npz files.")
    parser.add_argument("--out-dir", default="analysis/clusters", help="Output directory for cluster CSV files.")
    parser.add_argument("--aa", default=None, help="Optional single reference amino acid to cluster, e.g. A.")
    parser.add_argument(
        "--space",
        choices=("latent", "umap"),
        default="latent",
        help="Cluster autoencoder latent vectors or 2D UMAP coordinates. Latent is recommended.",
    )
    parser.add_argument("--min-cluster-size", type=int, default=50, help="HDBSCAN min_cluster_size.")
    parser.add_argument("--min-samples", type=int, default=None, help="HDBSCAN min_samples. Defaults to min_cluster_size.")
    args = parser.parse_args()

    embeddings_dir = Path(args.embeddings_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    paths = sorted(embeddings_dir.glob("pdu_embedding_*.npz"))
    if args.aa:
        paths = [embeddings_dir / f"pdu_embedding_{args.aa}.npz"]

    if not paths:
        raise SystemExit(f"No embedding files found in {embeddings_dir}")

    summary_rows = []
    for path in paths:
        if not path.exists():
            print(f"Skipping missing embedding file: {path}")
            continue

        aa = path.stem.replace("pdu_embedding_", "")
        data = np.load(path)
        pdu_ids = data["pdu_ids"].astype(np.int64)
        X = data["Z"] if args.space == "latent" else data["coords"]

        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=args.min_cluster_size,
            min_samples=args.min_samples,
            prediction_data=False,
        )
        labels = clusterer.fit_predict(X)
        probabilities = clusterer.probabilities_

        cluster_path = out_dir / f"pdu_clusters_{aa}_{args.space}.csv"
        write_cluster_csv(cluster_path, pdu_ids, labels, probabilities)

        n_noise = int(np.sum(labels == -1))
        cluster_labels = sorted(label for label in set(labels.tolist()) if label != -1)
        counts = {label: int(np.sum(labels == label)) for label in cluster_labels}
        summary_rows.append(
            {
                "aa": aa,
                "space": args.space,
                "n_pdus": len(pdu_ids),
                "n_clusters": len(cluster_labels),
                "n_noise": n_noise,
                "noise_fraction": n_noise / len(pdu_ids),
                "largest_cluster": max(counts.values()) if counts else 0,
                "cluster_file": str(cluster_path),
            }
        )
        print(f"{aa}: {len(cluster_labels)} clusters, {n_noise}/{len(pdu_ids)} noise points")

    write_summary_csv(out_dir / f"cluster_summary_{args.space}.csv", summary_rows)


def write_cluster_csv(path, pdu_ids, labels, probabilities):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pdu_id", "cluster", "cluster_probability"])
        for pdu_id, label, probability in zip(pdu_ids, labels, probabilities):
            writer.writerow([int(pdu_id), int(label), float(probability)])


def write_summary_csv(path, rows):
    with open(path, "w", newline="") as handle:
        fieldnames = [
            "aa",
            "space",
            "n_pdus",
            "n_clusters",
            "n_noise",
            "noise_fraction",
            "largest_cluster",
            "cluster_file",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
