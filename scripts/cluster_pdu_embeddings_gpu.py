#!/usr/bin/env python3
"""
GPU-accelerated PDU clustering with RAPIDS CUML fallback to CPU HDBSCAN.

Requires RAPIDS CUML for GPU acceleration:
  pip install cuml-cu12  # For CUDA 12.x

Automatically falls back to CPU HDBSCAN if GPU not available.
"""

import argparse
import csv
import logging
from pathlib import Path

import numpy as np
from sklearn.preprocessing import StandardScaler

from kmers.logging_utils import add_logging_args, configure_logging

logger = logging.getLogger("cluster_pdu_embeddings_gpu")

# Try to import RAPIDS, fall back to CPU
try:
    from cuml import HDBSCAN as cumlHDBSCAN
    CUML_AVAILABLE = True
    logger.info("RAPIDS CUML HDBSCAN available - will use GPU acceleration")
except ImportError:
    CUML_AVAILABLE = False
    import hdbscan
    logger.warning("RAPIDS CUML not available - falling back to CPU HDBSCAN. For 10-50x speedup, install: pip install cuml-cu12")


def main():
    parser = argparse.ArgumentParser(
        description="Cluster PDU embeddings per reference amino-acid class (GPU-accelerated or CPU fallback)."
    )
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
    parser.add_argument(
        "--min-samples",
        type=int,
        default=5,
        help="HDBSCAN min_samples. Lower values reduce noise; higher values are more conservative.",
    )
    parser.add_argument("--standardize", action="store_true", help="Standardize embedding dimensions before clustering.")
    parser.add_argument(
        "--force-cpu",
        action="store_true",
        help="Force CPU HDBSCAN even if RAPIDS is available (for testing).",
    )
    add_logging_args(parser)
    args = parser.parse_args()
    configure_logging(args.log_file, args.log_level)

    embeddings_dir = Path(args.embeddings_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    paths = sorted(embeddings_dir.glob("pdu_embedding_*.npz"))
    if args.aa:
        paths = [embeddings_dir / f"pdu_embedding_{args.aa}.npz"]

    if not paths:
        raise SystemExit(f"No embedding files found in {embeddings_dir}")

    use_gpu = CUML_AVAILABLE and not args.force_cpu
    logger.info("Using %s for HDBSCAN clustering", "GPU (RAPIDS CUML)" if use_gpu else "CPU")

    summary_rows = []
    for path in paths:
        if not path.exists():
            logger.warning("Skipping missing embedding file: %s", path)
            continue

        aa = path.stem.replace("pdu_embedding_", "")
        data = np.load(path)
        pdu_ids = data["pdu_ids"].astype(np.int64)
        X = data["Z"] if args.space == "latent" else data["coords"]

        if args.standardize:
            X = StandardScaler().fit_transform(X)
            X_device = X
        else:
            # GPU code expects float32
            X = X.astype(np.float32)
            X_device = X

        logger.info("Clustering %s: %d PDUs in %dD space", aa, len(pdu_ids), X.shape[1])

        # Cluster with appropriate backend
        if use_gpu:
            try:
                clusterer = cumlHDBSCAN(
                    min_cluster_size=args.min_cluster_size,
                    min_samples=args.min_samples,
                    verbose=logger.level <= logging.DEBUG,
                )
                labels = clusterer.fit_predict(X_device)
                probabilities = clusterer.probabilities_ if hasattr(clusterer, "probabilities_") else np.ones(len(labels))
                logger.info("GPU clustering completed for %s", aa)
            except Exception as e:
                logger.warning("GPU clustering failed for %s: %s. Falling back to CPU.", aa, e)
                use_gpu = False
                clusterer = hdbscan.HDBSCAN(
                    min_cluster_size=args.min_cluster_size,
                    min_samples=args.min_samples,
                    prediction_data=False,
                )
                labels = clusterer.fit_predict(X)
                probabilities = clusterer.probabilities_
        else:
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
                "min_cluster_size": args.min_cluster_size,
                "min_samples": args.min_samples,
                "standardized": args.standardize,
                "backend": "GPU (RAPIDS CUML)" if use_gpu else "CPU (hdbscan)",
                "cluster_file": str(cluster_path),
            }
        )
        logger.info(
            "%s: %d clusters, %d/%d noise points, backend=%s",
            aa,
            len(cluster_labels),
            n_noise,
            len(pdu_ids),
            "GPU" if use_gpu else "CPU",
        )

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
            "min_cluster_size",
            "min_samples",
            "standardized",
            "backend",
            "cluster_file",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
