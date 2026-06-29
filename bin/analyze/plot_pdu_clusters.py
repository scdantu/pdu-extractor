#!/usr/bin/env python3
import argparse
import logging
import os
from pathlib import Path

import pandas as pd

from kmers.logging_utils import add_logging_args, configure_logging

os.environ.setdefault("MPLCONFIGDIR", str(Path(".cache") / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(".cache")))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

logger = logging.getLogger("plot_pdu_clusters")


def main():
    parser = argparse.ArgumentParser(description="Plot PDU clusters on 2D embedding coordinates.")
    parser.add_argument("--embeddings-dir", default="analysis/embeddings", help="Directory with pdu_embedding_*.csv files.")
    parser.add_argument("--clusters-dir", default="analysis/clusters", help="Directory with pdu_clusters_*_*.csv files.")
    parser.add_argument("--out-dir", default="analysis/plots", help="Output directory for plots.")
    parser.add_argument("--aa", required=True, help="Reference amino acid to plot, e.g. A.")
    parser.add_argument("--space", choices=("latent", "umap"), default="umap", help="Cluster label source to overlay.")
    parser.add_argument("--point-size", type=float, default=5.0, help="Scatter point size.")
    parser.add_argument("--alpha", type=float, default=0.75, help="Scatter alpha.")
    add_logging_args(parser)
    args = parser.parse_args()
    configure_logging(args.log_file, args.log_level)

    embeddings_path = Path(args.embeddings_dir) / f"pdu_embedding_{args.aa}.csv"
    clusters_path = Path(args.clusters_dir) / f"pdu_clusters_{args.aa}_{args.space}.csv"
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    embeddings = pd.read_csv(embeddings_path)
    clusters = pd.read_csv(clusters_path)
    data = embeddings.merge(clusters, on="pdu_id", how="inner")
    if data.empty:
        raise SystemExit("No overlapping pdu_id values between embedding and cluster files.")

    fig, ax = plt.subplots(figsize=(9, 7), dpi=180)
    noise = data[data["cluster"] == -1]
    clustered = data[data["cluster"] != -1]

    if not noise.empty:
        ax.scatter(
            noise["x"],
            noise["y"],
            s=args.point_size,
            c="#d0d0d0",
            alpha=0.35,
            linewidths=0,
            label="noise",
        )

    scatter = ax.scatter(
        clustered["x"],
        clustered["y"],
        s=args.point_size,
        c=clustered["cluster"],
        cmap="tab20",
        alpha=args.alpha,
        linewidths=0,
    )
    if not clustered.empty:
        colorbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        colorbar.set_label("Cluster")

    ax.set_title(f"PDU clusters for reference residue {args.aa} ({args.space} HDBSCAN labels)")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.spines[["top", "right"]].set_visible(False)
    if not noise.empty:
        ax.legend(loc="best", frameon=False)

    output_path = out_dir / f"pdu_clusters_{args.aa}_{args.space}.png"
    fig.tight_layout()
    fig.savefig(output_path)
    logger.info("Wrote %s", output_path)


if __name__ == "__main__":
    main()
