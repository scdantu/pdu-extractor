#!/usr/bin/env python3
"""Visualize cluster separation in UMAP 2D projection."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

try:
    from umap import UMAP
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False
    logger.warning("UMAP not installed. Using random projection fallback.")

plt.style.use('seaborn-v0_8-darkgrid')

AAS = list('ACDEFGHIKLMNPQRSTVWY')
EMBEDDINGS_DIR = Path('analysis/embeddings')
CLUSTERS_DIR = Path('analysis/clusters_umap')
OUTPUT_DIR = Path('analysis/cluster_projections')
OUTPUT_DIR.mkdir(exist_ok=True)


def compute_umap_projection(embeddings, n_neighbors=15, min_dist=0.1):
    """Compute 2D UMAP projection from 16D embeddings."""
    if not HAS_UMAP:
        # Fallback: PCA-like projection
        logger.info("  Using random projection (UMAP not available)")
        u, s, vt = np.linalg.svd(embeddings, full_matrices=False)
        return u[:, :2]

    logger.info(f"  Computing UMAP (n_neighbors={n_neighbors}, min_dist={min_dist})...")
    umap_reducer = UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric='euclidean',
        random_state=42,
        n_jobs=-1
    )
    return umap_reducer.fit_transform(embeddings)


def plot_cluster_separation(aa):
    """Plot cluster separation for one amino acid."""
    embedding_file = EMBEDDINGS_DIR / f'pdu_embedding_{aa}.npz'
    cluster_file = CLUSTERS_DIR / f'pdu_clusters_{aa}_umap.csv'

    if not embedding_file.exists() or not cluster_file.exists():
        logger.warning(f"Missing files for {aa}")
        return

    logger.info(f"Processing {aa}...")

    # Load data
    embedding_data = np.load(embedding_file)

    # Check available keys
    if 'coords' in embedding_data:
        umap_2d = embedding_data['coords']  # Pre-computed UMAP coordinates
        logger.info(f"  Loaded pre-computed UMAP: {umap_2d.shape}")
    elif 'Z' in embedding_data:
        embeddings = embedding_data['Z']  # 16-dimensional latent space
        logger.info(f"  Loaded embeddings: {embeddings.shape}")
        umap_2d = compute_umap_projection(embeddings)
    else:
        logger.error(f"  Unknown keys in NPZ: {list(embedding_data.keys())}")
        return

    clusters_df = pd.read_csv(cluster_file)

    # Get PDU IDs from embedding file
    pdu_ids = embedding_data['pdu_ids']

    # Match PDU IDs to ensure alignment
    pdu_to_idx = {pdu_id: idx for idx, pdu_id in enumerate(pdu_ids)}

    # Map cluster data to embedding indices
    coords_indices = np.array([pdu_to_idx.get(pdu_id, -1) for pdu_id in clusters_df['pdu_id']])

    # Filter valid matches
    valid_mask = coords_indices >= 0
    valid_indices = coords_indices[valid_mask]
    umap_2d_matched = umap_2d[valid_indices]
    clusters_matched = clusters_df[valid_mask].copy()

    if len(valid_indices) == 0:
        logger.warning(f"  No PDU ID matches found for {aa}")
        return

    logger.info(f"  Matched {len(valid_indices):,} / {len(clusters_df):,} PDUs")

    # Separate assigned vs noise
    assigned_mask = clusters_matched['cluster'] != -1
    assigned_idx = np.where(assigned_mask)[0]
    noise_idx = np.where(~assigned_mask)[0]

    # Create figure with subplots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Plot 1: All clusters colored
    ax = axes[0]
    clusters = clusters_matched['cluster'].values
    n_clusters = int(clusters_matched['cluster'].max()) + 1

    # Plot noise points first (light gray, low alpha)
    if len(noise_idx) > 0:
        ax.scatter(
            umap_2d_matched[noise_idx, 0],
            umap_2d_matched[noise_idx, 1],
            c='lightgray', s=1, alpha=0.1, label=f'Noise ({len(noise_idx):,})'
        )

    # Plot assigned clusters with colors
    if len(assigned_idx) > 0:
        scatter = ax.scatter(
            umap_2d_matched[assigned_idx, 0],
            umap_2d_matched[assigned_idx, 1],
            c=clusters[assigned_idx],
            cmap='tab20b' if n_clusters < 20 else 'gist_ncar',
            s=5, alpha=0.6, edgecolors='none'
        )
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Cluster ID', fontsize=10)

    ax.set_title(f'{aa}: {n_clusters} clusters ({len(assigned_idx):,} assigned, {len(noise_idx):,} noise)',
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.grid(True, alpha=0.2)

    # Plot 2: Colored by confidence
    ax = axes[1]
    confidence = clusters_matched['cluster_probability'].values

    # Only plot assigned points
    if len(assigned_idx) > 0:
        scatter = ax.scatter(
            umap_2d_matched[assigned_idx, 0],
            umap_2d_matched[assigned_idx, 1],
            c=confidence[assigned_idx],
            cmap='RdYlGn', s=5, alpha=0.7, vmin=0, vmax=1, edgecolors='none'
        )
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Cluster Confidence', fontsize=10)

    ax.set_title(f'{aa}: Cluster Assignment Confidence', fontsize=12, fontweight='bold')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    output_file = OUTPUT_DIR / f'umap_clusters_{aa}.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    logger.info(f"  ✓ Saved: {output_file.name}")
    plt.close()

    return {
        'aa': aa,
        'n_clusters': n_clusters,
        'n_assigned': len(assigned_idx),
        'n_noise': len(noise_idx),
        'avg_confidence': confidence[assigned_idx].mean() if len(assigned_idx) > 0 else 0
    }


def main():
    logger.info("="*80)
    logger.info("Visualizing Cluster Separation in UMAP Space")
    logger.info("="*80)

    if not EMBEDDINGS_DIR.exists() or not CLUSTERS_DIR.exists():
        logger.error("Missing embeddings or clusters directories")
        return

    logger.info(f"\nGenerating projections for all amino acids...\n")

    stats = []
    for aa in AAS:
        try:
            result = plot_cluster_separation(aa)
            if result:
                stats.append(result)
        except Exception as e:
            logger.error(f"Error processing {aa}: {e}")

    # Summary
    logger.info("\n" + "="*80)
    logger.info("✓ Visualization complete!")
    logger.info("="*80)
    logger.info(f"\nOutputs saved to: {OUTPUT_DIR}")
    logger.info(f"\nGenerated {len(stats)} UMAP projections")

    if stats:
        summary = pd.DataFrame(stats)
        logger.info("\n" + summary.to_string(index=False))


if __name__ == "__main__":
    main()
