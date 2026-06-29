#!/usr/bin/env python3
"""Visualize PDU cluster data - sizes, distributions, quality metrics."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

AAS = list('ACDEFGHIKLMNPQRSTVWY')
CLUSTERS_DIR = Path('analysis/clusters_umap')
OUTPUT_DIR = Path('analysis/cluster_visualizations')
OUTPUT_DIR.mkdir(exist_ok=True)


def load_cluster_stats():
    """Load cluster statistics for all AAs."""
    stats = []

    for aa in AAS:
        cluster_file = CLUSTERS_DIR / f'pdu_clusters_{aa}_umap.csv'

        if not cluster_file.exists():
            logger.warning(f"Cluster file not found for {aa}: {cluster_file}")
            continue

        df = pd.read_csv(cluster_file)

        # Remove noise points (cluster == -1)
        assigned = df[df['cluster'] != -1]
        noise = df[df['cluster'] == -1]

        stats.append({
            'aa': aa,
            'total_pdus': len(df),
            'n_clusters': df['cluster'].max() + 1,
            'assigned_pdus': len(assigned),
            'noise_pdus': len(noise),
            'noise_pct': 100 * len(noise) / len(df),
            'avg_cluster_size': len(assigned) / (df['cluster'].max() + 1) if df['cluster'].max() >= 0 else 0,
            'min_cluster_size': assigned['cluster'].value_counts().min() if len(assigned) > 0 else 0,
            'max_cluster_size': assigned['cluster'].value_counts().max() if len(assigned) > 0 else 0,
            'avg_confidence': assigned['cluster_probability'].mean() if len(assigned) > 0 else 0,
            'min_confidence': assigned['cluster_probability'].min() if len(assigned) > 0 else 0,
        })

    return pd.DataFrame(stats)


def plot_cluster_overview():
    """Overview of clusters per AA."""
    stats = load_cluster_stats()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Number of clusters
    axes[0, 0].bar(stats['aa'], stats['n_clusters'], color='steelblue', alpha=0.7)
    axes[0, 0].set_title('Number of Clusters per AA', fontsize=12, fontweight='bold')
    axes[0, 0].set_ylabel('Cluster Count')
    axes[0, 0].set_xlabel('Amino Acid')
    axes[0, 0].grid(True, alpha=0.3)

    # Noise percentage
    axes[0, 1].bar(stats['aa'], stats['noise_pct'], color='coral', alpha=0.7)
    axes[0, 1].axhline(y=stats['noise_pct'].mean(), color='red', linestyle='--', label='Mean')
    axes[0, 1].set_title('Noise Fraction (cluster=-1)', fontsize=12, fontweight='bold')
    axes[0, 1].set_ylabel('Noise %')
    axes[0, 1].set_xlabel('Amino Acid')
    axes[0, 1].set_ylim([0, 100])
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    # Average cluster size
    axes[1, 0].bar(stats['aa'], stats['avg_cluster_size'], color='green', alpha=0.7)
    axes[1, 0].set_title('Average Cluster Size', fontsize=12, fontweight='bold')
    axes[1, 0].set_ylabel('PDUs per Cluster')
    axes[1, 0].set_xlabel('Amino Acid')
    axes[1, 0].set_yscale('log')
    axes[1, 0].grid(True, alpha=0.3)

    # Average confidence
    axes[1, 1].bar(stats['aa'], stats['avg_confidence'], color='purple', alpha=0.7)
    axes[1, 1].axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='50% threshold')
    axes[1, 1].set_title('Average Cluster Confidence', fontsize=12, fontweight='bold')
    axes[1, 1].set_ylabel('Mean Cluster Probability')
    axes[1, 1].set_xlabel('Amino Acid')
    axes[1, 1].set_ylim([0, 1])
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '01_cluster_overview.png', dpi=300, bbox_inches='tight')
    logger.info("✓ Saved: 01_cluster_overview.png")
    plt.close()


def plot_cluster_size_distribution():
    """Cluster size distribution for each AA."""
    selected_aas = ['L', 'A', 'G', 'V', 'P', 'R']

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()

    for idx, aa in enumerate(selected_aas):
        cluster_file = CLUSTERS_DIR / f'pdu_clusters_{aa}_umap.csv'

        if not cluster_file.exists():
            axes[idx].text(0.5, 0.5, f'{aa}: No data', ha='center', va='center')
            continue

        df = pd.read_csv(cluster_file)
        assigned = df[df['cluster'] != -1]

        if len(assigned) == 0:
            axes[idx].text(0.5, 0.5, f'{aa}: No clusters', ha='center', va='center')
            continue

        cluster_sizes = assigned['cluster'].value_counts()

        axes[idx].hist(cluster_sizes.values, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
        axes[idx].set_title(f'{aa}: {len(cluster_sizes)} clusters', fontweight='bold')
        axes[idx].set_xlabel('PDUs per Cluster')
        axes[idx].set_ylabel('Frequency')
        axes[idx].set_xscale('log')
        axes[idx].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '02_cluster_size_distribution.png', dpi=300, bbox_inches='tight')
    logger.info("✓ Saved: 02_cluster_size_distribution.png")
    plt.close()


def plot_confidence_distribution():
    """Distribution of cluster confidence scores."""
    selected_aas = ['L', 'A', 'G', 'V', 'P', 'R']

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()

    for idx, aa in enumerate(selected_aas):
        cluster_file = CLUSTERS_DIR / f'pdu_clusters_{aa}_umap.csv'

        if not cluster_file.exists():
            axes[idx].text(0.5, 0.5, f'{aa}: No data', ha='center', va='center')
            continue

        df = pd.read_csv(cluster_file)
        assigned = df[df['cluster'] != -1]

        if len(assigned) == 0:
            axes[idx].text(0.5, 0.5, f'{aa}: No clusters', ha='center', va='center')
            continue

        axes[idx].hist(assigned['cluster_probability'].values, bins=50, color='coral', alpha=0.7, edgecolor='black')
        axes[idx].set_title(f'{aa}: Mean={assigned["cluster_probability"].mean():.3f}', fontweight='bold')
        axes[idx].set_xlabel('Cluster Probability')
        axes[idx].set_ylabel('Frequency')
        axes[idx].set_xlim([0, 1])
        axes[idx].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '03_confidence_distribution.png', dpi=300, bbox_inches='tight')
    logger.info("✓ Saved: 03_confidence_distribution.png")
    plt.close()


def plot_cluster_quality():
    """Cluster quality: size vs confidence."""
    selected_aas = ['L', 'A', 'G', 'V', 'P', 'R']

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()

    for idx, aa in enumerate(selected_aas):
        cluster_file = CLUSTERS_DIR / f'pdu_clusters_{aa}_umap.csv'

        if not cluster_file.exists():
            axes[idx].text(0.5, 0.5, f'{aa}: No data', ha='center', va='center')
            continue

        df = pd.read_csv(cluster_file)
        assigned = df[df['cluster'] != -1]

        if len(assigned) == 0:
            axes[idx].text(0.5, 0.5, f'{aa}: No clusters', ha='center', va='center')
            continue

        # Get cluster sizes
        cluster_sizes = assigned['cluster'].value_counts()
        cluster_sizes_map = cluster_sizes.to_dict()
        assigned['size'] = assigned['cluster'].map(cluster_sizes_map)

        # Scatter: cluster size vs confidence
        scatter = axes[idx].scatter(
            assigned['size'],
            assigned['cluster_probability'],
            alpha=0.3, s=10, c=assigned['cluster_probability'], cmap='viridis'
        )
        axes[idx].set_xlabel('Cluster Size')
        axes[idx].set_ylabel('Cluster Probability')
        axes[idx].set_title(f'{aa}: Size vs Quality', fontweight='bold')
        axes[idx].set_xscale('log')
        axes[idx].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '04_cluster_quality_scatter.png', dpi=300, bbox_inches='tight')
    logger.info("✓ Saved: 04_cluster_quality_scatter.png")
    plt.close()


def generate_cluster_summary():
    """Generate summary statistics table."""
    stats = load_cluster_stats()

    summary = pd.DataFrame({
        'AA': stats['aa'],
        'Total PDUs': stats['total_pdus'].astype(int),
        'N Clusters': stats['n_clusters'].astype(int),
        'Assigned': stats['assigned_pdus'].astype(int),
        'Noise %': stats['noise_pct'].round(1),
        'Avg Size': stats['avg_cluster_size'].round(0).astype(int),
        'Size Range': stats['min_cluster_size'].astype(int).astype(str) + '-' + stats['max_cluster_size'].astype(int).astype(str),
        'Avg Conf': stats['avg_confidence'].round(3),
    })

    # Save to CSV
    summary.to_csv(OUTPUT_DIR / 'cluster_summary.csv', index=False)
    logger.info(f"✓ Saved: cluster_summary.csv")

    # Print table
    logger.info("\n" + "="*120)
    logger.info("CLUSTER SUMMARY STATISTICS")
    logger.info("="*120)
    logger.info(summary.to_string(index=False))
    logger.info("="*120)

    return summary


def main():
    logger.info("="*80)
    logger.info("Visualizing PDU Cluster Data")
    logger.info("="*80)

    if not CLUSTERS_DIR.exists():
        logger.error(f"Clusters directory not found: {CLUSTERS_DIR}")
        return

    logger.info(f"\n[1/5] Generating cluster summary...")
    summary = generate_cluster_summary()

    logger.info(f"\n[2/5] Plotting cluster overview...")
    plot_cluster_overview()

    logger.info(f"\n[3/5] Plotting cluster size distribution...")
    plot_cluster_size_distribution()

    logger.info(f"\n[4/5] Plotting confidence distribution...")
    plot_confidence_distribution()

    logger.info(f"\n[5/5] Plotting quality scatter plots...")
    plot_cluster_quality()

    logger.info("\n" + "="*80)
    logger.info("✓ Cluster visualization complete!")
    logger.info("="*80)
    logger.info(f"\nOutputs saved to: {OUTPUT_DIR}")
    logger.info(f"\nGenerated files:")
    logger.info(f"  - cluster_summary.csv")
    logger.info(f"  - 01_cluster_overview.png")
    logger.info(f"  - 02_cluster_size_distribution.png")
    logger.info(f"  - 03_confidence_distribution.png")
    logger.info(f"  - 04_cluster_quality_scatter.png")


if __name__ == "__main__":
    main()
