#!/usr/bin/env python3
"""
Analyze all cluster CSV outputs from ARCHER2 run.

Generates comprehensive summary across all 20 amino acids.

Usage:
  python3 analyze_cluster_outputs.py --clusters-dir analysis/clusters_umap
"""

import argparse
from pathlib import Path
import logging

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def analyze_clusters(csv_file):
    """Analyze single cluster CSV."""
    df = pd.read_csv(csv_file)

    n_pdus = len(df)
    n_clusters = len(set(df['cluster'])) - (1 if -1 in df['cluster'].values else 0)
    n_noise = (df['cluster'] == -1).sum()
    noise_pct = 100 * n_noise / n_pdus

    assigned = df[df['cluster'] != -1]
    mean_prob = assigned['cluster_probability'].mean() if len(assigned) > 0 else 0

    return {
        'n_pdus': n_pdus,
        'n_clusters': n_clusters,
        'n_noise': n_noise,
        'noise_pct': noise_pct,
        'mean_prob': mean_prob,
    }


def main():
    parser = argparse.ArgumentParser(description="Analyze cluster outputs")
    parser.add_argument("--clusters-dir", required=True, help="Clusters directory")

    args = parser.parse_args()

    clusters_dir = Path(args.clusters_dir)

    logger.info("\n" + "=" * 100)
    logger.info("CLUSTER ANALYSIS ACROSS ALL 20 AMINO ACIDS")
    logger.info("=" * 100)

    # Find all cluster CSVs
    csv_files = sorted(clusters_dir.glob("pdu_clusters_*_umap.csv"))

    logger.info(f"\nFound {len(csv_files)} cluster files\n")

    results = []
    for csv_file in csv_files:
        aa = csv_file.stem.split('_')[2]

        stats = analyze_clusters(csv_file)
        stats['aa'] = aa
        results.append(stats)

        logger.info(f"AA={aa}: ✓ {stats['n_pdus']:,} PDUs, {stats['n_clusters']} clusters, {stats['noise_pct']:.1f}% noise")

    results_df = pd.DataFrame(results)

    # Sort by cluster count
    results_df = results_df.sort_values('n_clusters', ascending=False)

    logger.info("\n" + "=" * 100)
    logger.info("SUMMARY TABLE (sorted by cluster count)")
    logger.info("=" * 100)
    logger.info(f"\n{'AA':<5} {'PDUs':<12} {'Clusters':<12} {'Noise':<12} {'Noise %':<10} {'Mean Prob':<12}")
    logger.info("-" * 100)

    for _, row in results_df.iterrows():
        logger.info(
            f"{row['aa']:<5} {int(row['n_pdus']):>11,} {int(row['n_clusters']):>11} "
            f"{int(row['n_noise']):>11,} {row['noise_pct']:>9.1f}% {row['mean_prob']:>11.3f}"
        )

    logger.info("-" * 100)

    # Statistics
    logger.info("\n" + "=" * 100)
    logger.info("STATISTICS")
    logger.info("=" * 100)

    logger.info(f"\nCluster count across AAs:")
    logger.info(f"  Mean: {results_df['n_clusters'].mean():.0f}")
    logger.info(f"  Median: {results_df['n_clusters'].median():.0f}")
    logger.info(f"  Min: {results_df['n_clusters'].min()} (AA={results_df.loc[results_df['n_clusters'].idxmin(), 'aa']})")
    logger.info(f"  Max: {results_df['n_clusters'].max()} (AA={results_df.loc[results_df['n_clusters'].idxmax(), 'aa']})")
    logger.info(f"  Std: {results_df['n_clusters'].std():.1f}")

    logger.info(f"\nNoise percentage across AAs:")
    logger.info(f"  Mean: {results_df['noise_pct'].mean():.1f}%")
    logger.info(f"  Median: {results_df['noise_pct'].median():.1f}%")
    logger.info(f"  Min: {results_df['noise_pct'].min():.1f}% (AA={results_df.loc[results_df['noise_pct'].idxmin(), 'aa']})")
    logger.info(f"  Max: {results_df['noise_pct'].max():.1f}% (AA={results_df.loc[results_df['noise_pct'].idxmax(), 'aa']})")

    logger.info(f"\nCluster probability (mean across PDUs):")
    logger.info(f"  Mean: {results_df['mean_prob'].mean():.3f}")
    logger.info(f"  Range: {results_df['mean_prob'].min():.3f} - {results_df['mean_prob'].max():.3f}")

    # Key observations
    logger.info("\n" + "=" * 100)
    logger.info("KEY OBSERVATIONS")
    logger.info("=" * 100)

    largest_cluster = results_df.loc[results_df['n_clusters'].idxmax()]
    smallest_cluster = results_df.loc[results_df['n_clusters'].idxmin()]

    logger.info(f"\n1. DIVERSITY OF CLUSTERING:")
    logger.info(f"   - Largest: AA={largest_cluster['aa']} with {int(largest_cluster['n_clusters'])} clusters ({int(largest_cluster['n_pdus']):,} PDUs)")
    logger.info(f"   - Smallest: AA={smallest_cluster['aa']} with {int(smallest_cluster['n_clusters'])} clusters ({int(smallest_cluster['n_pdus']):,} PDUs)")
    logger.info(f"   - Ratio: {largest_cluster['n_clusters'] / smallest_cluster['n_clusters']:.1f}x difference")

    logger.info(f"\n2. CLUSTER QUALITY:")
    high_prob = (results_df['mean_prob'] > 0.95).sum()
    logger.info(f"   - High confidence (mean prob >0.95): {high_prob} / {len(results_df)} AAs")
    logger.info(f"   - Average confidence: {results_df['mean_prob'].mean():.3f}")

    logger.info(f"\n3. NOISE LEVELS:")
    low_noise = (results_df['noise_pct'] < 10).sum()
    high_noise = (results_df['noise_pct'] > 20).sum()
    logger.info(f"   - Low noise (<10%): {low_noise} AAs")
    logger.info(f"   - High noise (>20%): {high_noise} AAs")

    logger.info(f"\n4. TOTAL COVERAGE:")
    total_pdus = results_df['n_pdus'].sum()
    total_clusters = results_df['n_clusters'].sum()
    total_noise = results_df['n_noise'].sum()
    logger.info(f"   - Total PDUs: {total_pdus:,}")
    logger.info(f"   - Total clusters: {total_clusters}")
    logger.info(f"   - Total noise points: {total_noise:,} ({100*total_noise/total_pdus:.1f}%)")

    # Save summary
    summary_file = clusters_dir.parent / "cluster_summary_all_aa.csv"
    results_df.to_csv(summary_file, index=False)
    logger.info(f"\n✓ Summary saved: {summary_file}")

    logger.info("\n" + "=" * 100 + "\n")


if __name__ == "__main__":
    main()
