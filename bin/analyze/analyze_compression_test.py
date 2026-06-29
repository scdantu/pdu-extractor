#!/usr/bin/env python3
"""Analyze compression hypothesis test results.

Compare 10Å vs 15Å across latent dimensions 8-256D
to determine if compression is the bottleneck.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

plt.style.use('seaborn-v0_8-darkgrid')

def analyze_results(summary_file):
    """Load and analyze compression test results."""

    if not Path(summary_file).exists():
        logger.error(f"Summary file not found: {summary_file}")
        return None

    df = pd.read_csv(summary_file)

    logger.info("="*80)
    logger.info("COMPRESSION HYPOTHESIS TEST RESULTS")
    logger.info("="*80)
    logger.info(f"\nDataset: {Path(summary_file).stem}")
    logger.info(f"\nResults:\n{df.to_string(index=False)}\n")

    # Analysis
    logger.info("="*80)
    logger.info("KEY FINDINGS")
    logger.info("="*80)

    # Compare by distance
    logger.info("\n1. DISTANCE COMPARISON (10Å vs 15Å):")
    for dist in df['Distance'].unique():
        subset = df[df['Distance'] == dist]
        logger.info(f"\n  {dist}Å:")
        logger.info(f"    Avg noise: {subset['NoisePercent'].mean():.1f}%")
        logger.info(f"    Avg confidence: {subset['AvgConf'].mean():.3f}")
        logger.info(f"    Best dim for noise: {subset.loc[subset['NoisePercent'].idxmin(), 'LatentDim']:.0f}D "
                   f"({subset['NoisePercent'].min():.1f}%)")

    # Compare by dimension
    logger.info("\n2. DIMENSION ANALYSIS (compression ratio):")
    logger.info(f"\n  900D → XD compression:")
    for dim in sorted(df['LatentDim'].unique()):
        subset = df[df['LatentDim'] == dim]
        compression = 900 / dim
        logger.info(f"    {dim:3d}D (compress {compression:6.1f}x): "
                   f"avg noise {subset['NoisePercent'].mean():.1f}%, "
                   f"conf {subset['AvgConf'].mean():.3f}")

    # Find optimal
    logger.info("\n3. HYPOTHESIS TEST:")
    best_row = df.loc[df['NoisePercent'].idxmin()]
    logger.info(f"\n  Best result: {best_row['Distance']:.0f}Å × {best_row['LatentDim']:.0f}D")
    logger.info(f"    Noise: {best_row['NoisePercent']:.1f}%")
    logger.info(f"    Confidence: {best_row['AvgConf']:.3f}")

    # Is compression the problem?
    improvement_8d = df[df['LatentDim'] == 8]['NoisePercent'].mean()
    improvement_256d = df[df['LatentDim'] == 256]['NoisePercent'].mean()

    logger.info(f"\n  Compression impact:")
    logger.info(f"    8D  (112.5x compression): {improvement_8d:.1f}% noise")
    logger.info(f"    256D (3.5x compression):  {improvement_256d:.1f}% noise")
    logger.info(f"    Difference: {abs(improvement_8d - improvement_256d):.1f}%")

    if abs(improvement_8d - improvement_256d) < 10:
        logger.info(f"\n  → CONCLUSION: Compression is NOT the problem")
        logger.info(f"    (Noise varies only {abs(improvement_8d - improvement_256d):.1f}% across dims)")
    else:
        logger.info(f"\n  → CONCLUSION: Compression IS significant")
        logger.info(f"    (Higher dims reduce noise by {abs(improvement_8d - improvement_256d):.1f}%)")

    # Distance impact
    dist_10_noise = df[df['Distance'] == 10]['NoisePercent'].mean()
    dist_15_noise = df[df['Distance'] == 15]['NoisePercent'].mean()

    logger.info(f"\n  Distance cutoff impact:")
    logger.info(f"    10Å: {dist_10_noise:.1f}% noise")
    logger.info(f"    15Å: {dist_15_noise:.1f}% noise")
    logger.info(f"    Difference: {abs(dist_10_noise - dist_15_noise):.1f}%")

    if abs(dist_10_noise - dist_15_noise) > 10:
        logger.info(f"\n  → CONCLUSION: Distance cutoff IS the problem")
        logger.info(f"    (10Å reduces noise by {abs(dist_10_noise - dist_15_noise):.1f}%)")
    else:
        logger.info(f"\n  → CONCLUSION: Distance cutoff is NOT critical")
        logger.info(f"    (Noise varies only {abs(dist_10_noise - dist_15_noise):.1f}%)")

    logger.info("="*80)

    return df


def plot_results(df, output_file):
    """Create comparison visualizations."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Noise vs dimension (by distance)
    ax = axes[0, 0]
    for dist in sorted(df['Distance'].unique()):
        subset = df[df['Distance'] == dist].sort_values('LatentDim')
        ax.plot(subset['LatentDim'], subset['NoisePercent'],
               marker='o', label=f"{dist}Å", linewidth=2)
    ax.set_xlabel('Latent Dimension')
    ax.set_ylabel('Noise %')
    ax.set_title('Noise vs Compression Ratio')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Confidence vs dimension (by distance)
    ax = axes[0, 1]
    for dist in sorted(df['Distance'].unique()):
        subset = df[df['Distance'] == dist].sort_values('LatentDim')
        ax.plot(subset['LatentDim'], subset['AvgConf'],
               marker='o', label=f"{dist}Å", linewidth=2)
    ax.set_xlabel('Latent Dimension')
    ax.set_ylabel('Average Confidence')
    ax.set_title('Cluster Quality vs Compression')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Heatmap: noise
    ax = axes[1, 0]
    pivot_noise = df.pivot(index='Distance', columns='LatentDim', values='NoisePercent')
    sns.heatmap(pivot_noise, annot=True, fmt='.1f', cmap='RdYlGn_r', ax=ax, cbar_kws={'label': 'Noise %'})
    ax.set_title('Noise Percentage Heatmap')

    # Heatmap: confidence
    ax = axes[1, 1]
    pivot_conf = df.pivot(index='Distance', columns='LatentDim', values='AvgConf')
    sns.heatmap(pivot_conf, annot=True, fmt='.3f', cmap='RdYlGn', ax=ax, cbar_kws={'label': 'Avg Confidence'})
    ax.set_title('Cluster Confidence Heatmap')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    logger.info(f"\n✓ Saved: {output_file}")
    plt.close()


if __name__ == "__main__":
    import sys

    aa = sys.argv[1] if len(sys.argv) > 1 else "L"
    summary_file = Path(f"analysis/compression_test/summary_{aa}.csv")

    df = analyze_results(str(summary_file))

    if df is not None:
        plot_file = summary_file.parent / f"compression_analysis_{aa}.png"
        plot_results(df, str(plot_file))
