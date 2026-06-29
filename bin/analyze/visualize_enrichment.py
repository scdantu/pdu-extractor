#!/usr/bin/env python3
"""Visualize PyCoM enrichment results across all amino acids."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

AAS = list('ACDEFGHIKLMNPQRSTVWY')
ENRICHMENT_DIR = Path('analysis/enrichment')


def load_enrichment_stats():
    """Extract enrichment statistics from contingency tables."""
    stats = []

    for aa in AAS:
        cath_file = ENRICHMENT_DIR / f'contingency_cath_{aa}.csv'
        ec_file = ENRICHMENT_DIR / f'contingency_ec_{aa}.csv'

        if not cath_file.exists() or not ec_file.exists():
            logger.warning(f"Missing enrichment files for {aa}")
            continue

        # Load contingency tables
        cath_contingency = pd.read_csv(cath_file, index_col=0)
        ec_contingency = pd.read_csv(ec_file, index_col=0)

        # Compute chi-square manually
        def chi2_stat(contingency):
            contingency = contingency.astype(float).values
            row_sums = contingency.sum(axis=1, keepdims=True)
            col_sums = contingency.sum(axis=0, keepdims=True)
            total = contingency.sum()
            expected = (row_sums * col_sums) / total
            expected[expected == 0] = 1
            chi2 = ((contingency - expected) ** 2 / expected).sum()
            return chi2

        stats.append({
            'aa': aa,
            'cath_chi2': chi2_stat(cath_contingency),
            'ec_chi2': chi2_stat(ec_contingency),
            'n_cath_classes': len(cath_contingency),
            'n_ec_classes': len(ec_contingency),
            'n_clusters': cath_contingency.shape[1]
        })

    return pd.DataFrame(stats)


def plot_chi2_comparison():
    """Compare chi-square values across AAs."""
    stats = load_enrichment_stats()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # CATH enrichment
    ax1.bar(stats['aa'], stats['cath_chi2'], color='steelblue', alpha=0.7)
    ax1.set_title('CATH Structural Enrichment (χ²)', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Amino Acid')
    ax1.set_ylabel('χ² Statistic')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3)

    # EC enzymatic enrichment
    ax2.bar(stats['aa'], stats['ec_chi2'], color='coral', alpha=0.7)
    ax2.set_title('EC Enzymatic Enrichment (χ²)', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Amino Acid')
    ax2.set_ylabel('χ² Statistic')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(ENRICHMENT_DIR / 'chi2_enrichment_comparison.png', dpi=300, bbox_inches='tight')
    logger.info(f"✓ Saved: chi2_enrichment_comparison.png")
    plt.close()


def plot_enrichment_heatmap():
    """Create heatmap of enrichment statistics."""
    stats = load_enrichment_stats()

    # Normalize chi-square values
    stats['cath_chi2_norm'] = stats['cath_chi2'] / stats['cath_chi2'].max()
    stats['ec_chi2_norm'] = stats['ec_chi2'] / stats['ec_chi2'].max()

    data = stats[['aa', 'cath_chi2_norm', 'ec_chi2_norm', 'n_clusters']].set_index('aa')
    data.columns = ['CATH Enrichment', 'EC Enrichment', 'N Clusters']

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(data.T, annot=True, fmt='.2f', cmap='YlOrRd', cbar_kws={'label': 'Normalized χ²'}, ax=ax)
    ax.set_title('Enrichment Signal Across Amino Acids', fontsize=14, fontweight='bold')
    ax.set_xlabel('Amino Acid')

    plt.tight_layout()
    plt.savefig(ENRICHMENT_DIR / 'enrichment_heatmap.png', dpi=300, bbox_inches='tight')
    logger.info(f"✓ Saved: enrichment_heatmap.png")
    plt.close()


def plot_cluster_diversity():
    """Show cluster diversity per AA."""
    stats = load_enrichment_stats()

    fig, ax = plt.subplots(figsize=(12, 5))

    x = np.arange(len(stats))
    width = 0.35

    ax.bar(x - width/2, stats['n_cath_classes'], width, label='CATH Classes', alpha=0.8)
    ax.bar(x + width/2, stats['n_clusters'], width, label='PDU Clusters', alpha=0.8)

    ax.set_xlabel('Amino Acid', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Structural Diversity: CATH Classes vs PDU Clusters', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(stats['aa'])
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(ENRICHMENT_DIR / 'cluster_diversity.png', dpi=300, bbox_inches='tight')
    logger.info(f"✓ Saved: cluster_diversity.png")
    plt.close()


def plot_top_cath_ec():
    """Plot top CATH and EC classes for selected AAs."""
    selected_aas = ['L', 'A', 'G', 'V', 'P']  # Most common/interesting

    fig, axes = plt.subplots(len(selected_aas), 2, figsize=(14, 3*len(selected_aas)))

    for idx, aa in enumerate(selected_aas):
        cath_file = ENRICHMENT_DIR / f'contingency_cath_{aa}.csv'
        ec_file = ENRICHMENT_DIR / f'contingency_ec_{aa}.csv'

        if not cath_file.exists() or not ec_file.exists():
            continue

        # CATH
        cath_contingency = pd.read_csv(cath_file, index_col=0)
        cath_totals = cath_contingency.sum(axis=1).sort_values(ascending=False)[:10]

        axes[idx, 0].barh(range(len(cath_totals)), cath_totals.values, color='steelblue', alpha=0.7)
        axes[idx, 0].set_yticks(range(len(cath_totals)))
        axes[idx, 0].set_yticklabels(cath_totals.index, fontsize=9)
        axes[idx, 0].set_xlabel('PDU Count')
        axes[idx, 0].set_title(f'{aa}: Top CATH Classes')
        axes[idx, 0].invert_yaxis()

        # EC
        ec_contingency = pd.read_csv(ec_file, index_col=0)
        ec_totals = ec_contingency.sum(axis=1).sort_values(ascending=False)[:10]

        axes[idx, 1].barh(range(len(ec_totals)), ec_totals.values, color='coral', alpha=0.7)
        axes[idx, 1].set_yticks(range(len(ec_totals)))
        axes[idx, 1].set_yticklabels(ec_totals.index, fontsize=9)
        axes[idx, 1].set_xlabel('PDU Count')
        axes[idx, 1].set_title(f'{aa}: Top EC Numbers')
        axes[idx, 1].invert_yaxis()

    plt.tight_layout()
    plt.savefig(ENRICHMENT_DIR / 'top_cath_ec_classes.png', dpi=300, bbox_inches='tight')
    logger.info(f"✓ Saved: top_cath_ec_classes.png")
    plt.close()


def generate_summary_table():
    """Generate summary statistics table."""
    stats = load_enrichment_stats()

    summary = pd.DataFrame({
        'AA': stats['aa'],
        'CATH χ²': stats['cath_chi2'].astype(int),
        'EC χ²': stats['ec_chi2'].astype(int),
        'CATH Classes': stats['n_cath_classes'],
        'EC Classes': stats['n_ec_classes'],
        'PDU Clusters': stats['n_clusters']
    })

    # Save to CSV
    summary.to_csv(ENRICHMENT_DIR / 'enrichment_summary.csv', index=False)
    logger.info(f"✓ Saved: enrichment_summary.csv")

    # Print table
    logger.info("\n" + "="*80)
    logger.info("ENRICHMENT SUMMARY")
    logger.info("="*80)
    logger.info(summary.to_string(index=False))
    logger.info("="*80)

    return summary


def main():
    logger.info("="*80)
    logger.info("Visualizing PyCoM Enrichment Results")
    logger.info("="*80)

    if not ENRICHMENT_DIR.exists():
        logger.error(f"Enrichment directory not found: {ENRICHMENT_DIR}")
        return

    logger.info(f"\n[1/5] Generating summary table...")
    summary = generate_summary_table()

    logger.info(f"\n[2/5] Plotting chi-square comparison...")
    plot_chi2_comparison()

    logger.info(f"\n[3/5] Creating enrichment heatmap...")
    plot_enrichment_heatmap()

    logger.info(f"\n[4/5] Plotting cluster diversity...")
    plot_cluster_diversity()

    logger.info(f"\n[5/5] Showing top CATH/EC classes...")
    plot_top_cath_ec()

    logger.info("\n" + "="*80)
    logger.info("✓ Visualization complete!")
    logger.info("="*80)
    logger.info(f"\nOutputs saved to: {ENRICHMENT_DIR}")
    logger.info(f"\nGenerated files:")
    logger.info(f"  - enrichment_summary.csv")
    logger.info(f"  - chi2_enrichment_comparison.png")
    logger.info(f"  - enrichment_heatmap.png")
    logger.info(f"  - cluster_diversity.png")
    logger.info(f"  - top_cath_ec_classes.png")


if __name__ == "__main__":
    main()
