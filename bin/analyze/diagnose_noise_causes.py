#!/usr/bin/env python3
"""Diagnose what causes high noise in clustering.

Hypotheses:
1. Structural diversity (CATH families) - too many different fold types nearby
2. Distance cutoff (15Å too large) - includes too much variability
3. Central residue properties - certain positions more noise-prone
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

plt.style.use('seaborn-v0_8-darkgrid')

AAS = list('ACDEFGHIKLMNPQRSTVWY')
EMBEDDINGS_DIR = Path('analysis/embeddings')
CLUSTERS_DIR = Path('analysis/clusters_umap')
FEATURES_DIR = Path('analysis/features')
OUTPUT_DIR = Path('analysis/noise_diagnosis')
OUTPUT_DIR.mkdir(exist_ok=True)


def analyze_noise_characteristics(aa):
    """Analyze what makes a PDU noisy vs assigned."""
    cluster_file = CLUSTERS_DIR / f'pdu_clusters_{aa}_umap.csv'
    embedding_file = EMBEDDINGS_DIR / f'pdu_embedding_{aa}.npz'
    features_file = FEATURES_DIR / f'pdu_features_{aa}.npz'

    if not cluster_file.exists():
        return None

    logger.info(f"Analyzing {aa}...")

    # Load data
    clusters_df = pd.read_csv(cluster_file)
    embedding_data = np.load(embedding_file)

    try:
        features_data = np.load(features_file)
        features_Z = features_data['Z']  # 900-dim features
        has_features = True
    except:
        has_features = False
        logger.warning(f"  No feature file for {aa}")
        features_Z = None

    embeddings_Z = embedding_data['Z']  # 16-dim latent
    pdu_ids = embedding_data['pdu_ids']

    # Separate assigned vs noise
    assigned_mask = clusters_df['cluster'] != -1
    noise_mask = ~assigned_mask

    # Match PDU IDs
    pdu_to_idx = {pdu_id: idx for idx, pdu_id in enumerate(pdu_ids)}
    match_idx = np.array([pdu_to_idx.get(pid, -1) for pid in clusters_df['pdu_id']])
    valid = match_idx >= 0

    match_idx = match_idx[valid]
    assigned_mask_valid = assigned_mask[valid]
    noise_mask_valid = noise_mask[valid]

    # Metrics on embeddings (16D latent space)
    logger.info(f"  Embeddings analysis...")

    assigned_emb = embeddings_Z[match_idx[assigned_mask_valid]]
    noise_emb = embeddings_Z[match_idx[noise_mask_valid]]

    # Distance to cluster center
    if len(assigned_emb) > 0:
        assigned_mean = assigned_emb.mean(axis=0)
        assigned_dist = np.linalg.norm(assigned_emb - assigned_mean, axis=1)
        assigned_mean_dist = assigned_dist.mean()
        assigned_std_dist = assigned_dist.std()
    else:
        assigned_mean_dist = assigned_std_dist = 0

    # Noise spread
    if len(noise_emb) > 0:
        noise_mean = noise_emb.mean(axis=0)
        noise_dist = np.linalg.norm(noise_emb - noise_mean, axis=1)
        noise_mean_dist = noise_dist.mean()
        noise_std_dist = noise_dist.std()
    else:
        noise_mean_dist = noise_std_dist = 0

    # Feature-space analysis (900D)
    feature_stats = None
    if has_features:
        logger.info(f"  Feature space analysis...")
        assigned_feat = features_Z[match_idx[assigned_mask_valid]]
        noise_feat = features_Z[match_idx[noise_mask_valid]]

        # Amino acid composition diversity (first 20 dims)
        assigned_aa_div = assigned_feat[:, :20].std(axis=1).mean()
        noise_aa_div = noise_feat[:, :20].std(axis=1).mean()

        # Secondary structure diversity (next 60 dims = 20 AAs × 3 SS types)
        assigned_ss_div = assigned_feat[:, 20:80].std(axis=1).mean()
        noise_ss_div = noise_feat[:, 20:80].std(axis=1).mean()

        # Distance shell diversity (remaining dims)
        assigned_dist_div = assigned_feat[:, 80:].std(axis=1).mean()
        noise_dist_div = noise_feat[:, 80:].std(axis=1).mean()

        feature_stats = {
            'aa_diversity_assigned': assigned_aa_div,
            'aa_diversity_noise': noise_aa_div,
            'ss_diversity_assigned': assigned_ss_div,
            'ss_diversity_noise': noise_ss_div,
            'dist_diversity_assigned': assigned_dist_div,
            'dist_diversity_noise': noise_dist_div,
        }

    return {
        'aa': aa,
        'n_assigned': assigned_mask_valid.sum(),
        'n_noise': noise_mask_valid.sum(),
        'noise_pct': 100 * noise_mask_valid.sum() / len(valid),
        'assigned_mean_dist': assigned_mean_dist,
        'assigned_std_dist': assigned_std_dist,
        'noise_mean_dist': noise_mean_dist,
        'noise_std_dist': noise_std_dist,
        'dist_ratio': noise_mean_dist / assigned_mean_dist if assigned_mean_dist > 0 else 0,
        **(feature_stats or {})
    }


def plot_embedding_space_analysis():
    """Plot assigned vs noise in embedding space."""
    stats = []
    for aa in AAS:
        result = analyze_noise_characteristics(aa)
        if result:
            stats.append(result)

    df = pd.DataFrame(stats)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Distance to centroid
    axes[0, 0].scatter(df['assigned_mean_dist'], df['noise_mean_dist'], s=100, alpha=0.6)
    axes[0, 0].plot([0, df[['assigned_mean_dist', 'noise_mean_dist']].max().max()],
                     [0, df[['assigned_mean_dist', 'noise_mean_dist']].max().max()],
                     'r--', alpha=0.5, label='Equal')
    axes[0, 0].set_xlabel('Assigned PDUs: Mean Distance to Centroid')
    axes[0, 0].set_ylabel('Noise PDUs: Mean Distance to Centroid')
    axes[0, 0].set_title('Embedding Space Spread')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # Distance ratio
    axes[0, 1].bar(df['aa'], df['dist_ratio'], color='coral', alpha=0.7)
    axes[0, 1].axhline(y=1, color='r', linestyle='--', alpha=0.5, label='Equal')
    axes[0, 1].set_title('Noise Spread vs Assigned Spread Ratio')
    axes[0, 1].set_ylabel('Noise Distance / Assigned Distance')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3, axis='y')

    # Noise percentage
    axes[1, 0].bar(df['aa'], df['noise_pct'], color='steelblue', alpha=0.7)
    axes[1, 0].set_title('Noise Percentage')
    axes[1, 0].set_ylabel('% Noise')
    axes[1, 0].set_ylim([0, 100])
    axes[1, 0].grid(True, alpha=0.3, axis='y')

    # Amino acid diversity (if available)
    if 'aa_diversity_assigned' in df.columns:
        div_cols = ['aa_diversity_assigned', 'aa_diversity_noise']
        if all(col in df.columns for col in div_cols):
            x = np.arange(len(df))
            width = 0.35
            axes[1, 1].bar(x - width/2, df['aa_diversity_assigned'], width,
                          label='Assigned', alpha=0.7)
            axes[1, 1].bar(x + width/2, df['aa_diversity_noise'], width,
                          label='Noise', alpha=0.7)
            axes[1, 1].set_xticks(x)
            axes[1, 1].set_xticklabels(df['aa'])
            axes[1, 1].set_title('Structural Diversity (AA types)')
            axes[1, 1].set_ylabel('Mean Feature Std Dev')
            axes[1, 1].legend()
            axes[1, 1].grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'noise_analysis.png', dpi=300, bbox_inches='tight')
    logger.info("✓ Saved: noise_analysis.png")
    plt.close()

    # Summary table
    summary = df[['aa', 'noise_pct', 'assigned_mean_dist', 'noise_mean_dist', 'dist_ratio']].copy()
    summary.columns = ['AA', 'Noise %', 'Assigned Dist', 'Noise Dist', 'Ratio']
    summary = summary.round(3)

    logger.info("\n" + "="*80)
    logger.info("NOISE ANALYSIS SUMMARY")
    logger.info("="*80)
    logger.info(summary.to_string(index=False))
    logger.info("="*80)

    summary.to_csv(OUTPUT_DIR / 'noise_analysis.csv', index=False)


def main():
    logger.info("="*80)
    logger.info("Diagnosing Noise Causes")
    logger.info("="*80)
    logger.info("\nHypotheses:")
    logger.info("1. Structural diversity (CATH) - too many fold types near central AA")
    logger.info("2. Distance cutoff (15Å) - too broad, includes too much variability")
    logger.info("3. Central AA properties - certain positions inherently noisier")
    logger.info("\n")

    plot_embedding_space_analysis()

    logger.info("\n" + "="*80)
    logger.info("Key Insights:")
    logger.info("="*80)
    logger.info("→ Compare 'Assigned Dist' vs 'Noise Dist':")
    logger.info("  - If similar: noise = sparse regions (good)")
    logger.info("  - If noise >> assigned: too much variability (need tighter cutoff)")
    logger.info("\n→ Check 'Ratio' column:")
    logger.info("  - Ratio > 1.5: noise is MORE spread out (suggests 15Å is too large)")
    logger.info("  - Ratio ≈ 1.0: noise = natural sparsity (15Å cutoff OK)")
    logger.info("="*80)


if __name__ == "__main__":
    main()
