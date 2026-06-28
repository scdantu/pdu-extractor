#!/usr/bin/env python3
"""
Validate PDU clusters against central residue's own secondary structure.

Question: Does the secondary structure of the CENTRAL residue predict
which cluster its local environment belongs to?

This proves clusters capture real structural context, not just neighbor stats.
"""

import argparse
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt
import seaborn as sns


def get_central_ss(db_path: str, pdu_ids: list) -> dict:
    """
    Query database for each PDU's central residue secondary structure.

    The central residue is at distance_angstrom ≈ 0 (the query residue itself).
    """
    conn = sqlite3.connect(db_path)

    # Get central residue SS for each PDU
    # The central residue is the one at distance ~0 from itself
    placeholders = ",".join("?" for _ in pdu_ids)
    query = f"""
    SELECT pr.pdu_id, pr.secondary_structure
    FROM pdu_residue pr
    WHERE pr.pdu_id IN ({placeholders})
      AND pr.distance_angstrom < 0.5
    """

    result = conn.execute(query, pdu_ids).fetchall()
    conn.close()

    # Build dict: pdu_id → secondary_structure
    central_ss_map = {pdu_id: ss for pdu_id, ss in result}

    return central_ss_map


def main():
    parser = argparse.ArgumentParser(
        description="Validate clusters against central residue secondary structure"
    )
    parser.add_argument("--clusters-csv", required=True, help="Cluster CSV from HDBSCAN")
    parser.add_argument("--db", required=True, help="PDU SQLite database")
    parser.add_argument("--aa", required=True, help="Amino acid (e.g., L)")
    parser.add_argument("--out-dir", default="analysis/validation", help="Output directory")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Read cluster assignments
    clusters_df = pd.read_csv(args.clusters_csv)
    print(f"Loaded {len(clusters_df)} cluster assignments from {args.clusters_csv}")

    # Get central SS for all PDUs
    pdu_ids = clusters_df['pdu_id'].tolist()
    central_ss_map = get_central_ss(args.db, pdu_ids)
    print(f"Retrieved central SS for {len(central_ss_map)} PDUs")

    # Map central SS to clusters
    clusters_df['central_ss'] = clusters_df['pdu_id'].map(central_ss_map)
    clusters_df = clusters_df.dropna(subset=['central_ss'])

    print(f"Mapped to {len(clusters_df)} PDUs with central SS data")

    # Filter out noise points (-1)
    assigned_df = clusters_df[clusters_df['cluster'] != -1]
    noise_df = clusters_df[clusters_df['cluster'] == -1]

    print(f"Assigned: {len(assigned_df)}, Noise: {len(noise_df)}")

    # Build contingency table
    contingency = pd.crosstab(
        assigned_df['central_ss'],
        assigned_df['cluster'],
        margins=False
    )

    print("\n" + "=" * 70)
    print(f"Contingency Table: Central SS × Cluster (AA={args.aa})")
    print("=" * 70)
    print(f"Shape: {contingency.shape[0]} SS types × {contingency.shape[1]} clusters")

    # Chi-square test
    chi2, p_value, dof, expected = chi2_contingency(contingency.values)

    print(f"\nChi-Square Test:")
    print(f"  χ² statistic: {chi2:.2f}")
    print(f"  p-value: {p_value:.2e}")
    print(f"  DOF: {dof}")

    if p_value < 0.001:
        print(f"  ✓ HIGHLY SIGNIFICANT: Central SS strongly predicts cluster")
    elif p_value < 0.05:
        print(f"  ✓ SIGNIFICANT: Central SS predicts cluster")
    else:
        print(f"  ✗ NOT SIGNIFICANT: Central SS does not predict cluster")

    # Cramér's V (effect size)
    n = contingency.values.sum()
    min_dim = min(contingency.shape[0] - 1, contingency.shape[1] - 1)
    cramers_v = np.sqrt(chi2 / (n * min_dim)) if min_dim > 0 else 0
    print(f"\nCramér's V (effect size): {cramers_v:.3f}")
    if cramers_v < 0.1:
        print(f"  Weak effect")
    elif cramers_v < 0.3:
        print(f"  Moderate effect")
    else:
        print(f"  Strong effect")

    # Distribution per SS type
    print(f"\nDistribution by Central SS:")
    for ss_type in sorted(assigned_df['central_ss'].unique()):
        count = (assigned_df['central_ss'] == ss_type).sum()
        pct = 100 * count / len(assigned_df)
        print(f"  {ss_type}: {count:,} PDUs ({pct:.1f}%)")

    # Save results
    results = {
        'aa': args.aa,
        'n_pdus': len(assigned_df),
        'n_clusters': contingency.shape[1],
        'chi2': chi2,
        'p_value': p_value,
        'dof': dof,
        'cramers_v': cramers_v,
        'significant': p_value < 0.05
    }

    results_df = pd.DataFrame([results])
    results_df.to_csv(out_dir / f"validation_central_ss_{args.aa}.csv", index=False)
    print(f"\nResults saved to {out_dir / f'validation_central_ss_{args.aa}.csv'}")

    # Save contingency table
    contingency.to_csv(out_dir / f"contingency_ss_cluster_{args.aa}.csv")
    print(f"Contingency table saved to {out_dir / f'contingency_ss_cluster_{args.aa}.csv'}")

    # Visualize
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    # Heatmap
    contingency_norm = contingency.div(contingency.sum(axis=1), axis=0)
    sns.heatmap(contingency_norm, cmap='YlOrRd', ax=axes[0], cbar_kws={'label': 'Fraction'})
    axes[0].set_title(f'AA={args.aa}: Central SS → Cluster Assignment\n(χ²={chi2:.1f}, p={p_value:.2e}, V={cramers_v:.3f})')
    axes[0].set_xlabel('Cluster ID')
    axes[0].set_ylabel('Central Secondary Structure')

    # Count by SS type
    ss_counts = assigned_df['central_ss'].value_counts().sort_index()
    axes[1].bar(ss_counts.index, ss_counts.values, color='steelblue')
    axes[1].set_title(f'PDU Distribution by Central SS (AA={args.aa})')
    axes[1].set_xlabel('Central Secondary Structure')
    axes[1].set_ylabel('Count')
    axes[1].set_yscale('log')

    plt.tight_layout()
    plt.savefig(out_dir / f"validation_central_ss_{args.aa}.png", dpi=150, bbox_inches='tight')
    print(f"Figure saved to {out_dir / f'validation_central_ss_{args.aa}.png'}")

    print("\n✓ Validation complete")


if __name__ == "__main__":
    main()
