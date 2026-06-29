#!/usr/bin/env python3
"""
Analyze PDU cluster profiles: What does each cluster's neighborhood look like?

Run locally on your machine where the database is accessible.
"""

import sqlite3
import argparse
from pathlib import Path
from collections import Counter
import json

import numpy as np
import pandas as pd


def analyze_cluster_profiles(db_path: str, aa: str = "L", top_clusters: int = 10):
    """
    Extract neighborhood profiles for each cluster.

    For each cluster: What amino acids appear? What SS types? What distances?
    """

    conn = sqlite3.connect(db_path)

    # Get cluster assignments from the CSV (if available)
    # For now, we'll just look at clusters from clustering results

    print(f"\n{'='*70}")
    print(f"CLUSTER PROFILE ANALYSIS: {aa}")
    print(f"{'='*70}")

    # 1. Get basic stats
    pdu_count = conn.execute(
        "SELECT COUNT(*) FROM pdu WHERE reference_residue_one_letter = ?",
        (aa,)
    ).fetchone()[0]

    print(f"\nTotal PDUs: {pdu_count:,}")

    # 2. Sample a few clusters and show their neighborhood patterns
    # Get clusters with enough data
    clusters_query = """
    SELECT id FROM pdu
    WHERE reference_residue_one_letter = ?
    LIMIT 3
    """

    sample_pdu_ids = [row[0] for row in conn.execute(clusters_query, (aa,)).fetchall()]

    print(f"\nAnalyzing {len(sample_pdu_ids)} sample PDUs...\n")

    for pdu_id in sample_pdu_ids:
        # Get neighbors for this PDU
        neighbors_query = """
        SELECT residue_one_letter, secondary_structure, distance_angstrom
        FROM pdu_residue
        WHERE pdu_id = ?
        ORDER BY distance_angstrom
        """

        neighbors = conn.execute(neighbors_query, (pdu_id,)).fetchall()

        if not neighbors:
            continue

        # Analyze this PDU's neighborhood
        neighbor_aas = [n[0] for n in neighbors]
        neighbor_ss = [n[1] for n in neighbors]
        neighbor_dists = [n[2] for n in neighbors]

        aa_composition = Counter(neighbor_aas)
        ss_composition = Counter(neighbor_ss)

        close_neighbors = sum(1 for d in neighbor_dists if d < 5)
        mid_neighbors = sum(1 for d in neighbor_dists if 5 <= d < 10)
        far_neighbors = sum(1 for d in neighbor_dists if 10 <= d <= 15)

        print(f"PDU {pdu_id}:")
        print(f"  Total neighbors: {len(neighbors)}")
        print(f"  Neighbor AA composition:")
        for aa_type, count in aa_composition.most_common(5):
            pct = 100 * count / len(neighbors)
            print(f"    {aa_type}: {count} ({pct:.1f}%)")
        print(f"  Neighbor SS composition:")
        for ss_type, count in ss_composition.most_common():
            pct = 100 * count / len(neighbors)
            print(f"    {ss_type}: {count} ({pct:.1f}%)")
        print(f"  Distance distribution:")
        print(f"    Close (0-5Å): {close_neighbors} ({100*close_neighbors/len(neighbors):.1f}%)")
        print(f"    Mid (5-10Å): {mid_neighbors} ({100*mid_neighbors/len(neighbors):.1f}%)")
        print(f"    Far (10-15Å): {far_neighbors} ({100*far_neighbors/len(neighbors):.1f}%)")
        print()

    # 3. Global statistics
    print(f"{'='*70}")
    print(f"GLOBAL NEIGHBORHOOD STATISTICS")
    print(f"{'='*70}\n")

    global_query = """
    SELECT
        residue_one_letter,
        secondary_structure,
        COUNT(*) as count
    FROM pdu_residue pr
    WHERE pr.pdu_id IN (
        SELECT id FROM pdu WHERE reference_residue_one_letter = ?
    )
    GROUP BY residue_one_letter, secondary_structure
    ORDER BY count DESC
    """

    results = conn.execute(global_query, (aa,)).fetchall()

    total_neighbors = sum(r[2] for r in results)

    print("Top 15 (AA, SS) pairs:")
    for i, (aa_type, ss_type, count) in enumerate(results[:15], 1):
        pct = 100 * count / total_neighbors
        print(f"  {i}. {aa_type}-{ss_type}: {count:,} ({pct:.2f}%)")

    # 4. Distance distribution
    print(f"\nDistance distribution:")
    dist_query = """
    SELECT
        CASE
            WHEN distance_angstrom < 2 THEN '0-2'
            WHEN distance_angstrom < 3 THEN '2-3'
            WHEN distance_angstrom < 5 THEN '3-5'
            WHEN distance_angstrom < 7 THEN '5-7'
            WHEN distance_angstrom < 10 THEN '7-10'
            WHEN distance_angstrom < 12 THEN '10-12'
            WHEN distance_angstrom <= 15 THEN '12-15'
            ELSE '>15'
        END as dist_bin,
        COUNT(*) as count
    FROM pdu_residue pr
    WHERE pr.pdu_id IN (
        SELECT id FROM pdu WHERE reference_residue_one_letter = ?
    )
    GROUP BY dist_bin
    ORDER BY
        CAST(SUBSTR(dist_bin, 1, INSTR(dist_bin, '-')-1) AS FLOAT)
    """

    dist_results = conn.execute(dist_query, (aa,)).fetchall()
    total_dist = sum(r[1] for r in dist_results)

    for dist_bin, count in dist_results:
        pct = 100 * count / total_dist
        bar_width = int(pct / 2)
        bar = "█" * bar_width
        print(f"  {dist_bin:>6}Å: {count:>9,} ({pct:>5.1f}%) {bar}")

    conn.close()

    print(f"\n{'='*70}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze PDU cluster neighborhood profiles"
    )
    parser.add_argument(
        "--db",
        required=True,
        help="Path to PDU SQLite database (e.g., /path/to/pdus_L.sqlite)"
    )
    parser.add_argument(
        "--aa",
        default="L",
        help="Amino acid to analyze (default: L)"
    )

    args = parser.parse_args()

    if not Path(args.db).exists():
        print(f"Error: Database not found at {args.db}")
        exit(1)

    analyze_cluster_profiles(args.db, args.aa)
