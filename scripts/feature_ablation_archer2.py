#!/usr/bin/env python3
"""
Feature ablation study: Extract different feature subsets from SQLite.

Tests which features matter most for clustering quality.

Feature combinations:
  - full: All 900 features (20 AA × 3 SS × 15 distance shells)
  - aa_only: 20 features (amino acid composition)
  - ss_only: 3 features (secondary structure)
  - dist_only: 16 features (distance distribution)
  - aa_ss: 20+3 = 23 features
  - aa_dist: 20+16 = 36 features
  - ss_dist: 3+16 = 19 features
  - custom: User-specified subset

Usage:
  python3 feature_ablation_archer2.py \
    --db /path/to/pdus_L.sqlite \
    --aa L \
    --subset full \
    --out-dir analysis/ablation
"""

import sqlite3
import argparse
from pathlib import Path
from collections import defaultdict
import logging

import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class FeatureExtractor:
    """Extract feature subsets from PDU SQLite database."""

    def __init__(self, db_path, aa):
        self.db_path = db_path
        self.aa = aa
        self.conn = sqlite3.connect(db_path)
        self.cursor = self.conn.cursor()

        # Get metadata
        self.cursor.execute("SELECT COUNT(*) FROM pdu WHERE reference_residue_one_letter = ?", (aa,))
        self.n_pdus = self.cursor.fetchone()[0]

        logger.info(f"Database: {db_path}")
        logger.info(f"AA: {aa}, PDUs: {self.n_pdus:,}")

    def extract_full(self):
        """Extract all 900 features (20 AA × 3 SS × 15 distance bins)."""
        logger.info("Extracting full 900-dim features...")

        features = np.zeros((self.n_pdus, 900), dtype=np.float32)
        pdu_ids = np.zeros(self.n_pdus, dtype=np.int64)

        # 20 AAs
        aas = 'ACDEFGHIKLMNPQRSTVWY'
        ss_types = ['H', 'E', 'C']  # Helix, Extended, Coil
        distance_bins = [(0, 2), (2, 3), (3, 5), (5, 7), (7, 10), (10, 12), (12, 15)]
        distance_bins += [(0, 15)]  # Also include 0-15Å aggregate

        # Query all neighbors for this AA
        self.cursor.execute(
            "SELECT pr.pdu_id, pr.residue_one_letter, pr.secondary_structure, pr.distance_angstrom "
            "FROM pdu_residue pr "
            "INNER JOIN pdu p ON pr.pdu_id = p.id "
            "WHERE p.reference_residue_one_letter = ? "
            "ORDER BY pr.pdu_id",
            (self.aa,)
        )

        current_pdu_idx = 0
        current_pdu_id = None
        pdu_neighbors = defaultdict(list)

        for pdu_id, neighbor_aa, ss, dist in self.cursor.fetchall():
            pdu_neighbors[pdu_id].append((neighbor_aa, ss, dist))

        # Build feature matrix
        for pdu_idx, (pdu_id, neighbors) in enumerate(sorted(pdu_neighbors.items())):
            pdu_ids[pdu_idx] = pdu_id

            feature_idx = 0
            for aa_idx, aa_letter in enumerate(aas):
                for ss_idx, ss_type in enumerate(ss_types):
                    for dist_min, dist_max in distance_bins:
                        count = sum(
                            1 for n_aa, n_ss, n_dist in neighbors
                            if n_aa == aa_letter and n_ss == ss_type
                            and dist_min <= n_dist < dist_max
                        )
                        features[pdu_idx, feature_idx] = count
                        feature_idx += 1

        logger.info(f"Shape: {features.shape}")
        return features, pdu_ids

    def extract_aa_only(self):
        """Extract only amino acid composition (20 features)."""
        logger.info("Extracting AA-only features...")

        features = np.zeros((self.n_pdus, 20), dtype=np.float32)
        pdu_ids = np.zeros(self.n_pdus, dtype=np.int64)
        aas = 'ACDEFGHIKLMNPQRSTVWY'

        self.cursor.execute(
            "SELECT pr.pdu_id, pr.residue_one_letter "
            "FROM pdu_residue pr "
            "INNER JOIN pdu p ON pr.pdu_id = p.id "
            "WHERE p.reference_residue_one_letter = ? "
            "ORDER BY pr.pdu_id",
            (self.aa,)
        )

        pdu_neighbors = defaultdict(list)
        for pdu_id, neighbor_aa in self.cursor.fetchall():
            pdu_neighbors[pdu_id].append(neighbor_aa)

        for pdu_idx, (pdu_id, neighbors) in enumerate(sorted(pdu_neighbors.items())):
            pdu_ids[pdu_idx] = pdu_id
            for aa_idx, aa_letter in enumerate(aas):
                features[pdu_idx, aa_idx] = neighbors.count(aa_letter)

        logger.info(f"Shape: {features.shape}")
        return features, pdu_ids

    def extract_ss_only(self):
        """Extract only secondary structure (3 features)."""
        logger.info("Extracting SS-only features...")

        features = np.zeros((self.n_pdus, 3), dtype=np.float32)
        pdu_ids = np.zeros(self.n_pdus, dtype=np.int64)
        ss_types = ['H', 'E', 'C']

        self.cursor.execute(
            "SELECT pr.pdu_id, pr.secondary_structure "
            "FROM pdu_residue pr "
            "INNER JOIN pdu p ON pr.pdu_id = p.id "
            "WHERE p.reference_residue_one_letter = ? "
            "ORDER BY pr.pdu_id",
            (self.aa,)
        )

        pdu_neighbors = defaultdict(list)
        for pdu_id, ss in self.cursor.fetchall():
            pdu_neighbors[pdu_id].append(ss)

        for pdu_idx, (pdu_id, neighbors) in enumerate(sorted(pdu_neighbors.items())):
            pdu_ids[pdu_idx] = pdu_id
            for ss_idx, ss_type in enumerate(ss_types):
                features[pdu_idx, ss_idx] = neighbors.count(ss_type)

        logger.info(f"Shape: {features.shape}")
        return features, pdu_ids

    def extract_dist_only(self):
        """Extract only distance distribution (16 features)."""
        logger.info("Extracting distance-only features...")

        features = np.zeros((self.n_pdus, 16), dtype=np.float32)
        pdu_ids = np.zeros(self.n_pdus, dtype=np.int64)
        distance_bins = [(0, 2), (2, 3), (3, 5), (5, 7), (7, 10), (10, 12), (12, 15)]
        distance_bins += [(0, 15)]  # Aggregate

        self.cursor.execute(
            "SELECT pr.pdu_id, pr.distance_angstrom "
            "FROM pdu_residue pr "
            "INNER JOIN pdu p ON pr.pdu_id = p.id "
            "WHERE p.reference_residue_one_letter = ? "
            "ORDER BY pr.pdu_id",
            (self.aa,)
        )

        pdu_neighbors = defaultdict(list)
        for pdu_id, dist in self.cursor.fetchall():
            pdu_neighbors[pdu_id].append(dist)

        for pdu_idx, (pdu_id, neighbors) in enumerate(sorted(pdu_neighbors.items())):
            pdu_ids[pdu_idx] = pdu_id
            for dist_idx, (dist_min, dist_max) in enumerate(distance_bins):
                count = sum(1 for d in neighbors if dist_min <= d < dist_max)
                features[pdu_idx, dist_idx] = count

        logger.info(f"Shape: {features.shape}")
        return features, pdu_ids

    def close(self):
        self.conn.close()


def main():
    parser = argparse.ArgumentParser(description="Feature ablation: extract subsets")
    parser.add_argument("--db", required=True, help="SQLite database path")
    parser.add_argument("--aa", required=True, help="Amino acid")
    parser.add_argument(
        "--subset",
        default="full",
        choices=["full", "aa_only", "ss_only", "dist_only"],
        help="Feature subset"
    )
    parser.add_argument("--out-dir", default="analysis/ablation", help="Output directory")

    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info(f"Feature Ablation: {args.subset}")
    logger.info("=" * 80)

    extractor = FeatureExtractor(args.db, args.aa)

    # Extract features
    if args.subset == "full":
        features, pdu_ids = extractor.extract_full()
    elif args.subset == "aa_only":
        features, pdu_ids = extractor.extract_aa_only()
    elif args.subset == "ss_only":
        features, pdu_ids = extractor.extract_ss_only()
    elif args.subset == "dist_only":
        features, pdu_ids = extractor.extract_dist_only()

    extractor.close()

    # Save
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / f"pdu_features_{args.aa}_{args.subset}.npz"
    np.savez_compressed(
        out_file,
        X=features,
        pdu_ids=pdu_ids,
        subset=args.subset
    )

    logger.info(f"\n✓ Saved: {out_file}")
    logger.info(f"Shape: {features.shape}")
    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
