#!/usr/bin/env python3
"""
Cluster enrichment analysis using PyCoM PDB mapping (direct, no UniProt).

Input:
  1. Cluster CSV (pdu_id, cluster, cluster_probability)
  2. Database (to map PDU ID → PDB ID)
  3. PyCoM mapping CSV (PDB ID → CATH/EC/GO)

Output:
  - Contingency tables: Cluster × CATH, Cluster × EC, etc.
  - Enrichment statistics

Usage:
  python3 analyze_clusters_pycom_direct.py \
    --clusters analysis/clusters_umap/pdu_clusters_L_umap.csv \
    --db pdus_L.sqlite \
    --pycom-mapping pdb_pycom_mapping.csv \
    --aa L \
    --out-dir analysis/pycom_enrichment
"""

import sqlite3
import argparse
from pathlib import Path
from collections import defaultdict
import logging

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def query_pdu_to_pdb(db_path, aa):
    """Map PDU IDs to PDB IDs from database."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute(
        "SELECT id, pdb_id FROM pdu WHERE reference_residue_one_letter = ?",
        (aa,)
    )

    pdu_to_pdb = {row[0]: row[1] for row in cursor.fetchall()}
    conn.close()

    return pdu_to_pdb


def main():
    parser = argparse.ArgumentParser(description="Cluster enrichment via PyCoM PDB mapping")
    parser.add_argument("--clusters", required=True, help="Cluster CSV file")
    parser.add_argument("--db", required=True, help="PDU SQLite database")
    parser.add_argument("--pycom-mapping", required=True, help="PyCoM PDB mapping CSV")
    parser.add_argument("--aa", required=True, help="Amino acid")
    parser.add_argument("--out-dir", default="analysis/pycom_enrichment", help="Output directory")

    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info(f"PyCoM Enrichment Analysis: {args.aa}")
    logger.info("=" * 80)

    # Load data
    logger.info(f"\n[1/3] Loading data...")

    clusters_df = pd.read_csv(args.clusters)
    logger.info(f"  Clusters: {len(clusters_df):,} PDUs")

    pycom_df = pd.read_csv(args.pycom_mapping)
    logger.info(f"  PyCoM mapping: {len(pycom_df):,} PDB entries")

    pdu_to_pdb = query_pdu_to_pdb(args.db, args.aa)
    logger.info(f"  Database mapping: {len(pdu_to_pdb):,} PDUs → PDB")

    # Map PDU → PDB → PyCoM
    logger.info(f"\n[2/3] Mapping PDU → PDB → PyCoM...")

    clusters_df['pdb_id'] = clusters_df['pdu_id'].map(pdu_to_pdb)
    clusters_df = clusters_df.dropna(subset=['pdb_id'])
    logger.info(f"  Mapped: {len(clusters_df):,} / {len(clusters_df) + (len(pdu_to_pdb) - len(clusters_df)):,}")

    # Map to PyCoM metadata
    pycom_map = {row['pdb_id']: row for _, row in pycom_df.iterrows()}

    clusters_df['cath'] = clusters_df['pdb_id'].apply(lambda x: pycom_map.get(str(x).upper(), {}).get('cath'))
    clusters_df['enzyme_ec'] = clusters_df['pdb_id'].apply(lambda x: pycom_map.get(str(x).upper(), {}).get('enzyme_ec'))
    clusters_df['biological_process'] = clusters_df['pdb_id'].apply(lambda x: pycom_map.get(str(x).upper(), {}).get('biological_process'))
    clusters_df['molecular_function'] = clusters_df['pdb_id'].apply(lambda x: pycom_map.get(str(x).upper(), {}).get('molecular_function'))
    clusters_df['cellular_component'] = clusters_df['pdb_id'].apply(lambda x: pycom_map.get(str(x).upper(), {}).get('cellular_component'))

    # Analysis
    logger.info(f"\n[3/3] Enrichment analysis...")

    assigned_df = clusters_df[clusters_df['cluster'] != -1].copy()
    n_clusters = assigned_df['cluster'].nunique()

    logger.info(f"\n### CATH Enrichment")
    cath_df = assigned_df.dropna(subset=['cath'])
    logger.info(f"PDUs with CATH: {len(cath_df):,} / {len(assigned_df):,}")

    if len(cath_df) > 0:
        cath_counts = cath_df['cath'].value_counts()
        print(f"\nTop CATH classes:")
        for cath, count in cath_counts.head(10).items():
            pct = 100 * count / len(cath_df)
            print(f"  {cath}: {count:,} ({pct:.1f}%)")

        contingency_cath = pd.crosstab(cath_df['cath'], cath_df['cluster'])
        chi2, p_val, _, _ = chi2_contingency(contingency_cath.values)
        logger.info(f"  χ² = {chi2:.1f}, p = {p_val:.2e}")

        contingency_cath.to_csv(Path(args.out_dir) / f"contingency_cath_{args.aa}.csv")

    logger.info(f"\n### EC Number Enrichment")
    ec_df = assigned_df.dropna(subset=['enzyme_ec'])
    logger.info(f"PDUs with EC: {len(ec_df):,} / {len(assigned_df):,}")

    if len(ec_df) > 0:
        ec_counts = ec_df['enzyme_ec'].value_counts()
        print(f"\nTop EC classes:")
        for ec, count in ec_counts.head(10).items():
            pct = 100 * count / len(ec_df)
            print(f"  {ec}: {count:,} ({pct:.1f}%)")

        contingency_ec = pd.crosstab(ec_df['enzyme_ec'], ec_df['cluster'])
        chi2, p_val, _, _ = chi2_contingency(contingency_ec.values)
        logger.info(f"  χ² = {chi2:.1f}, p = {p_val:.2e}")

        contingency_ec.to_csv(Path(args.out_dir) / f"contingency_ec_{args.aa}.csv")

    logger.info(f"\n### GO Enrichment")
    bp_df = assigned_df.dropna(subset=['biological_process'])
    logger.info(f"PDUs with GO terms: {len(bp_df):,} / {len(assigned_df):,}")

    if len(bp_df) > 0:
        bp_terms = bp_df['biological_process'].value_counts()
        print(f"\nTop biological processes:")
        for term, count in bp_terms.head(10).items():
            pct = 100 * count / len(bp_df)
            print(f"  {term}: {count:,} ({pct:.1f}%)")

    # Save enriched clusters
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    clusters_df.to_csv(Path(args.out_dir) / f"clusters_enriched_{args.aa}.csv", index=False)

    logger.info(f"\n✓ Results saved to {args.out_dir}")
    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
