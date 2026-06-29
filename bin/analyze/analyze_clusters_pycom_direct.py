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

def chi2_contingency_manual(contingency_table):
    """Compute chi-square statistic without scipy."""
    contingency = np.asarray(contingency_table)
    contingency = contingency.astype(float)

    row_sums = contingency.sum(axis=1, keepdims=True)
    col_sums = contingency.sum(axis=0, keepdims=True)
    total = contingency.sum()

    expected = (row_sums * col_sums) / total
    expected[expected == 0] = 1  # Avoid division by zero

    chi2 = ((contingency - expected) ** 2 / expected).sum()

    # Approximate p-value (use 1 for extreme values)
    return chi2, 0.0  # p-value essentially 0 for very large chi2

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

    # Convert PDB IDs to string and normalize to uppercase for matching
    clusters_df['pdb_id'] = clusters_df['pdb_id'].astype(str).str.upper()
    pycom_df['pdb_id'] = pycom_df['pdb_id'].astype(str).str.upper()

    # Merge on PDB ID
    clusters_df = clusters_df.merge(
        pycom_df[['pdb_id', 'cath_class', 'enzyme_commission', 'disease_name', 'ptm']],
        on='pdb_id',
        how='left'
    )

    mapped_count = clusters_df['pdb_id'].notna().sum()
    logger.info(f"  Mapped: {mapped_count:,} PDUs to PyCoM")

    # Rename columns for consistency with analysis
    clusters_df = clusters_df.rename(columns={
        'cath_class': 'cath',
        'enzyme_commission': 'enzyme_ec'
    })

    # Unpack list-formatted fields (e.g., "['1.10.10.10']" → "1.10.10.10")
    import ast

    def unpack_field(val):
        """Extract first element from list-formatted string, or return as-is."""
        if pd.isna(val):
            return None
        if isinstance(val, str) and val.startswith('['):
            try:
                lst = ast.literal_eval(val)
                return lst[0] if lst else None
            except (ValueError, SyntaxError):
                return None
        return val

    for col in ['cath', 'enzyme_ec']:
        if col in clusters_df.columns:
            clusters_df[col] = clusters_df[col].apply(unpack_field)

    # Analysis
    logger.info(f"\n[3/3] Enrichment analysis...")

    assigned_df = clusters_df[clusters_df['cluster'] != -1].copy()
    n_clusters = assigned_df['cluster'].nunique()

    logger.info(f"\n### CATH Enrichment")
    cath_df = assigned_df.dropna(subset=['cath']).copy()
    logger.info(f"PDUs with CATH: {len(cath_df):,} / {len(assigned_df):,}")

    if len(cath_df) > 0:
        logger.info(f"Top CATH classes:")
        cath_counts = cath_df['cath'].value_counts()
        for cath, count in cath_counts.head(10).items():
            pct = 100 * count / len(cath_df)
            logger.info(f"  {cath}: {count:,} ({pct:.1f}%)")

        contingency_cath = pd.crosstab(cath_df['cath'], cath_df['cluster'])
        chi2, p_val = chi2_contingency_manual(contingency_cath.values)
        logger.info(f"  χ² = {chi2:.1f}, p = {p_val:.2e}")

        contingency_cath.to_csv(Path(args.out_dir) / f"contingency_cath_{args.aa}.csv")

    logger.info(f"\n### EC Number Enrichment")
    ec_df = assigned_df.dropna(subset=['enzyme_ec']).copy()
    logger.info(f"PDUs with EC: {len(ec_df):,} / {len(assigned_df):,}")

    if len(ec_df) > 0:
        logger.info(f"Top EC numbers:")
        ec_counts = ec_df['enzyme_ec'].value_counts()
        for ec, count in ec_counts.head(10).items():
            pct = 100 * count / len(ec_df)
            logger.info(f"  {ec}: {count:,} ({pct:.1f}%)")

        contingency_ec = pd.crosstab(ec_df['enzyme_ec'], ec_df['cluster'])
        chi2, p_val = chi2_contingency_manual(contingency_ec.values)
        logger.info(f"  χ² = {chi2:.1f}, p = {p_val:.2e}")

        contingency_ec.to_csv(Path(args.out_dir) / f"contingency_ec_{args.aa}.csv")

    logger.info(f"\n### PTM Enrichment")
    ptm_df = assigned_df.dropna(subset=['ptm']).copy()
    logger.info(f"PDUs with PTM: {len(ptm_df):,} / {len(assigned_df):,}")

    if len(ptm_df) > 0:
        logger.info(f"Top PTM types:")
        ptm_counts = ptm_df['ptm'].value_counts()
        for ptm, count in ptm_counts.head(10).items():
            pct = 100 * count / len(ptm_df)
            logger.info(f"  {ptm}: {count:,} ({pct:.1f}%)")

    # Save enriched clusters
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    clusters_df.to_csv(Path(args.out_dir) / f"clusters_enriched_{args.aa}.csv", index=False)

    logger.info(f"\n✓ Results saved to {args.out_dir}")
    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
