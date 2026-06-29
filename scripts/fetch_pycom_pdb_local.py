#!/usr/bin/env python3
"""
Build PDB → CATH/EC/GO mapping using LOCAL PyCoM database.

Downloads the PyCoM database locally, queries for all entries with PDB structures,
adds biological features (CATH, EC, GO, Pfam), and creates a reference CSV.

Setup:
  1. Download PyCoM database from https://pycom.brunel.ac.uk/database.html
  2. Set DB_PATH and MAT_PATH below

Usage:
  python3 fetch_pycom_pdb_local.py \
    --db-path /path/to/pycom.db \
    --mat-path /path/to/pycom.mat \
    --out-file pdb_pycom_mapping.csv
"""

import argparse
from pathlib import Path
import logging

import pandas as pd

try:
    from pycom import PyCom, ProteinParams
    HAS_PYCOM = True
except ImportError:
    HAS_PYCOM = False
    print("Error: PyCoM not installed")
    print("Install with: pip install pycom")
    exit(1)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Build PDB → CATH/EC/GO mapping from local PyCoM")
    parser.add_argument("--db-path", required=True, help="Path to pycom.db")
    parser.add_argument("--mat-path", required=True, help="Path to pycom.mat")
    parser.add_argument("--out-file", default="pdb_pycom_mapping.csv", help="Output CSV file")

    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info("Build PDB → CATH/EC/GO mapping from PyCoM database")
    logger.info("=" * 80)

    # Verify database files exist
    if not Path(args.db_path).exists():
        logger.error(f"Database file not found: {args.db_path}")
        exit(1)

    if not Path(args.mat_path).exists():
        logger.error(f"Matrix file not found: {args.mat_path}")
        exit(1)

    logger.info(f"DB: {args.db_path}")
    logger.info(f"Matrix: {args.mat_path}")

    # Initialize PyCoM
    logger.info("\n[1/3] Initializing PyCoM...")
    obj_pycom = PyCom(db_path=args.db_path, mat_path=args.mat_path)
    logger.info("  ✓ PyCoM initialized")

    # Query for all entries with PDB structures
    logger.info("\n[2/3] Querying for entries with PDB structures...")
    query_params = {ProteinParams.HAS_PDB: True}

    df = obj_pycom.find(query_params)
    logger.info(f"  ✓ Found {len(df):,} entries with PDB structures")

    # Add biological features
    logger.info("\n[3/3] Adding biological features...")
    obj_loader = obj_pycom.get_data_loader()

    logger.info("  - Adding CATH...")
    df = obj_loader.add_cath_class(df, force_single_entry=False)

    logger.info("  - Adding EC...")
    df = obj_loader.add_enzyme_commission(df, force_single_entry=False)

    logger.info("  - Adding PDB IDs...")
    df = obj_loader.add_pdbs(df, force_single_entry=False)

    logger.info("  - Adding Pfam...")
    df = obj_loader.add_domain(df, force_single_entry=False)

    logger.info("  - Adding diseases...")
    df = obj_loader.add_diseases(df, force_single_entry=False)

    # Select relevant columns
    columns_to_keep = [
        'uniprot_id',
        'pdb_id',
        'cath_class',
        'enzyme_commission',
        'domain',
        'disease_name',
        'organism_id',
        'sequence_length',
        'helix_frac',
        'turn_frac',
        'strand_frac',
    ]

    # Check which columns exist
    available_cols = [col for col in columns_to_keep if col in df.columns]
    df_filtered = df[available_cols].copy()

    # Save to CSV
    df_filtered.to_csv(args.out_file, index=False)

    logger.info(f"\n✓ Saved: {args.out_file}")

    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total entries: {len(df_filtered):,}")
    logger.info(f"With CATH: {df_filtered['cath_class'].notna().sum():,}")
    logger.info(f"With EC: {df_filtered['enzyme_commission'].notna().sum():,}")
    logger.info(f"With Pfam: {df_filtered['domain'].notna().sum():,}")
    logger.info(f"With disease: {df_filtered['disease_name'].notna().sum():,}")

    logger.info(f"\nColumns in output:")
    for col in df_filtered.columns:
        logger.info(f"  - {col}")

    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
