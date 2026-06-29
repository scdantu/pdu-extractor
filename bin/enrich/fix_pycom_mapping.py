#!/usr/bin/env python3
"""Fix PyCoM mapping CSV - unpack PDB ID lists into individual rows."""

import pandas as pd
import ast
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def fix_pycom_mapping(input_csv, output_csv):
    """Expand PDB ID lists into individual rows."""
    logger.info(f"Loading: {input_csv}")
    df = pd.read_csv(input_csv)
    logger.info(f"Input: {len(df):,} rows")

    # Check if pdb_id is string-encoded list
    sample = df['pdb_id'].iloc[0]
    is_list_str = isinstance(sample, str) and sample.startswith('[')

    if not is_list_str:
        logger.warning("PDB IDs already unpacked or in different format")
        return

    logger.info("Unpacking PDB ID lists...")

    new_rows = []
    for _, row in df.iterrows():
        try:
            # Parse the string representation of list
            pdb_ids = ast.literal_eval(row['pdb_id'])

            # Create one row per PDB ID
            for pdb_id in pdb_ids:
                new_row = row.copy()
                new_row['pdb_id'] = pdb_id.upper()  # Normalize to uppercase
                new_rows.append(new_row)

        except Exception as e:
            logger.warning(f"Could not parse row: {e}")

    df_fixed = pd.DataFrame(new_rows)
    logger.info(f"Output: {len(df_fixed):,} rows ({len(df_fixed)/len(df):.1f}x expansion)")

    # Save
    df_fixed.to_csv(output_csv, index=False)
    logger.info(f"✓ Saved: {output_csv}")

    # Summary
    logger.info(f"\nSummary:")
    logger.info(f"  Unique PDB IDs: {df_fixed['pdb_id'].nunique():,}")
    logger.info(f"  With CATH: {df_fixed['cath_class'].notna().sum():,}")
    logger.info(f"  With EC: {df_fixed['enzyme_commission'].notna().sum():,}")
    logger.info(f"  With PTM: {df_fixed['ptm'].notna().sum():,}")
    logger.info(f"  With disease: {df_fixed['disease_name'].notna().sum():,}")


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        input_file = "pycom/pycom_pdb_mapping.csv"
        output_file = "pycom/pycom_pdb_mapping_fixed.csv"
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.replace('.csv', '_fixed.csv')

    fix_pycom_mapping(input_file, output_file)
