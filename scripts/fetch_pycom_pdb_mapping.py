#!/usr/bin/env python3
"""
Fetch PyCoM entries with PDB IDs and build local reference mapping.

Query PyCoM API for all proteins with PDB entries, extract:
  - PDB ID
  - CATH classification
  - EC number
  - GO terms (biological process, molecular function, cellular component)
  - Pfam domains
  - Disease associations

Build local CSV that maps PDB ID → metadata.

Usage:
  python3 fetch_pycom_pdb_mapping.py --out-file pdb_pycom_mapping.csv
"""

import argparse
from pathlib import Path
import logging
import time

import requests
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

PYCOM_API = "https://pycom.brunel.ac.uk/api"


def query_pycom_with_pdb(batch_size=100, max_pages=None):
    """
    Query PyCoM for entries with PDB structures.

    Yields batch of results per page.
    """
    page = 1
    total_fetched = 0

    while True:
        if max_pages and page > max_pages:
            logger.info(f"Reached max_pages={max_pages}")
            break

        logger.info(f"Fetching page {page}...")

        try:
            params = {
                'has_pdb': 'true',
                'per_page': batch_size,
                'page': page,
            }

            resp = requests.get(f"{PYCOM_API}/find", params=params, timeout=30)
            resp.raise_for_status()

            data = resp.json()

            if 'results' not in data or not data['results']:
                logger.info(f"No more results at page {page}")
                break

            results = data['results']
            logger.info(f"  Page {page}: {len(results)} entries")

            yield results

            total_fetched += len(results)
            page += 1

            # Be nice to the API
            time.sleep(0.5)

        except Exception as e:
            logger.error(f"Error fetching page {page}: {e}")
            break

    logger.info(f"Total entries fetched: {total_fetched}")


def extract_pdb_entries(entry):
    """Extract PDB-specific metadata from PyCoM entry."""

    # PyCoM stores PDB info differently - check structure
    pdb_id = entry.get('pdb')

    if not pdb_id:
        return None

    result = {
        'pdb_id': pdb_id,
        'uniprot_id': entry.get('uniprot_id'),
        'cath': entry.get('cath'),
        'enzyme_ec': entry.get('enzyme'),
        'pfam': ','.join(entry.get('domain', [])) if entry.get('domain') else None,
        'biological_process': ','.join(entry.get('biological_process', [])) if entry.get('biological_process') else None,
        'molecular_function': ','.join(entry.get('molecular_function', [])) if entry.get('molecular_function') else None,
        'cellular_component': ','.join(entry.get('cellular_component', [])) if entry.get('cellular_component') else None,
        'disease': ','.join(entry.get('disease', [])) if entry.get('disease') else None,
        'organism': entry.get('organism'),
        'length': entry.get('length'),
    }

    return result


def main():
    parser = argparse.ArgumentParser(description="Fetch PyCoM entries with PDB IDs")
    parser.add_argument("--out-file", default="pdb_pycom_mapping.csv", help="Output CSV file")
    parser.add_argument("--max-pages", type=int, default=None, help="Max pages to fetch (for testing)")
    parser.add_argument("--batch-size", type=int, default=100, help="Entries per page")

    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info("Fetch PyCoM entries with PDB IDs")
    logger.info("=" * 80)
    logger.info(f"Output: {args.out_file}")

    # Fetch from PyCoM
    all_entries = []

    for page_results in query_pycom_with_pdb(args.batch_size, args.max_pages):
        for entry in page_results:
            pdb_entry = extract_pdb_entries(entry)
            if pdb_entry:
                all_entries.append(pdb_entry)

    logger.info(f"\n✓ Extracted {len(all_entries)} PDB entries")

    if not all_entries:
        logger.error("No entries found! Check PyCoM API or network connection.")
        return

    # Save to CSV
    df = pd.DataFrame(all_entries)
    df.to_csv(args.out_file, index=False)

    logger.info(f"Saved: {args.out_file}")

    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total PDB entries: {len(df)}")
    logger.info(f"Unique CATH classes: {df['cath'].notna().sum()}")
    logger.info(f"Unique EC numbers: {df['enzyme_ec'].notna().sum()}")
    logger.info(f"With GO terms: {df['biological_process'].notna().sum()}")
    logger.info(f"With disease info: {df['disease'].notna().sum()}")

    logger.info(f"\nColumns in output:")
    for col in df.columns:
        logger.info(f"  - {col}")

    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
