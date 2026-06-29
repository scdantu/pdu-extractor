#!/usr/bin/env python3
"""
Analyze CATH and EC classification enrichment in PDU clusters.

Maps PDB IDs → CATH domains → EC numbers → checks if clusters separate by function/structure.

Run from mdpdu directory:
  python3 pdu-extractor/scripts/analyze_cath_ec_enrichment.py \
    --clusters analysis/clusters_umap/pdu_clusters_L_umap.csv \
    --db /path/to/pdu_db/pdus_L.sqlite \
    --aa L \
    --out-dir analysis/cath_ec_enrichment
"""

import argparse
import sqlite3
from pathlib import Path
from collections import Counter, defaultdict
import json
import sys

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False
    print("⚠️  requests not available; using offline CATH cache only")


class PDBToCATH:
    """Map PDB IDs to CATH domains."""

    def __init__(self, cache_file=None):
        self.cache = {}
        self.cache_file = cache_file
        if cache_file and Path(cache_file).exists():
            self._load_cache()

    def _load_cache(self):
        """Load cached CATH mappings."""
        try:
            with open(self.cache_file) as f:
                self.cache = json.load(f)
            print(f"Loaded CATH cache: {len(self.cache)} PDB IDs")
        except Exception as e:
            print(f"Error loading cache: {e}")

    def _save_cache(self):
        """Save CATH mappings to cache."""
        if self.cache_file:
            with open(self.cache_file, 'w') as f:
                json.dump(self.cache, f, indent=2)

    def get_cath_domain(self, pdb_id):
        """
        Get CATH classification for PDB ID.

        Returns:
        - dict with keys: cath_id, sf_id (superfamily), sf_name
        - None if not found
        """
        if pdb_id in self.cache:
            return self.cache[pdb_id]

        # Try API if requests available
        if HAS_REQUESTS:
            try:
                url = f"https://www.cathdb.info/api/rest/pdb_entry/{pdb_id.lower()}"
                resp = requests.get(url, timeout=5)
                if resp.status_code == 200:
                    data = resp.json()
                    if data and 'data' in data and data['data']:
                        entry = data['data'][0]
                        result = {
                            'cath_id': entry.get('cath_id'),
                            'sf_id': entry.get('superfamily_id'),
                            'sf_name': entry.get('superfamily_description')
                        }
                        self.cache[pdb_id] = result
                        return result
            except Exception as e:
                pass

        return None


class PDBToEC:
    """Map PDB IDs to EC numbers."""

    def __init__(self, cache_file=None):
        self.cache = {}
        self.cache_file = cache_file
        if cache_file and Path(cache_file).exists():
            self._load_cache()

    def _load_cache(self):
        """Load cached EC mappings."""
        try:
            with open(self.cache_file) as f:
                self.cache = json.load(f)
            print(f"Loaded EC cache: {len(self.cache)} PDB IDs")
        except Exception as e:
            print(f"Error loading cache: {e}")

    def _save_cache(self):
        """Save EC mappings to cache."""
        if self.cache_file:
            with open(self.cache_file, 'w') as f:
                json.dump(self.cache, f, indent=2)

    def get_ec_number(self, pdb_id):
        """
        Get EC number for PDB ID.

        Returns: list of EC numbers or empty list
        """
        if pdb_id in self.cache:
            return self.cache[pdb_id]

        # Try API if requests available
        if HAS_REQUESTS:
            try:
                url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/enzyme/{pdb_id.lower()}"
                resp = requests.get(url, timeout=5)
                if resp.status_code == 200:
                    data = resp.json()
                    if pdb_id.lower() in data:
                        ec_list = [e['ec_number'] for e in data[pdb_id.lower()]]
                        self.cache[pdb_id] = ec_list
                        return ec_list
            except Exception as e:
                pass

        return []


def query_pdu_pdb_ids(db_path, aa):
    """Extract all PDB IDs for given amino acid."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Get unique PDB IDs
    cursor.execute(
        "SELECT DISTINCT pdb_id FROM pdu WHERE reference_residue_one_letter = ?",
        (aa,)
    )
    pdb_ids = [row[0] for row in cursor.fetchall()]
    conn.close()

    return pdb_ids


def main():
    parser = argparse.ArgumentParser(
        description="Analyze CATH/EC enrichment in PDU clusters"
    )
    parser.add_argument("--clusters", required=True, help="Cluster CSV file")
    parser.add_argument("--db", required=True, help="PDU SQLite database")
    parser.add_argument("--aa", default="L", help="Amino acid")
    parser.add_argument("--out-dir", default="analysis/cath_ec_enrichment", help="Output directory")
    parser.add_argument("--cath-cache", default=None, help="CATH cache file")
    parser.add_argument("--ec-cache", default=None, help="EC cache file")

    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 80)
    print(f"CATH/EC ENRICHMENT ANALYSIS: {args.aa}")
    print("=" * 80)

    # Load clusters
    print(f"\n[1/5] Loading clusters...")
    clusters_df = pd.read_csv(args.clusters)
    print(f"  ✓ Loaded {len(clusters_df)} PDU cluster assignments")

    # Get unique PDB IDs
    print(f"\n[2/5] Extracting PDB IDs from database...")
    pdb_ids = query_pdu_pdb_ids(args.db, args.aa)
    print(f"  ✓ Found {len(pdb_ids)} unique PDB IDs for {args.aa}")

    # Map PDUs to their PDB IDs (from database)
    print(f"\n[3/5] Building PDU→PDB mapping...")
    conn = sqlite3.connect(args.db)
    cursor = conn.cursor()
    cursor.execute(
        "SELECT id, pdb_id FROM pdu WHERE reference_residue_one_letter = ?",
        (args.aa,)
    )
    pdu_to_pdb = {row[0]: row[1] for row in cursor.fetchall()}
    conn.close()
    print(f"  ✓ Mapped {len(pdu_to_pdb)} PDUs to PDB IDs")

    # Add PDB IDs to cluster dataframe
    clusters_df['pdb_id'] = clusters_df['pdu_id'].map(pdu_to_pdb)
    clusters_df = clusters_df.dropna(subset=['pdb_id'])

    # Initialize mappers (offline mode without API keys)
    cath_mapper = PDBToCATH(cache_file=args.cath_cache)
    ec_mapper = PDBToEC(cache_file=args.ec_cache)

    print(f"\n[4/5] Analyzing CATH/EC distribution...")

    # Extract CATH/EC for each unique PDB
    pdb_to_cath = {}
    pdb_to_ec = {}

    unique_pdbs = clusters_df['pdb_id'].unique()
    print(f"  Processing {len(unique_pdbs)} unique PDB IDs...")

    for i, pdb_id in enumerate(unique_pdbs):
        if (i + 1) % 100 == 0:
            print(f"    {i+1}/{len(unique_pdbs)}")

        cath = cath_mapper.get_cath_domain(pdb_id)
        if cath:
            pdb_to_cath[pdb_id] = cath

        ec = ec_mapper.get_ec_number(pdb_id)
        if ec:
            pdb_to_ec[pdb_id] = ec

    cath_mapper._save_cache()
    ec_mapper._save_cache()

    print(f"  ✓ CATH: {len(pdb_to_cath)}/{len(unique_pdbs)} PDB IDs mapped")
    print(f"  ✓ EC:   {len(pdb_to_ec)}/{len(unique_pdbs)} PDB IDs mapped")

    # Map to clusters
    clusters_df['cath_sf'] = clusters_df['pdb_id'].apply(
        lambda x: pdb_to_cath.get(x, {}).get('sf_name') if x in pdb_to_cath else None
    )
    clusters_df['ec'] = clusters_df['pdb_id'].apply(
        lambda x: pdb_to_ec.get(x, [None])[0] if x in pdb_to_ec else None
    )

    # Remove noise points
    assigned_df = clusters_df[clusters_df['cluster'] != -1].copy()

    print(f"\n[5/5] Statistical analysis...")

    # CATH enrichment (if we have CATH data)
    cath_assigned = assigned_df.dropna(subset=['cath_sf'])
    if len(cath_assigned) > 0:
        print(f"\n### CATH Superfamily Enrichment")
        print(f"PDUs with CATH data: {len(cath_assigned)}/{len(assigned_df)}")

        contingency_cath = pd.crosstab(
            cath_assigned['cath_sf'],
            cath_assigned['cluster'],
            margins=False
        )

        if contingency_cath.shape[0] > 1 and contingency_cath.shape[1] > 1:
            print(f"  Contingency table: {contingency_cath.shape[0]} CATH families × {contingency_cath.shape[1]} clusters")
            contingency_cath.to_csv(out_dir / f"contingency_cath_{args.aa}.csv")

        # Show top CATH superfamilies
        cath_counts = cath_assigned['cath_sf'].value_counts()
        print(f"\nTop CATH superfamilies:")
        for sf, count in cath_counts.head(10).items():
            pct = 100 * count / len(cath_assigned)
            print(f"  {sf}: {count:,} ({pct:.1f}%)")

    # EC enrichment (if we have EC data)
    ec_assigned = assigned_df.dropna(subset=['ec'])
    if len(ec_assigned) > 0:
        print(f"\n### EC Number Enrichment")
        print(f"PDUs with EC data: {len(ec_assigned)}/{len(assigned_df)}")

        contingency_ec = pd.crosstab(
            ec_assigned['ec'],
            ec_assigned['cluster'],
            margins=False
        )

        if contingency_ec.shape[0] > 1 and contingency_ec.shape[1] > 1:
            print(f"  Contingency table: {contingency_ec.shape[0]} EC numbers × {contingency_ec.shape[1]} clusters")
            contingency_ec.to_csv(out_dir / f"contingency_ec_{args.aa}.csv")

        # Show top EC classes
        ec_counts = ec_assigned['ec'].value_counts()
        print(f"\nTop EC classes:")
        for ec, count in ec_counts.head(10).items():
            pct = 100 * count / len(ec_assigned)
            print(f"  {ec}: {count:,} ({pct:.1f}%)")

    # Save detailed results
    clusters_df.to_csv(out_dir / f"clusters_with_cath_ec_{args.aa}.csv", index=False)
    print(f"\nResults saved to {out_dir}/")

    print("\n" + "=" * 80)
    print("✓ Analysis complete")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
