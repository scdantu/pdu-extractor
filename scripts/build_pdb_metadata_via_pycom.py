#!/usr/bin/env python3
"""
Build PDB metadata using PyCoM (Brunel's local resource) instead of external APIs.

Much faster and more reliable than hitting PDBe/CATH/EC APIs.

Workflow:
  1. Extract unique PDB IDs from all 20 AA databases
  2. Map PDB → UniProt (cached locally, one-time fetch)
  3. Query PyCoM API with UniProt IDs → get CATH, EC, GO, Pfam, etc.
  4. Save to JSON reference file

Run once, reuse for all downstream analysis.
"""

import sqlite3
import json
from pathlib import Path
from collections import defaultdict
import argparse
import time

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False
    print("⚠️  requests not available")


PYCOM_API_BASE = "https://pycom.brunel.ac.uk/api"


class PyComMetadataFetcher:
    """Fetch protein metadata from PyCoM API using UniProt IDs."""

    def __init__(self, cache_file=None):
        self.cache = {}
        self.cache_file = cache_file
        self.fetch_count = 0
        self.cache_hits = 0

        if cache_file and Path(cache_file).exists():
            self._load_cache()

    def _load_cache(self):
        """Load existing metadata cache."""
        try:
            with open(self.cache_file) as f:
                self.cache = json.load(f)
            print(f"✓ Loaded cache: {len(self.cache)} entries")
        except Exception as e:
            print(f"⚠️  Error loading cache: {e}")

    def _save_cache(self):
        """Save metadata cache."""
        if self.cache_file:
            with open(self.cache_file, 'w') as f:
                json.dump(self.cache, f, indent=2)

    def map_pdb_to_uniprot(self, pdb_id, cache_file=None):
        """
        Map PDB ID to UniProt ID using PDBe API.

        Cache separately to avoid refetching.
        """
        pdb_id = pdb_id.lower()

        if cache_file and Path(cache_file).exists():
            with open(cache_file) as f:
                pdb_uniprot_map = json.load(f)
                if pdb_id in pdb_uniprot_map:
                    return pdb_uniprot_map[pdb_id]
        else:
            pdb_uniprot_map = {}

        try:
            url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/uniprot/{pdb_id}"
            resp = requests.get(url, timeout=10)

            if resp.status_code == 200:
                data = resp.json()
                if pdb_id in data and data[pdb_id]:
                    uniprot_id = data[pdb_id][0].get('uniprot_id')
                    if uniprot_id:
                        pdb_uniprot_map[pdb_id] = uniprot_id

                        # Save mapping
                        if cache_file:
                            with open(cache_file, 'w') as f:
                                json.dump(pdb_uniprot_map, f, indent=2)

                        return uniprot_id
        except Exception as e:
            pass

        return None

    def fetch_from_pycom(self, uniprot_id):
        """
        Query PyCoM API for a UniProt ID.

        Returns metadata with CATH, EC, GO, Pfam, ligand, cofactor, PTM, disease.
        """
        if uniprot_id in self.cache:
            self.cache_hits += 1
            return self.cache[uniprot_id]

        result = {
            'uniprot_id': uniprot_id,
            'cath': None,
            'enzyme_ec': None,
            'biological_process': [],
            'molecular_function': [],
            'cellular_component': [],
            'domain': [],
            'ligand': [],
            'cofactor': [],
            'ptm': [],
            'disease': [],
            'organism': None,
            'length': None,
        }

        try:
            # Query PyCoM API with uniprot_id
            url = f"{PYCOM_API_BASE}/find"
            params = {
                'uniprot_id': uniprot_id,
                'per_page': 1  # Should return 0 or 1 result
            }

            resp = requests.get(url, params=params, timeout=10)

            if resp.status_code == 200:
                data = resp.json()

                if 'results' in data and data['results']:
                    protein = data['results'][0]

                    # Extract fields from PyCoM response
                    result['length'] = protein.get('length')
                    result['organism'] = protein.get('organism')
                    result['cath'] = protein.get('cath')
                    result['enzyme_ec'] = protein.get('enzyme')

                    # GO terms
                    if 'biological_process' in protein and protein['biological_process']:
                        result['biological_process'] = protein['biological_process']

                    if 'molecular_function' in protein and protein['molecular_function']:
                        result['molecular_function'] = protein['molecular_function']

                    if 'cellular_component' in protein and protein['cellular_component']:
                        result['cellular_component'] = protein['cellular_component']

                    # Pfam domains
                    if 'domain' in protein and protein['domain']:
                        result['domain'] = protein['domain']

                    # Ligands, cofactors, PTMs, diseases
                    if 'ligand' in protein and protein['ligand']:
                        result['ligand'] = protein['ligand']

                    if 'cofactor' in protein and protein['cofactor']:
                        result['cofactor'] = protein['cofactor']

                    if 'ptm' in protein and protein['ptm']:
                        result['ptm'] = protein['ptm']

                    if 'disease' in protein and protein['disease']:
                        result['disease'] = protein['disease']

                self.fetch_count += 1

        except Exception as e:
            print(f"⚠️  Error fetching {uniprot_id}: {e}")

        self.cache[uniprot_id] = result
        return result


def extract_unique_pdb_ids(db_dir, aas=None):
    """Extract unique PDB IDs from all amino acid databases."""
    if aas is None:
        aas = [chr(65 + i) for i in range(20)]  # A-T
        aas = [a for a in aas if a not in ['B', 'J', 'O', 'U', 'X', 'Z']]

    pdb_ids = set()
    pdb_to_aas = defaultdict(list)

    for aa in aas:
        db_path = Path(db_dir) / f"pdus_{aa}.sqlite"

        if not db_path.exists():
            print(f"⚠️  {db_path} not found, skipping")
            continue

        try:
            conn = sqlite3.connect(str(db_path))
            cursor = conn.cursor()
            cursor.execute("SELECT DISTINCT pdb_id FROM pdu WHERE pdb_id IS NOT NULL")

            for row in cursor.fetchall():
                pdb_id = row[0]
                if pdb_id:
                    pdb_ids.add(pdb_id)
                    pdb_to_aas[pdb_id].append(aa)

            conn.close()
            print(f"✓ AA={aa}: extracted PDB IDs")
        except Exception as e:
            print(f"⚠️  Error reading {db_path}: {e}")

    return sorted(pdb_ids), pdb_to_aas


def main():
    parser = argparse.ArgumentParser(
        description="Build PDB metadata using PyCoM API"
    )
    parser.add_argument(
        "--db-dir",
        required=True,
        help="Directory containing per-AA sqlite databases"
    )
    parser.add_argument(
        "--out-file",
        default="pdb_metadata_pycom.json",
        help="Output metadata file (JSON)"
    )
    parser.add_argument(
        "--pdb-uniprot-cache",
        default=".pdb_uniprot_map.json",
        help="Cache file for PDB→UniProt mappings"
    )

    args = parser.parse_args()

    print("\n" + "=" * 80)
    print("BUILD PDB METADATA USING PyCoM API")
    print("=" * 80)

    # Step 1: Extract unique PDB IDs
    print(f"\n[1/3] Extracting unique PDB IDs from all amino acid databases...")
    pdb_ids, pdb_to_aas = extract_unique_pdb_ids(args.db_dir)
    print(f"\n✓ Found {len(pdb_ids)} unique PDB IDs across 20 amino acids")

    # Step 2: Map PDB → UniProt
    print(f"\n[2/3] Mapping PDB IDs to UniProt IDs...")
    fetcher = PyComMetadataFetcher(cache_file=args.out_file)

    pdb_uniprot_cache = {}
    uniprot_ids = []

    for i, pdb_id in enumerate(pdb_ids):
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(pdb_ids)}")

        uniprot_id = fetcher.map_pdb_to_uniprot(pdb_id, cache_file=args.pdb_uniprot_cache)

        if uniprot_id:
            pdb_uniprot_cache[pdb_id] = uniprot_id
            uniprot_ids.append(uniprot_id)

    print(f"✓ Mapped {len(uniprot_ids)} PDB IDs to UniProt")

    # Step 3: Query PyCoM for each UniProt ID
    print(f"\n[3/3] Querying PyCoM API for metadata...")

    fetcher = PyComMetadataFetcher(cache_file=args.out_file)

    metadata_by_uniprot = {}

    for i, uniprot_id in enumerate(uniprot_ids):
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(uniprot_ids)} ({100*i/len(uniprot_ids):.1f}%)")
            if (i + 1) % 500 == 0:
                fetcher._save_cache()
                print(f"    [checkpoint: {fetcher.cache_hits} cache hits, {fetcher.fetch_count} new fetches]")

        metadata = fetcher.fetch_from_pycom(uniprot_id)
        metadata_by_uniprot[uniprot_id] = metadata

    fetcher._save_cache()

    # Build final output: PDB ID → metadata (via UniProt)
    final_metadata = {}
    for pdb_id, uniprot_id in pdb_uniprot_cache.items():
        if uniprot_id in metadata_by_uniprot:
            final_metadata[pdb_id] = metadata_by_uniprot[uniprot_id]
            final_metadata[pdb_id]['pdb_id'] = pdb_id
            final_metadata[pdb_id]['amino_acids'] = pdb_to_aas[pdb_id]

    # Save
    with open(args.out_file, 'w') as f:
        json.dump(final_metadata, f, indent=2)

    print(f"\n✓ Saved metadata for {len(final_metadata)} PDB entries")
    print(f"  File: {args.out_file}")
    print(f"  Size: {Path(args.out_file).stat().st_size / 1e6:.1f} MB")
    print(f"  PyCoM queries: {fetcher.fetch_count} new, {fetcher.cache_hits} cached")

    print("\n" + "=" * 80)
    print("Next: Use this reference in cluster enrichment analysis")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
