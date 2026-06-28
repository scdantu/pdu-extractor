#!/usr/bin/env python3
"""
Build a unified PDB metadata reference file from all 20 amino acid databases.

Extracts unique PDB IDs across all AAs, then maps to:
  - CATH (structural families)
  - EC (enzyme classification)
  - GO (gene ontology)
  - Pfam (protein families)
  - UniProt (general annotation)

Run once, use for all downstream analysis.

Output: pdb_metadata_reference.json (fast lookup)
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
    print("⚠️  requests not available; will use offline mode only")


class PDBMetadataFetcher:
    """Fetch and cache PDB metadata from various sources."""

    def __init__(self, cache_file=None, use_api=True):
        self.cache = {}
        self.cache_file = cache_file
        self.use_api = use_api and HAS_REQUESTS
        self.fetch_count = 0
        self.cache_hits = 0

        if cache_file and Path(cache_file).exists():
            self._load_cache()

    def _load_cache(self):
        """Load existing metadata cache."""
        try:
            with open(self.cache_file) as f:
                self.cache = json.load(f)
            print(f"✓ Loaded cache: {len(self.cache)} PDB entries")
        except Exception as e:
            print(f"⚠️  Error loading cache: {e}")

    def _save_cache(self):
        """Save metadata cache."""
        if self.cache_file:
            with open(self.cache_file, 'w') as f:
                json.dump(self.cache, f, indent=2)
            print(f"✓ Saved cache: {len(self.cache)} entries ({self.fetch_count} new fetches)")

    def fetch_metadata(self, pdb_id):
        """
        Fetch metadata for a PDB ID.

        Returns dict with keys:
          - cath_sf_id: CATH superfamily ID
          - cath_sf_name: CATH superfamily name
          - ec: List of EC numbers
          - go: List of GO terms
          - pfam: List of Pfam IDs
          - uniprot: UniProt ID
        """
        pdb_id = pdb_id.lower()

        if pdb_id in self.cache:
            self.cache_hits += 1
            return self.cache[pdb_id]

        result = {
            'pdb_id': pdb_id,
            'cath_sf_id': None,
            'cath_sf_name': None,
            'ec': [],
            'go': [],
            'pfam': [],
            'uniprot': None,
            'resolution': None,
        }

        if not self.use_api:
            self.cache[pdb_id] = result
            return result

        # Fetch from PDBe API (comprehensive)
        try:
            # PDBe API for CATH, EC, GO, Pfam
            url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/{pdb_id}"
            resp = requests.get(url, timeout=10)

            if resp.status_code == 200:
                data = resp.json()
                if pdb_id in data:
                    entry = data[pdb_id][0] if data[pdb_id] else {}

                    # Resolution
                    if 'resolution' in entry:
                        result['resolution'] = entry['resolution']

                    # Get chain-level annotation
                    if 'molecules' in entry:
                        for mol in entry['molecules']:
                            # Pfam
                            if 'pfam_accession' in mol:
                                result['pfam'].extend(mol['pfam_accession'])

            # CATH API
            url = f"https://www.cathdb.info/api/rest/pdb_entry/{pdb_id}"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                if data and 'data' in data and data['data']:
                    entry = data['data'][0]
                    result['cath_sf_id'] = entry.get('superfamily_id')
                    result['cath_sf_name'] = entry.get('superfamily_description')

            # EC via PDBe
            url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/enzyme/{pdb_id}"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                if pdb_id in data:
                    for enzyme in data[pdb_id]:
                        if 'ec_number' in enzyme:
                            result['ec'].append(enzyme['ec_number'])

            # GO via PDBe
            url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/go/{pdb_id}"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                if pdb_id in data:
                    for go_entry in data[pdb_id]:
                        if 'go_id' in go_entry:
                            result['go'].append({
                                'id': go_entry['go_id'],
                                'name': go_entry.get('name', ''),
                                'category': go_entry.get('category', '')
                            })

            # UniProt via PDBe
            url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/uniprot/{pdb_id}"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                if pdb_id in data:
                    uniprot_entries = data[pdb_id]
                    if uniprot_entries:
                        result['uniprot'] = uniprot_entries[0].get('uniprot_id')

            self.fetch_count += 1

        except Exception as e:
            print(f"⚠️  Error fetching {pdb_id}: {e}")

        # Clean up empty lists
        if not result['ec']:
            result['ec'] = None
        if not result['go']:
            result['go'] = None
        if not result['pfam']:
            result['pfam'] = None

        self.cache[pdb_id] = result
        return result


def extract_unique_pdb_ids(db_dir, aas=None):
    """Extract unique PDB IDs from all amino acid databases."""
    if aas is None:
        aas = [chr(65 + i) for i in range(20)]  # A-T
        aas = [a for a in aas if a not in ['B', 'J', 'O', 'U', 'X', 'Z']]  # Valid AAs

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
        description="Build unified PDB metadata reference"
    )
    parser.add_argument(
        "--db-dir",
        required=True,
        help="Directory containing per-AA sqlite databases"
    )
    parser.add_argument(
        "--out-file",
        default="pdb_metadata_reference.json",
        help="Output file (JSON)"
    )
    parser.add_argument(
        "--api",
        action="store_true",
        default=True,
        help="Use APIs to fetch metadata (default: True)"
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Offline mode (cache only, no API calls)"
    )

    args = parser.parse_args()

    print("\n" + "=" * 80)
    print("BUILD UNIFIED PDB METADATA REFERENCE")
    print("=" * 80)

    # Extract unique PDB IDs from all 20 AAs
    print(f"\n[1/3] Extracting unique PDB IDs from all amino acid databases...")
    pdb_ids, pdb_to_aas = extract_unique_pdb_ids(args.db_dir)
    print(f"\n✓ Found {len(pdb_ids)} unique PDB IDs across 20 amino acids")
    print(f"  (Range: {min(len(v) for v in pdb_to_aas.values())} to {max(len(v) for v in pdb_to_aas.values())} AAs per PDB)")

    # Initialize fetcher
    print(f"\n[2/3] Fetching metadata for each PDB...")
    fetcher = PDBMetadataFetcher(
        cache_file=args.out_file,
        use_api=args.api and not args.offline
    )

    # Fetch metadata
    for i, pdb_id in enumerate(pdb_ids):
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(pdb_ids)} ({100*i/len(pdb_ids):.1f}%)")
            if (i + 1) % 500 == 0:
                fetcher._save_cache()
                print(f"    [checkpoint saved: {fetcher.cache_hits} cache hits, {fetcher.fetch_count} new fetches]")

        metadata = fetcher.fetch_metadata(pdb_id)
        # Add amino acid context
        metadata['amino_acids'] = pdb_to_aas[pdb_id]

    fetcher._save_cache()

    # Save as JSON
    print(f"\n[3/3] Saving metadata to {args.out_file}...")
    with open(args.out_file, 'w') as f:
        json.dump(fetcher.cache, f, indent=2)

    print(f"\n✓ Reference saved: {len(fetcher.cache)} PDB entries")
    print(f"  Cache hits: {fetcher.cache_hits}")
    print(f"  New fetches: {fetcher.fetch_count}")
    print(f"  File size: {Path(args.out_file).stat().st_size / 1e6:.1f} MB")

    print("\n" + "=" * 80)
    print("Next: Use this reference in cluster enrichment analysis")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
