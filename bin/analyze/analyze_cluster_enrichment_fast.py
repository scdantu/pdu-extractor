#!/usr/bin/env python3
"""
Fast cluster enrichment analysis using PyCoM metadata reference.

Uses pre-built pdb_metadata_pycom.json (no API calls, instant analysis).

Run after: python3 scripts/build_pdb_metadata_via_pycom.py

Example:
  python3 scripts/analyze_cluster_enrichment_fast.py \
    --clusters analysis/clusters_umap/pdu_clusters_L_umap.csv \
    --db /path/to/pdus_L.sqlite \
    --metadata pdb_metadata_pycom.json \
    --aa L
"""

import sqlite3
import json
from pathlib import Path
from collections import Counter, defaultdict
import argparse

import numpy as np
import pandas as pd


def load_metadata(metadata_file):
    """Load PyCoM metadata reference."""
    with open(metadata_file) as f:
        metadata = json.load(f)
    print(f"Loaded metadata for {len(metadata)} PDB entries via PyCoM")
    return metadata


def query_pdu_to_pdb(db_path, aa):
    """Map PDU IDs to PDB IDs."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(
        "SELECT id, pdb_id FROM pdu WHERE reference_residue_one_letter = ?",
        (aa,)
    )
    pdu_to_pdb = {row[0]: row[1] for row in cursor.fetchall()}
    conn.close()
    return pdu_to_pdb


def analyze_enrichment(clusters_df, pdu_to_pdb, metadata, aa):
    """Analyze cluster enrichment in CATH, EC, GO via PyCoM."""

    # Map PDUs to PDB IDs
    clusters_df['pdb_id'] = clusters_df['pdu_id'].map(pdu_to_pdb)
    clusters_df = clusters_df.dropna(subset=['pdb_id'])

    # Map PDB IDs to PyCoM metadata
    clusters_df['cath'] = clusters_df['pdb_id'].apply(
        lambda x: metadata.get(x.lower(), {}).get('cath')
    )
    clusters_df['enzyme_ec'] = clusters_df['pdb_id'].apply(
        lambda x: metadata.get(x.lower(), {}).get('enzyme_ec')
    )
    clusters_df['biological_process'] = clusters_df['pdb_id'].apply(
        lambda x: metadata.get(x.lower(), {}).get('biological_process')
    )
    clusters_df['molecular_function'] = clusters_df['pdb_id'].apply(
        lambda x: metadata.get(x.lower(), {}).get('molecular_function')
    )
    clusters_df['cellular_component'] = clusters_df['pdb_id'].apply(
        lambda x: metadata.get(x.lower(), {}).get('cellular_component')
    )
    clusters_df['domain'] = clusters_df['pdb_id'].apply(
        lambda x: metadata.get(x.lower(), {}).get('domain')
    )
    clusters_df['ligand'] = clusters_df['pdb_id'].apply(
        lambda x: metadata.get(x.lower(), {}).get('ligand')
    )
    clusters_df['cofactor'] = clusters_df['pdb_id'].apply(
        lambda x: metadata.get(x.lower(), {}).get('cofactor')
    )
    clusters_df['ptm'] = clusters_df['pdb_id'].apply(
        lambda x: metadata.get(x.lower(), {}).get('ptm')
    )
    clusters_df['disease'] = clusters_df['pdb_id'].apply(
        lambda x: metadata.get(x.lower(), {}).get('disease')
    )
    clusters_df['uniprot'] = clusters_df['pdb_id'].apply(
        lambda x: metadata.get(x.lower(), {}).get('uniprot_id')
    )

    # Remove noise
    assigned_df = clusters_df[clusters_df['cluster'] != -1].copy()

    print(f"\n{'='*80}")
    print(f"CLUSTER ENRICHMENT: AA={aa}")
    print(f"{'='*80}")
    print(f"Total PDUs: {len(clusters_df):,}")
    print(f"Assigned (not noise): {len(assigned_df):,}")
    print(f"Clusters: {assigned_df['cluster'].nunique()}")

    # CATH enrichment (from PyCoM)
    cath_assigned = assigned_df.dropna(subset=['cath'])
    if len(cath_assigned) > 0:
        print(f"\n### CATH Superfamily (via PyCoM)")
        print(f"PDUs with CATH data: {len(cath_assigned):,}/{len(assigned_df):,} ({100*len(cath_assigned)/len(assigned_df):.1f}%)")

        # Top CATH families
        cath_counts = cath_assigned['cath'].value_counts()
        print(f"\nTop CATH superfamilies:")
        for sf, count in cath_counts.head(10).items():
            pct = 100 * count / len(cath_assigned)
            print(f"  {sf}: {count:,} ({pct:.1f}%)")

        # Distribution per cluster
        contingency_cath = pd.crosstab(
            cath_assigned['cath'],
            cath_assigned['cluster']
        )
        print(f"\nCATH × Cluster contingency: {contingency_cath.shape[0]} families × {contingency_cath.shape[1]} clusters")

        # Save
        contingency_cath.to_csv(f"contingency_cath_{aa}.csv")
        print(f"  Saved: contingency_cath_{aa}.csv")

    # EC enrichment (from PyCoM)
    ec_assigned = assigned_df.dropna(subset=['enzyme_ec'])
    if len(ec_assigned) > 0:
        print(f"\n### EC Numbers (via PyCoM)")
        print(f"PDUs with EC data: {len(ec_assigned):,}/{len(assigned_df):,} ({100*len(ec_assigned)/len(assigned_df):.1f}%)")

        # Flatten EC lists
        ec_flat = []
        cluster_flat = []
        for idx, row in ec_assigned.iterrows():
            if row['enzyme_ec']:
                ec_entries = row['enzyme_ec'] if isinstance(row['enzyme_ec'], list) else [row['enzyme_ec']]
                for ec in ec_entries:
                    ec_flat.append(ec)
                    cluster_flat.append(row['cluster'])

        if ec_flat:
            ec_df = pd.DataFrame({'ec': ec_flat, 'cluster': cluster_flat})
            ec_counts = ec_df['ec'].value_counts()

            print(f"\nTop EC classes:")
            for ec, count in ec_counts.head(10).items():
                pct = 100 * count / len(ec_df)
                print(f"  {ec}: {count:,} ({pct:.1f}%)")

            # Save
            contingency_ec = pd.crosstab(ec_df['ec'], ec_df['cluster'])
            contingency_ec.to_csv(f"contingency_ec_{aa}.csv")
            print(f"  Saved: contingency_ec_{aa}.csv")

    # GO enrichment (from PyCoM)
    # Biological Process
    bp_assigned = assigned_df.dropna(subset=['biological_process'])
    if len(bp_assigned) > 0:
        print(f"\n### Gene Ontology - Biological Process (via PyCoM)")
        print(f"PDUs with BP data: {len(bp_assigned):,}/{len(assigned_df):,} ({100*len(bp_assigned)/len(assigned_df):.1f}%)")

        bp_flat = []
        cluster_flat = []
        for idx, row in bp_assigned.iterrows():
            if row['biological_process']:
                bp_list = row['biological_process'] if isinstance(row['biological_process'], list) else [row['biological_process']]
                for bp in bp_list:
                    bp_flat.append(bp)
                    cluster_flat.append(row['cluster'])

        if bp_flat:
            bp_df = pd.DataFrame({'bp': bp_flat, 'cluster': cluster_flat})
            bp_counts = bp_df['bp'].value_counts()

            print(f"\nTop biological processes:")
            for bp, count in bp_counts.head(10).items():
                pct = 100 * count / len(bp_df)
                print(f"  {bp}: {count:,} ({pct:.1f}%)")

    # Molecular Function
    mf_assigned = assigned_df.dropna(subset=['molecular_function'])
    if len(mf_assigned) > 0:
        print(f"\n### Gene Ontology - Molecular Function (via PyCoM)")
        print(f"PDUs with MF data: {len(mf_assigned):,}/{len(assigned_df):,} ({100*len(mf_assigned)/len(assigned_df):.1f}%)")

        mf_flat = []
        cluster_flat = []
        for idx, row in mf_assigned.iterrows():
            if row['molecular_function']:
                mf_list = row['molecular_function'] if isinstance(row['molecular_function'], list) else [row['molecular_function']]
                for mf in mf_list:
                    mf_flat.append(mf)
                    cluster_flat.append(row['cluster'])

        if mf_flat:
            mf_df = pd.DataFrame({'mf': mf_flat, 'cluster': cluster_flat})
            mf_counts = mf_df['mf'].value_counts()

            print(f"\nTop molecular functions:")
            for mf, count in mf_counts.head(10).items():
                pct = 100 * count / len(mf_df)
                print(f"  {mf}: {count:,} ({pct:.1f}%)")

    # Domain enrichment (from PyCoM)
    domain_assigned = assigned_df.dropna(subset=['domain'])
    if len(domain_assigned) > 0:
        print(f"\n### Protein Domains (via PyCoM)")
        print(f"PDUs with domain data: {len(domain_assigned):,}/{len(assigned_df):,} ({100*len(domain_assigned)/len(assigned_df):.1f}%)")

        domain_flat = []
        cluster_flat = []
        for idx, row in domain_assigned.iterrows():
            if row['domain']:
                domain_list = row['domain'] if isinstance(row['domain'], list) else [row['domain']]
                for domain in domain_list:
                    domain_flat.append(domain)
                    cluster_flat.append(row['cluster'])

        if domain_flat:
            domain_df = pd.DataFrame({'domain': domain_flat, 'cluster': cluster_flat})
            domain_counts = domain_df['domain'].value_counts()

            print(f"\nTop protein domains:")
            for domain, count in domain_counts.head(10).items():
                pct = 100 * count / len(domain_df)
                print(f"  {domain}: {count:,} ({pct:.1f}%)")

    # Save full enriched data
    clusters_df.to_csv(f"clusters_enriched_{aa}.csv", index=False)
    print(f"\nSaved: clusters_enriched_{aa}.csv")

    return clusters_df


def main():
    parser = argparse.ArgumentParser(
        description="Fast cluster enrichment analysis (uses pre-built metadata reference)"
    )
    parser.add_argument("--clusters", required=True, help="Cluster CSV")
    parser.add_argument("--db", required=True, help="PDU database")
    parser.add_argument("--metadata", required=True, help="PDB metadata reference JSON")
    parser.add_argument("--aa", default="L", help="Amino acid")

    args = parser.parse_args()

    # Load
    print(f"\n[1/3] Loading cluster assignments...")
    clusters_df = pd.read_csv(args.clusters)
    print(f"  ✓ {len(clusters_df):,} PDUs")

    print(f"\n[2/3] Mapping PDU→PDB from database...")
    pdu_to_pdb = query_pdu_to_pdb(args.db, args.aa)
    print(f"  ✓ {len(pdu_to_pdb):,} PDUs mapped")

    print(f"\n[3/3] Loading metadata reference...")
    metadata = load_metadata(args.metadata)

    # Analyze
    analyze_enrichment(clusters_df, pdu_to_pdb, metadata, args.aa)

    print(f"\n{'='*80}")
    print("✓ Analysis complete")
    print(f"{'='*80}\n")


if __name__ == "__main__":
    main()
