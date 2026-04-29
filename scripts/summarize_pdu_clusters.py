#!/usr/bin/env python3
import argparse
import csv
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Summarize PDU clusters by joining cluster labels to SQLite PDU data.")
    parser.add_argument("--db", default="pdu_output/pdus.sqlite", help="PDU SQLite database.")
    parser.add_argument("--clusters-dir", default="analysis/clusters_umap", help="Directory with pdu_clusters_*_*.csv.")
    parser.add_argument("--out-dir", default="analysis/cluster_summaries", help="Output directory.")
    parser.add_argument("--aa", default=None, help="Optional single reference amino acid.")
    parser.add_argument("--space", choices=("latent", "umap"), default="umap", help="Cluster label source.")
    parser.add_argument("--top-n", type=int, default=8, help="Number of top residues/SS labels/examples to report.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cluster_paths = sorted(Path(args.clusters_dir).glob(f"pdu_clusters_*_{args.space}.csv"))
    if args.aa:
        cluster_paths = [Path(args.clusters_dir) / f"pdu_clusters_{args.aa}_{args.space}.csv"]
    if not cluster_paths:
        raise SystemExit(f"No cluster files found in {args.clusters_dir}")

    conn = sqlite3.connect(args.db)
    all_summary_rows = []
    for cluster_path in cluster_paths:
        if not cluster_path.exists():
            print(f"Skipping missing cluster file: {cluster_path}")
            continue
        aa = cluster_path.name.replace("pdu_clusters_", "").replace(f"_{args.space}.csv", "")
        clusters = pd.read_csv(cluster_path)
        clusters = clusters[clusters["cluster"] != -1].copy()
        if clusters.empty:
            print(f"{aa}: no non-noise clusters")
            continue

        pdu_to_cluster = dict(zip(clusters["pdu_id"].astype(int), clusters["cluster"].astype(int)))
        pdu_ids = sorted(pdu_to_cluster)
        cluster_stats = build_cluster_stats(conn, pdu_ids, pdu_to_cluster)

        summary_path = out_dir / f"cluster_summary_{aa}_{args.space}.csv"
        examples_path = out_dir / f"cluster_examples_{aa}_{args.space}.csv"
        write_summary(summary_path, aa, args.space, cluster_stats, args.top_n)
        write_examples(examples_path, aa, args.space, cluster_stats, args.top_n)

        for cluster_id, stats in sorted(cluster_stats.items()):
            all_summary_rows.append(summary_row(aa, args.space, cluster_id, stats, args.top_n))
        print(f"{aa}: wrote {len(cluster_stats)} cluster summaries")

    write_combined_summary(out_dir / f"cluster_summary_all_{args.space}.csv", all_summary_rows)


def build_cluster_stats(conn, pdu_ids, pdu_to_cluster):
    cluster_stats = defaultdict(new_cluster_stats)
    placeholders = ",".join("?" for _ in pdu_ids)

    pdu_rows = conn.execute(
        f"""
        SELECT id, pdb_id, reference_chain_id, reference_residue_number,
               reference_residue_name, reference_residue_one_letter
        FROM pdu
        WHERE id IN ({placeholders})
        """,
        pdu_ids,
    )
    for pdu_id, pdb_id, chain_id, residue_number, residue_name, residue_one_letter in pdu_rows:
        cluster_id = pdu_to_cluster[pdu_id]
        stats = cluster_stats[cluster_id]
        stats["pdu_ids"].append(pdu_id)
        stats["examples"].append(
            {
                "pdu_id": pdu_id,
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "residue_number": residue_number,
                "residue_name": residue_name,
                "residue_one_letter": residue_one_letter,
            }
        )

    residue_rows = conn.execute(
        f"""
        SELECT pdu_id, residue_one_letter, residue_name, secondary_structure, distance_angstrom
        FROM pdu_residue
        WHERE pdu_id IN ({placeholders})
        """,
        pdu_ids,
    )
    for pdu_id, residue_one_letter, residue_name, secondary_structure, distance in residue_rows:
        cluster_id = pdu_to_cluster[pdu_id]
        stats = cluster_stats[cluster_id]
        stats["neighbor_count"] += 1
        stats["neighbor_residues"][residue_one_letter] += 1
        stats["neighbor_residue_names"][residue_name] += 1
        stats["secondary_structure"][secondary_structure or "C"] += 1
        stats["distance_sum"] += float(distance)

    return cluster_stats


def new_cluster_stats():
    return {
        "pdu_ids": [],
        "examples": [],
        "neighbor_count": 0,
        "neighbor_residues": Counter(),
        "neighbor_residue_names": Counter(),
        "secondary_structure": Counter(),
        "distance_sum": 0.0,
    }


def write_summary(path, aa, space, cluster_stats, top_n):
    rows = [summary_row(aa, space, cluster_id, stats, top_n) for cluster_id, stats in sorted(cluster_stats.items())]
    write_dict_rows(path, rows)


def summary_row(aa, space, cluster_id, stats, top_n):
    pdu_count = len(stats["pdu_ids"])
    neighbor_count = stats["neighbor_count"]
    return {
        "aa": aa,
        "space": space,
        "cluster": cluster_id,
        "pdu_count": pdu_count,
        "mean_neighbors": neighbor_count / pdu_count if pdu_count else 0.0,
        "mean_distance": stats["distance_sum"] / neighbor_count if neighbor_count else 0.0,
        "top_secondary_structure": format_counter(stats["secondary_structure"], top_n),
        "top_neighbor_residues": format_counter(stats["neighbor_residues"], top_n),
        "example_pdu_ids": ";".join(str(example["pdu_id"]) for example in stats["examples"][:top_n]),
        "example_sites": ";".join(format_site(example) for example in stats["examples"][:top_n]),
    }


def write_examples(path, aa, space, cluster_stats, top_n):
    rows = []
    for cluster_id, stats in sorted(cluster_stats.items()):
        for example in stats["examples"][:top_n]:
            rows.append(
                {
                    "aa": aa,
                    "space": space,
                    "cluster": cluster_id,
                    **example,
                }
            )
    write_dict_rows(path, rows)


def write_combined_summary(path, rows):
    write_dict_rows(path, rows)


def write_dict_rows(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def format_counter(counter, top_n):
    total = sum(counter.values())
    if total == 0:
        return ""
    return ";".join(f"{key}:{value / total:.3f}" for key, value in counter.most_common(top_n))


def format_site(example):
    chain = example["chain_id"] or "."
    return f"{example['pdb_id']}:{chain}:{example['residue_number']}:{example['residue_name']}"


if __name__ == "__main__":
    main()
