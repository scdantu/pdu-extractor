#!/usr/bin/env python3
import argparse
import csv
import logging
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from kmers.logging_utils import add_logging_args, configure_logging

logger = logging.getLogger("summarize_cluster_annotations")


def main():
    parser = argparse.ArgumentParser(description="Summarize CATH/EC/GO/function annotations per PDU cluster.")
    parser.add_argument("--db", default="pdu_output/pdus.sqlite", help="PDU SQLite database.")
    parser.add_argument("--clusters-dir", default="analysis/clusters_umap", help="Directory with pdu_clusters_*_*.csv.")
    parser.add_argument("--annotations", required=True, help="CSV with PDB-level, chain-level, or residue-range annotations.")
    parser.add_argument("--out-dir", default="analysis/annotation_summaries", help="Output directory.")
    parser.add_argument("--aa", default=None, help="Optional single reference amino acid.")
    parser.add_argument("--space", choices=("latent", "umap"), default="umap", help="Cluster label source.")
    parser.add_argument("--top-n", type=int, default=10, help="Number of top annotations to report per field.")
    add_logging_args(parser)
    args = parser.parse_args()
    configure_logging(args.log_file, args.log_level)

    annotation_rows, annotation_fields = load_annotations(args.annotations)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cluster_paths = sorted(Path(args.clusters_dir).glob(f"pdu_clusters_*_{args.space}.csv"))
    if args.aa:
        cluster_paths = [Path(args.clusters_dir) / f"pdu_clusters_{args.aa}_{args.space}.csv"]
    if not cluster_paths:
        raise SystemExit(f"No cluster files found in {args.clusters_dir}")

    conn = sqlite3.connect(args.db)
    all_rows = []
    for cluster_path in cluster_paths:
        if not cluster_path.exists():
            logger.warning("Skipping missing cluster file: %s", cluster_path)
            continue
        aa = cluster_path.name.replace("pdu_clusters_", "").replace(f"_{args.space}.csv", "")
        clusters = pd.read_csv(cluster_path)
        clusters = clusters[clusters["cluster"] != -1].copy()
        if clusters.empty:
            continue
        pdu_to_cluster = dict(zip(clusters["pdu_id"].astype(int), clusters["cluster"].astype(int)))
        pdu_rows = fetch_pdu_sites(conn, sorted(pdu_to_cluster))
        summary = summarize_annotations(pdu_rows, pdu_to_cluster, annotation_rows, annotation_fields)

        rows = [format_summary_row(aa, args.space, cluster_id, stats, annotation_fields, args.top_n)
                for cluster_id, stats in sorted(summary.items())]
        write_dict_rows(out_dir / f"annotation_summary_{aa}_{args.space}.csv", rows)
        all_rows.extend(rows)
        logger.info("%s: wrote annotation summaries for %s clusters", aa, len(rows))

    write_dict_rows(out_dir / f"annotation_summary_all_{args.space}.csv", all_rows)


def load_annotations(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        annotation_fields = [
            field for field in reader.fieldnames
            if field not in {"pdb_id", "pdb", "entry_id", "chain_id", "chain", "start_residue", "end_residue"}
        ]
        for row in reader:
            pdb_id = first_present(row, "pdb_id", "pdb", "entry_id")
            if not pdb_id:
                continue
            chain_id = first_present(row, "chain_id", "chain")
            start = first_present(row, "start_residue", "start", "residue_start")
            end = first_present(row, "end_residue", "end", "residue_end")
            rows.append(
                {
                    "pdb_id": str(pdb_id).upper(),
                    "chain_id": normalize_optional(chain_id),
                    "start_residue": int(start) if start not in (None, "") else None,
                    "end_residue": int(end) if end not in (None, "") else None,
                    "annotations": {
                        field: split_values(row.get(field, ""))
                        for field in annotation_fields
                    },
                }
            )
    return rows, annotation_fields


def fetch_pdu_sites(conn, pdu_ids):
    placeholders = ",".join("?" for _ in pdu_ids)
    return conn.execute(
        f"""
        SELECT id, pdb_id, reference_chain_id, reference_residue_number
        FROM pdu
        WHERE id IN ({placeholders})
        """,
        pdu_ids,
    ).fetchall()


def summarize_annotations(pdu_rows, pdu_to_cluster, annotation_rows, annotation_fields):
    annotations_by_pdb = defaultdict(list)
    for annotation in annotation_rows:
        annotations_by_pdb[annotation["pdb_id"]].append(annotation)

    summary = defaultdict(lambda: {"pdu_count": 0, "matched_pdu_count": 0, "fields": defaultdict(Counter)})
    for pdu_id, pdb_id, chain_id, residue_number in pdu_rows:
        cluster_id = pdu_to_cluster[pdu_id]
        stats = summary[cluster_id]
        stats["pdu_count"] += 1
        matched_annotations = [
            annotation for annotation in annotations_by_pdb.get(str(pdb_id).upper(), [])
            if annotation_matches(annotation, chain_id, residue_number)
        ]
        if matched_annotations:
            stats["matched_pdu_count"] += 1
        for annotation in matched_annotations:
            for field in annotation_fields:
                for value in annotation["annotations"].get(field, []):
                    stats["fields"][field][value] += 1
    return summary


def annotation_matches(annotation, chain_id, residue_number):
    annotation_chain = annotation["chain_id"]
    if annotation_chain and annotation_chain != normalize_optional(chain_id):
        return False
    start = annotation["start_residue"]
    end = annotation["end_residue"]
    if start is None or end is None:
        return True
    if residue_number is None:
        return False
    return start <= int(residue_number) <= end


def format_summary_row(aa, space, cluster_id, stats, annotation_fields, top_n):
    pdu_count = stats["pdu_count"]
    row = {
        "aa": aa,
        "space": space,
        "cluster": cluster_id,
        "pdu_count": pdu_count,
        "matched_pdu_count": stats["matched_pdu_count"],
        "matched_fraction": stats["matched_pdu_count"] / pdu_count if pdu_count else 0.0,
    }
    for field in annotation_fields:
        row[f"top_{field}"] = format_counter(stats["fields"][field], top_n)
    return row


def first_present(row, *names):
    for name in names:
        if name in row and row[name] not in (None, ""):
            return row[name]
    return None


def normalize_optional(value):
    if value in (None, "", "."):
        return None
    return str(value).strip()


def split_values(value):
    if value in (None, ""):
        return []
    values = []
    for part in str(value).replace("|", ";").split(";"):
        part = part.strip()
        if part:
            values.append(part)
    return values


def format_counter(counter, top_n):
    total = sum(counter.values())
    if total == 0:
        return ""
    return ";".join(f"{key}:{value / total:.3f}" for key, value in counter.most_common(top_n))


def write_dict_rows(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
