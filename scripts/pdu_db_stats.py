#!/usr/bin/env python3
import argparse
import csv
import logging
import sqlite3
from pathlib import Path

from kmers.logging_utils import add_logging_args, configure_logging
from kmers.residue_classes import residue_class


logger = logging.getLogger("pdu_db_stats")


def main():
    parser = argparse.ArgumentParser(description="Generate summary statistics from a PDU SQLite database.")
    parser.add_argument("--db", default="pdu_output/pdus.sqlite", help="PDU SQLite database.")
    parser.add_argument("--out-dir", default="analysis/db_stats", help="Directory for CSV/Markdown outputs.")
    parser.add_argument(
        "--include-residue-stats",
        action="store_true",
        help="Also scan pdu_residue for neighbor, secondary-structure, and distance-shell stats.",
    )
    parser.add_argument("--bin-width", type=float, default=1.0, help="Distance shell width for residue stats.")
    parser.add_argument("--top-pdbs", type=int, default=50, help="Number of largest PDB entries to write.")
    add_logging_args(parser)
    args = parser.parse_args()
    configure_logging(args.log_file, args.log_level)

    db_path = Path(args.db)
    if not db_path.exists() or db_path.stat().st_size == 0:
        raise SystemExit(f"Database does not exist or is empty: {db_path}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    stats = collect_basic_stats(conn, db_path)
    write_basic_outputs(conn, out_dir, args.top_pdbs)

    if args.include_residue_stats:
        logger.info("Collecting pdu_residue stats; this can take time for full PDB-scale databases")
        stats.update(collect_residue_stats(conn, out_dir, args.bin_width))

    summary_path = out_dir / "summary.md"
    write_summary(summary_path, stats, args.include_residue_stats)
    logger.info("Wrote %s", summary_path)
    print_summary(stats, args.include_residue_stats)


def collect_basic_stats(conn, db_path):
    logger.info("Collecting basic PDU stats from %s", db_path)
    pdu_count = scalar(conn, "SELECT COUNT(*) FROM pdu")
    distinct_pdb_count = scalar(conn, "SELECT COUNT(DISTINCT pdb_id) FROM pdu")
    unique_pdu_hash_count = scalar(conn, "SELECT COUNT(DISTINCT pdu_hash) FROM pdu")
    annotated_family_count = scalar(conn, "SELECT COUNT(DISTINCT family_id) FROM pdu WHERE family_id IS NOT NULL")
    return {
        "db_path": str(db_path),
        "db_size_gb": db_path.stat().st_size / (1024 ** 3),
        "pdu_count": pdu_count,
        "distinct_pdb_count": distinct_pdb_count,
        "unique_pdu_hash_count": unique_pdu_hash_count,
        "annotated_family_count": annotated_family_count,
        "mean_pdus_per_pdb": pdu_count / distinct_pdb_count if distinct_pdb_count else 0.0,
    }


def write_basic_outputs(conn, out_dir, top_pdbs):
    write_query(
        conn,
        out_dir / "pdu_count_by_reference_aa.csv",
        """
        SELECT reference_residue_one_letter AS reference_aa, COUNT(*) AS pdu_count
        FROM pdu
        GROUP BY reference_residue_one_letter
        ORDER BY pdu_count DESC
        """,
    )
    write_query(
        conn,
        out_dir / "pdu_count_by_pdb_top.csv",
        """
        SELECT pdb_id, COUNT(*) AS pdu_count
        FROM pdu
        GROUP BY pdb_id
        ORDER BY pdu_count DESC
        LIMIT ?
        """,
        (top_pdbs,),
    )
    write_query(
        conn,
        out_dir / "pdu_per_pdb_distribution.csv",
        """
        SELECT
            MIN(pdu_count) AS min_pdus_per_pdb,
            AVG(pdu_count) AS mean_pdus_per_pdb,
            MAX(pdu_count) AS max_pdus_per_pdb
        FROM (
            SELECT pdb_id, COUNT(*) AS pdu_count
            FROM pdu
            GROUP BY pdb_id
        )
        """,
    )


def collect_residue_stats(conn, out_dir, bin_width):
    pdu_residue_count = scalar(conn, "SELECT COUNT(*) FROM pdu_residue")
    pdu_count = scalar(conn, "SELECT COUNT(*) FROM pdu")
    write_query(
        conn,
        out_dir / "neighbor_secondary_structure.csv",
        """
        SELECT secondary_structure, COUNT(*) AS n
        FROM pdu_residue
        GROUP BY secondary_structure
        ORDER BY n DESC
        """,
    )
    residue_rows = query_rows(
        conn,
        """
        SELECT residue_one_letter, COUNT(*) AS n
        FROM pdu_residue
        GROUP BY residue_one_letter
        ORDER BY n DESC
        """,
    )
    write_rows(out_dir / "neighbor_residue_counts.csv", ["residue_one_letter", "n"], residue_rows)

    class_counts = {}
    for residue, count in residue_rows:
        class_counts[residue_class(residue)] = class_counts.get(residue_class(residue), 0) + count
    write_rows(
        out_dir / "neighbor_residue_class_counts.csv",
        ["residue_class", "n"],
        sorted(class_counts.items(), key=lambda item: item[1], reverse=True),
    )
    write_query(
        conn,
        out_dir / "neighbor_distance_shell_counts.csv",
        """
        SELECT
            CAST(distance_angstrom / ? AS INTEGER) * ? AS distance_start,
            (CAST(distance_angstrom / ? AS INTEGER) + 1) * ? AS distance_end,
            COUNT(*) AS n
        FROM pdu_residue
        WHERE distance_angstrom > 0
        GROUP BY CAST(distance_angstrom / ? AS INTEGER)
        ORDER BY distance_start
        """,
        (bin_width, bin_width, bin_width, bin_width, bin_width),
    )
    return {
        "pdu_residue_count": pdu_residue_count,
        "mean_neighbors_per_pdu": pdu_residue_count / pdu_count if pdu_count else 0.0,
    }


def write_summary(path, stats, include_residue_stats):
    lines = [
        "# PDU Database Summary",
        "",
        f"- Database: `{stats['db_path']}`",
        f"- Size: {stats['db_size_gb']:.2f} GB",
        f"- PDB entries represented: {stats['distinct_pdb_count']:,}",
        f"- PDUs: {stats['pdu_count']:,}",
        f"- Unique PDU hashes: {stats['unique_pdu_hash_count']:,}",
        f"- Mean PDUs/PDB: {stats['mean_pdus_per_pdb']:.2f}",
        f"- Annotated families: {stats['annotated_family_count']:,}",
    ]
    if include_residue_stats:
        lines.extend(
            [
                f"- PDU-residue rows: {stats['pdu_residue_count']:,}",
                f"- Mean neighbors/PDU: {stats['mean_neighbors_per_pdu']:.2f}",
            ]
        )
    path.write_text("\n".join(lines) + "\n")


def print_summary(stats, include_residue_stats):
    print(f"PDB entries represented: {stats['distinct_pdb_count']:,}")
    print(f"PDUs: {stats['pdu_count']:,}")
    print(f"Unique PDU hashes: {stats['unique_pdu_hash_count']:,}")
    print(f"Mean PDUs/PDB: {stats['mean_pdus_per_pdb']:.2f}")
    if include_residue_stats:
        print(f"PDU-residue rows: {stats['pdu_residue_count']:,}")
        print(f"Mean neighbors/PDU: {stats['mean_neighbors_per_pdu']:.2f}")


def scalar(conn, query, params=()):
    return conn.execute(query, params).fetchone()[0]


def query_rows(conn, query, params=()):
    return conn.execute(query, params).fetchall()


def write_query(conn, path, query, params=()):
    cursor = conn.execute(query, params)
    headers = [description[0] for description in cursor.description]
    write_rows(path, headers, cursor.fetchall())
    logger.info("Wrote %s", path)


def write_rows(path, headers, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
