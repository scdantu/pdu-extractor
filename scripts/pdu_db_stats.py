#!/usr/bin/env python3
import argparse
import csv
import logging
import sqlite3
from collections import Counter
from pathlib import Path

from kmers.logging_utils import add_logging_args, configure_logging
from kmers.residue_classes import residue_class


logger = logging.getLogger("pdu_db_stats")


def main():
    parser = argparse.ArgumentParser(description="Generate summary statistics from one PDU SQLite DB or a DB directory.")
    parser.add_argument("--db", default="pdu_output/pdus.sqlite", help="PDU SQLite database.")
    parser.add_argument("--db-dir", default=None, help="Directory containing per-AA SQLite databases.")
    parser.add_argument("--glob", default="*.sqlite", help="Glob used with --db-dir.")
    parser.add_argument("--out-dir", default="analysis/db_stats", help="Directory for CSV/Markdown outputs.")
    parser.add_argument(
        "--include-residue-stats",
        action="store_true",
        help="Also scan pdu_residue for neighbor, secondary-structure, and distance-shell stats.",
    )
    parser.add_argument("--bin-width", type=float, default=1.0, help="Distance shell width for residue stats.")
    parser.add_argument("--top-pdbs", type=int, default=50, help="Number of largest PDB entries to write.")
    parser.add_argument(
        "--exact-combined-hashes",
        action="store_true",
        help="For --db-dir, compute exact distinct PDU hashes across DBs. This can use substantial memory.",
    )
    add_logging_args(parser)
    args = parser.parse_args()
    configure_logging(args.log_file, args.log_level)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.db_dir:
        db_paths = sorted(Path(args.db_dir).glob(args.glob))
        db_paths = [path for path in db_paths if path.is_file() and path.stat().st_size > 0]
        if not db_paths:
            raise SystemExit(f"No non-empty SQLite databases found in {args.db_dir} matching {args.glob}")
        stats = collect_directory_stats(
            db_paths,
            out_dir,
            args.include_residue_stats,
            args.bin_width,
            args.top_pdbs,
            args.exact_combined_hashes,
        )
    else:
        db_path = Path(args.db)
        if not db_path.exists() or db_path.stat().st_size == 0:
            raise SystemExit(f"Database does not exist or is empty: {db_path}")

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


def collect_directory_stats(db_paths, out_dir, include_residue_stats, bin_width, top_pdbs, exact_combined_hashes):
    logger.info("Collecting combined stats from %s databases", len(db_paths))
    database_rows = []
    pdu_by_aa = Counter()
    pdu_by_pdb = Counter()
    all_pdb_ids = set()
    all_hashes = set() if exact_combined_hashes else None
    residue_count = 0
    ss_counts = Counter()
    neighbor_residue_counts = Counter()
    distance_shell_counts = Counter()

    combined = {
        "db_path": f"{len(db_paths)} databases",
        "db_size_gb": 0.0,
        "pdu_count": 0,
        "distinct_pdb_count": 0,
        "unique_pdu_hash_count": 0,
        "annotated_family_count": 0,
        "mean_pdus_per_pdb": 0.0,
    }

    for db_path in db_paths:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        basic = collect_basic_stats(conn, db_path)
        database_rows.append(
            [
                str(db_path),
                basic["db_size_gb"],
                basic["distinct_pdb_count"],
                basic["pdu_count"],
                basic["unique_pdu_hash_count"],
                basic["mean_pdus_per_pdb"],
            ]
        )
        combined["db_size_gb"] += basic["db_size_gb"]
        combined["pdu_count"] += basic["pdu_count"]
        combined["unique_pdu_hash_count"] += basic["unique_pdu_hash_count"]
        combined["annotated_family_count"] += basic["annotated_family_count"]

        pdu_by_aa.update(dict(query_rows(conn, """
            SELECT reference_residue_one_letter, COUNT(*)
            FROM pdu
            GROUP BY reference_residue_one_letter
        """)))
        pdu_by_pdb.update(dict(query_rows(conn, """
            SELECT pdb_id, COUNT(*)
            FROM pdu
            GROUP BY pdb_id
        """)))
        all_pdb_ids.update(row[0] for row in query_rows(conn, "SELECT DISTINCT pdb_id FROM pdu"))
        if all_hashes is not None:
            all_hashes.update(row[0] for row in query_rows(conn, "SELECT DISTINCT pdu_hash FROM pdu"))

        if include_residue_stats:
            residue_count += scalar(conn, "SELECT COUNT(*) FROM pdu_residue")
            ss_counts.update(dict(query_rows(conn, """
                SELECT secondary_structure, COUNT(*)
                FROM pdu_residue
                GROUP BY secondary_structure
            """)))
            neighbor_residue_counts.update(dict(query_rows(conn, """
                SELECT residue_one_letter, COUNT(*)
                FROM pdu_residue
                GROUP BY residue_one_letter
            """)))
            distance_shell_counts.update(dict(query_rows(conn, """
                SELECT CAST(distance_angstrom / ? AS INTEGER) * ? AS distance_start, COUNT(*)
                FROM pdu_residue
                WHERE distance_angstrom > 0
                GROUP BY CAST(distance_angstrom / ? AS INTEGER)
            """, (bin_width, bin_width, bin_width))))

    combined["distinct_pdb_count"] = len(all_pdb_ids)
    combined["mean_pdus_per_pdb"] = combined["pdu_count"] / combined["distinct_pdb_count"] if combined["distinct_pdb_count"] else 0.0
    if all_hashes is not None:
        combined["unique_pdu_hash_count"] = len(all_hashes)
    else:
        combined["unique_pdu_hash_count"] = f"{combined['unique_pdu_hash_count']} (sum across DBs)"

    write_rows(
        out_dir / "database_summary.csv",
        ["db_path", "db_size_gb", "distinct_pdb_count", "pdu_count", "unique_pdu_hash_count", "mean_pdus_per_pdb"],
        database_rows,
    )
    write_rows(out_dir / "pdu_count_by_reference_aa.csv", ["reference_aa", "pdu_count"], pdu_by_aa.most_common())
    write_rows(out_dir / "pdu_count_by_pdb_top.csv", ["pdb_id", "pdu_count"], pdu_by_pdb.most_common(top_pdbs))
    write_rows(
        out_dir / "pdu_per_pdb_distribution.csv",
        ["min_pdus_per_pdb", "mean_pdus_per_pdb", "max_pdus_per_pdb"],
        [[min(pdu_by_pdb.values()), combined["mean_pdus_per_pdb"], max(pdu_by_pdb.values())]] if pdu_by_pdb else [[0, 0, 0]],
    )

    if include_residue_stats:
        combined["pdu_residue_count"] = residue_count
        combined["mean_neighbors_per_pdu"] = residue_count / combined["pdu_count"] if combined["pdu_count"] else 0.0
        write_rows(out_dir / "neighbor_secondary_structure.csv", ["secondary_structure", "n"], ss_counts.most_common())
        write_rows(out_dir / "neighbor_residue_counts.csv", ["residue_one_letter", "n"], neighbor_residue_counts.most_common())
        class_counts = Counter()
        for residue, count in neighbor_residue_counts.items():
            class_counts[residue_class(residue)] += count
        write_rows(out_dir / "neighbor_residue_class_counts.csv", ["residue_class", "n"], class_counts.most_common())
        write_rows(
            out_dir / "neighbor_distance_shell_counts.csv",
            ["distance_start", "distance_end", "n"],
            [[start, start + bin_width, count] for start, count in sorted(distance_shell_counts.items())],
        )

    return combined


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
        f"- PDB entries represented: {format_count(stats['distinct_pdb_count'])}",
        f"- PDUs: {format_count(stats['pdu_count'])}",
        f"- Unique PDU hashes: {format_count(stats['unique_pdu_hash_count'])}",
        f"- Mean PDUs/PDB: {stats['mean_pdus_per_pdb']:.2f}",
        f"- Annotated families: {format_count(stats['annotated_family_count'])}",
    ]
    if include_residue_stats:
        lines.extend(
            [
                f"- PDU-residue rows: {format_count(stats['pdu_residue_count'])}",
                f"- Mean neighbors/PDU: {stats['mean_neighbors_per_pdu']:.2f}",
            ]
        )
    path.write_text("\n".join(lines) + "\n")


def print_summary(stats, include_residue_stats):
    print(f"PDB entries represented: {format_count(stats['distinct_pdb_count'])}")
    print(f"PDUs: {format_count(stats['pdu_count'])}")
    print(f"Unique PDU hashes: {format_count(stats['unique_pdu_hash_count'])}")
    print(f"Mean PDUs/PDB: {stats['mean_pdus_per_pdb']:.2f}")
    if include_residue_stats:
        print(f"PDU-residue rows: {format_count(stats['pdu_residue_count'])}")
        print(f"Mean neighbors/PDU: {stats['mean_neighbors_per_pdu']:.2f}")


def format_count(value):
    return f"{value:,}" if isinstance(value, int) else str(value)


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
