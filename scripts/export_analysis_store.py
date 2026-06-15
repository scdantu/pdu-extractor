#!/usr/bin/env python3
import argparse
import logging
import sqlite3
from pathlib import Path

from kmers.logging_utils import add_logging_args, configure_logging
from kmers.residue_classes import AA_ORDER


logger = logging.getLogger("export_analysis_store")


def main():
    parser = argparse.ArgumentParser(description="Export a large PDU SQLite DB into analysis-friendly stores.")
    parser.add_argument("--db", required=True, help="Source PDU SQLite database.")
    parser.add_argument("--out-dir", required=True, help="Output directory.")
    parser.add_argument("--aa", default=None, help="Optional single reference amino acid to export.")
    parser.add_argument(
        "--format",
        choices=("split-sqlite", "parquet"),
        default="split-sqlite",
        help="Export per-AA SQLite databases or partitioned Parquet files.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
    add_logging_args(parser)
    args = parser.parse_args()
    configure_logging(args.log_file, args.log_level)

    db_path = Path(args.db)
    if not db_path.exists() or db_path.stat().st_size == 0:
        raise SystemExit(f"Source database does not exist or is empty: {db_path}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    aa_list = [args.aa] if args.aa else list(AA_ORDER)
    if args.format == "split-sqlite":
        export_split_sqlite(db_path, out_dir, aa_list, args.overwrite)
    else:
        export_parquet(db_path, out_dir, aa_list, args.overwrite)


def export_split_sqlite(db_path, out_dir, aa_list, overwrite):
    source = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    for aa in aa_list:
        output_path = out_dir / f"pdus_{aa}.sqlite"
        if output_path.exists():
            if overwrite:
                output_path.unlink()
            else:
                logger.info("%s exists; skipping", output_path)
                continue

        logger.info("Exporting %s to %s", aa, output_path)
        source.execute("ATTACH DATABASE ? AS out", (str(output_path),))
        try:
            source.executescript(
                f"""
                CREATE TABLE out.pdu AS
                SELECT *
                FROM main.pdu
                WHERE reference_residue_one_letter = '{aa}';

                CREATE TABLE out.pdu_residue AS
                SELECT r.*
                FROM main.pdu_residue r
                JOIN out.pdu p ON p.id = r.pdu_id;

                CREATE INDEX out.idx_pdu_id ON pdu(id);
                CREATE INDEX out.idx_pdu_hash ON pdu(pdu_hash);
                CREATE INDEX out.idx_pdu_pdb ON pdu(pdb_id);
                CREATE INDEX out.idx_pdu_residue_pdu_id ON pdu_residue(pdu_id);
                CREATE INDEX out.idx_pdu_residue_distance ON pdu_residue(distance_angstrom);
                """
            )
            source.execute("DETACH DATABASE out")
        except Exception:
            source.execute("DETACH DATABASE out")
            raise
        logger.info("Wrote %s", output_path)


def export_parquet(db_path, out_dir, aa_list, overwrite):
    try:
        import duckdb
    except ModuleNotFoundError as exc:
        raise SystemExit("Parquet export requires duckdb. Install with: conda install -c conda-forge duckdb") from exc

    conn = duckdb.connect()
    conn.execute("INSTALL sqlite")
    conn.execute("LOAD sqlite")
    conn.execute("ATTACH ? AS pdu_db (TYPE sqlite, READ_ONLY)", [str(db_path)])

    for aa in aa_list:
        aa_dir = out_dir / f"reference_aa={aa}"
        if aa_dir.exists() and any(aa_dir.iterdir()):
            if overwrite:
                for path in aa_dir.glob("*.parquet"):
                    path.unlink()
            else:
                logger.info("%s exists; skipping", aa_dir)
                continue
        aa_dir.mkdir(parents=True, exist_ok=True)
        pdu_path = aa_dir / "pdu.parquet"
        residue_path = aa_dir / "pdu_residue.parquet"
        logger.info("Exporting %s to %s", aa, aa_dir)
        conn.execute(
            """
            COPY (
                SELECT *
                FROM pdu_db.pdu
                WHERE reference_residue_one_letter = ?
            ) TO ? (FORMAT PARQUET, COMPRESSION ZSTD)
            """,
            [aa, str(pdu_path)],
        )
        conn.execute(
            """
            COPY (
                SELECT r.*
                FROM pdu_db.pdu_residue r
                JOIN pdu_db.pdu p ON p.id = r.pdu_id
                WHERE p.reference_residue_one_letter = ?
            ) TO ? (FORMAT PARQUET, COMPRESSION ZSTD)
            """,
            [aa, str(residue_path)],
        )
        logger.info("Wrote %s and %s", pdu_path, residue_path)


if __name__ == "__main__":
    main()
