"""Database utilities for PDU extraction."""

import sqlite3
from pathlib import Path
from typing import List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


class PDUDatabase:
    """Wrapper for PDU SQLite database."""

    def __init__(self, db_path: str):
        """Initialize database connection.

        Args:
            db_path: Path to PDU SQLite database.
        """
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {db_path}")
        self.conn = None

    def connect(self):
        """Establish database connection."""
        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row

    def close(self):
        """Close database connection."""
        if self.conn:
            self.conn.close()

    def __enter__(self):
        """Context manager entry."""
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def get_aa_counts(self, aa: Optional[str] = None) -> List[Tuple[str, int]]:
        """Get PDU counts for each amino acid.

        Args:
            aa: Optional single amino acid to query.

        Returns:
            List of (amino_acid, count) tuples.
        """
        if not self.conn:
            self.connect()

        query = """
            SELECT reference_residue_one_letter, COUNT(*)
            FROM pdu
            GROUP BY reference_residue_one_letter
            ORDER BY reference_residue_one_letter
        """
        rows = self.conn.execute(query).fetchall()

        if aa:
            rows = [(a, c) for a, c in rows if a == aa]

        return rows

    def get_pdu_ids_for_aa(self, aa: str) -> List[int]:
        """Get all PDU IDs for a specific amino acid.

        Args:
            aa: Amino acid one-letter code.

        Returns:
            List of PDU IDs.
        """
        if not self.conn:
            self.connect()

        query = """
            SELECT id FROM pdu
            WHERE reference_residue_one_letter = ?
            ORDER BY id
        """
        rows = self.conn.execute(query, (aa,)).fetchall()
        return [row[0] for row in rows]

    def get_pdu_count(self, aa: Optional[str] = None) -> int:
        """Get total PDU count.

        Args:
            aa: Optional amino acid to filter by.

        Returns:
            Total PDU count.
        """
        if not self.conn:
            self.connect()

        if aa:
            query = "SELECT COUNT(*) FROM pdu WHERE reference_residue_one_letter = ?"
            return self.conn.execute(query, (aa,)).fetchone()[0]
        else:
            query = "SELECT COUNT(*) FROM pdu"
            return self.conn.execute(query).fetchone()[0]

    def get_pdu_neighbors(
        self,
        aa: str,
        max_distance: Optional[float] = None,
    ) -> List[Tuple]:
        """Get all PDU-neighbor relationships for an amino acid.

        Args:
            aa: Amino acid one-letter code.
            max_distance: Maximum distance to include (Angstroms).

        Returns:
            List of (pdu_id, residue_one_letter, secondary_structure, distance) tuples.
        """
        if not self.conn:
            self.connect()

        query = """
            SELECT pr.pdu_id, pr.residue_one_letter, pr.secondary_structure, pr.distance_angstrom
            FROM pdu_residue pr
            INNER JOIN pdu p ON pr.pdu_id = p.id
            WHERE p.reference_residue_one_letter = ?
        """
        params = [aa]

        if max_distance is not None:
            query += " AND pr.distance_angstrom <= ?"
            params.append(max_distance)

        rows = self.conn.execute(query, params).fetchall()
        return rows

    def get_schema(self) -> dict:
        """Get database schema information.

        Returns:
            Dictionary with table and column information.
        """
        if not self.conn:
            self.connect()

        schema = {}
        tables = self.conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()

        for (table_name,) in tables:
            columns = self.conn.execute(f"PRAGMA table_info({table_name})").fetchall()
            schema[table_name] = [col[1] for col in columns]

        return schema

    def get_db_stats(self) -> dict:
        """Get database statistics.

        Returns:
            Dictionary with database statistics.
        """
        if not self.conn:
            self.connect()

        stats = {}

        # Total PDUs
        stats["total_pdus"] = self.get_pdu_count()

        # PDUs per amino acid
        aa_counts = self.get_aa_counts()
        stats["pdus_per_aa"] = {aa: count for aa, count in aa_counts}

        # Total residue-PDU relationships
        total_neighbors = self.conn.execute(
            "SELECT COUNT(*) FROM pdu_residue"
        ).fetchone()[0]
        stats["total_neighbors"] = total_neighbors

        return stats


class BatchPDUFetcher:
    """Efficiently fetch large batches of PDU neighbors."""

    def __init__(self, db_path: str, batch_size: int = 1000):
        """Initialize batch fetcher.

        Args:
            db_path: Path to PDU database.
            batch_size: Number of PDUs to fetch at a time.
        """
        self.db = PDUDatabase(db_path)
        self.batch_size = batch_size
        self.db.connect()

    def __del__(self):
        """Cleanup on deletion."""
        self.db.close()

    def fetch_neighbors_for_pdu_batch(
        self,
        pdu_ids: List[int],
        max_distance: Optional[float] = None,
    ) -> dict:
        """Fetch neighbors for a batch of PDU IDs.

        Args:
            pdu_ids: List of PDU IDs.
            max_distance: Maximum distance to include.

        Returns:
            Dictionary mapping pdu_id -> list of neighbors.
        """
        neighbors = {pdu_id: [] for pdu_id in pdu_ids}

        # Build safe query with proper parameterization
        placeholders = ",".join("?" * len(pdu_ids))
        query = f"""
            SELECT pdu_id, residue_one_letter, secondary_structure, distance_angstrom
            FROM pdu_residue
            WHERE pdu_id IN ({placeholders})
        """
        params = pdu_ids

        if max_distance is not None:
            query += " AND distance_angstrom <= ?"
            params.append(max_distance)

        rows = self.db.conn.execute(query, params).fetchall()

        for row in rows:
            pdu_id = row[0]
            neighbors[pdu_id].append(row[1:])

        return neighbors
