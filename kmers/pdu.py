import csv
import gzip
import hashlib
import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path

from kmers.neighbors import query_radius
from kmers.pdb_data import PDBData


@dataclass(frozen=True)
class ResidueAnnotation:
    secondary_structure: str | None = None
    sasa: float | None = None
    family_id: str | None = None


class AnnotationStore:
    """
    Residue-level annotations keyed by PDB id, chain, residue number, and insertion code.
    CSV input is intentionally permissive so mdCATH exports can be normalized without
    rewriting this pipeline.
    """

    def __init__(self, csv_path: str | None = None):
        self._annotations = {}
        if csv_path:
            self._load_csv(csv_path)

    def _load_csv(self, csv_path: str):
        with open(csv_path, newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                pdb_id = _first_present(row, "pdb_id", "pdb", "structure_id", "entry_id")
                chain_id = _first_present(row, "chain_id", "chain", "auth_chain_id", "label_asym_id")
                residue_number = _first_present(row, "residue_number", "resnum", "residue_id", "auth_seq_id", "label_seq_id")
                if not pdb_id or chain_id is None or residue_number is None:
                    continue

                insertion_code = _first_present(row, "insertion_code", "icode", "pdbx_pdb_ins_code") or ""
                secondary_structure = _first_present(row, "secondary_structure", "ss", "dssp", "secstruct")
                sasa_value = _first_present(row, "sasa", "sasa_value", "solvent_accessible_surface_area", "rsa")
                family_id = _first_present(row, "family_id", "cath_id", "cath_domain", "superfamily", "protein_family")

                self._annotations[_annotation_key(pdb_id, chain_id, residue_number, insertion_code)] = ResidueAnnotation(
                    secondary_structure=secondary_structure,
                    sasa=_to_float_or_none(sasa_value),
                    family_id=family_id,
                )

    def get(self, pdb_id, chain_id, residue_number, insertion_code="") -> ResidueAnnotation:
        return self._annotations.get(
            _annotation_key(pdb_id, chain_id, residue_number, insertion_code),
            ResidueAnnotation(),
        )


class PDUWriter:
    def __init__(self, db_path: str):
        self.db_path = db_path
        Path(db_path).parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(db_path)
        self._create_schema()

    def close(self):
        self.conn.close()

    def _create_schema(self):
        self.conn.executescript(
            """
            PRAGMA journal_mode=WAL;

            CREATE TABLE IF NOT EXISTS pdu (
                id INTEGER PRIMARY KEY,
                pdu_hash TEXT NOT NULL,
                pdb_id TEXT NOT NULL,
                family_id TEXT,
                reference_residue_index INTEGER NOT NULL,
                reference_chain_id TEXT,
                reference_residue_number INTEGER,
                reference_insertion_code TEXT,
                reference_residue_name TEXT,
                reference_residue_one_letter TEXT,
                radius_angstrom REAL NOT NULL,
                UNIQUE (pdb_id, reference_residue_index)
            );

            CREATE TABLE IF NOT EXISTS pdu_residue (
                pdu_id INTEGER NOT NULL,
                neighbor_rank INTEGER NOT NULL,
                residue_index INTEGER NOT NULL,
                chain_id TEXT,
                residue_number INTEGER,
                insertion_code TEXT,
                residue_name TEXT,
                residue_one_letter TEXT NOT NULL,
                distance_angstrom REAL NOT NULL,
                secondary_structure TEXT,
                sasa REAL,
                PRIMARY KEY (pdu_id, neighbor_rank),
                FOREIGN KEY (pdu_id) REFERENCES pdu(id) ON DELETE CASCADE
            );

            CREATE INDEX IF NOT EXISTS idx_pdu_hash ON pdu (pdu_hash);
            CREATE INDEX IF NOT EXISTS idx_pdu_family ON pdu (family_id);
            CREATE INDEX IF NOT EXISTS idx_pdu_residue_lookup
                ON pdu_residue (chain_id, residue_number, insertion_code);

            CREATE VIEW IF NOT EXISTS unique_pdu_counts AS
            SELECT pdu_hash, COUNT(*) AS occurrences, COUNT(DISTINCT family_id) AS family_count
            FROM pdu
            GROUP BY pdu_hash;
            """
        )
        self.conn.commit()

    def write_pdus(self, pdus: list[dict]):
        with self.conn:
            for pdu in pdus:
                existing = self.conn.execute(
                    """
                    SELECT id FROM pdu
                    WHERE pdb_id = ? AND reference_residue_index = ?
                    """,
                    (pdu["pdb_id"], pdu["reference_residue_index"]),
                ).fetchone()
                if existing:
                    pdu_id = existing[0]
                    self.conn.execute(
                        """
                        UPDATE pdu
                        SET pdu_hash = ?, family_id = ?, reference_chain_id = ?,
                            reference_residue_number = ?, reference_insertion_code = ?,
                            reference_residue_name = ?, reference_residue_one_letter = ?,
                            radius_angstrom = ?
                        WHERE id = ?
                        """,
                        (
                            pdu["pdu_hash"],
                            pdu["family_id"],
                            pdu["reference_chain_id"],
                            pdu["reference_residue_number"],
                            pdu["reference_insertion_code"],
                            pdu["reference_residue_name"],
                            pdu["reference_residue_one_letter"],
                            pdu["radius_angstrom"],
                            pdu_id,
                        ),
                    )
                    self.conn.execute("DELETE FROM pdu_residue WHERE pdu_id = ?", (pdu_id,))
                else:
                    cursor = self.conn.execute(
                        """
                        INSERT INTO pdu (
                            pdu_hash, pdb_id, family_id, reference_residue_index,
                            reference_chain_id, reference_residue_number, reference_insertion_code,
                            reference_residue_name, reference_residue_one_letter, radius_angstrom
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            pdu["pdu_hash"],
                            pdu["pdb_id"],
                            pdu["family_id"],
                            pdu["reference_residue_index"],
                            pdu["reference_chain_id"],
                            pdu["reference_residue_number"],
                            pdu["reference_insertion_code"],
                            pdu["reference_residue_name"],
                            pdu["reference_residue_one_letter"],
                            pdu["radius_angstrom"],
                        ),
                    )
                    pdu_id = cursor.lastrowid
                self.conn.executemany(
                    """
                    INSERT INTO pdu_residue (
                        pdu_id, neighbor_rank, residue_index, chain_id, residue_number,
                        insertion_code, residue_name, residue_one_letter, distance_angstrom,
                        secondary_structure, sasa
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    [
                        (
                            pdu_id,
                            residue["neighbor_rank"],
                            residue["residue_index"],
                            residue["chain_id"],
                            residue["residue_number"],
                            residue["insertion_code"],
                            residue["residue_name"],
                            residue["residue_one_letter"],
                            residue["distance_angstrom"],
                            residue["secondary_structure"],
                            residue["sasa"],
                        )
                        for residue in pdu["residues"]
                    ],
                )


def parse_pdb_secondary_structure(gz_file: str | Path) -> dict:
    secondary_structure = {}
    with gzip.open(gz_file, "rt", errors="replace") as handle:
        for line in handle:
            record = line[:6]
            if record.startswith("HELIX "):
                chain_id = line[19].strip()
                start = _safe_int(line[21:25])
                start_icode = line[25].strip()
                end = _safe_int(line[33:37])
                end_icode = line[37].strip()
                _mark_range(secondary_structure, chain_id, start, start_icode, end, end_icode, "H")
            elif record.startswith("SHEET "):
                chain_id = line[21].strip()
                start = _safe_int(line[22:26])
                start_icode = line[26].strip()
                end = _safe_int(line[33:37])
                end_icode = line[37].strip()
                _mark_range(secondary_structure, chain_id, start, start_icode, end, end_icode, "E")
    return secondary_structure


def calculate_pdus(
    pdb_data: PDBData,
    pdb_secondary_structure: dict | None = None,
    annotation_store: AnnotationStore | None = None,
    radius_angstrom: float = 15.0,
    distance_precision: int = 3,
    sasa_precision: int = 3,
) -> list[dict]:
    residues = pdb_data.residue_list
    coordinates = pdb_data.coordinates
    indices, distances = query_radius(coordinates, radius_angstrom, sort_results=True)
    pdb_secondary_structure = pdb_secondary_structure or {}
    annotation_store = annotation_store or AnnotationStore()

    pdus = []
    for reference_idx, (neighbor_indices, neighbor_distances) in enumerate(zip(indices, distances)):
        reference_family = None
        pdu_residues = []
        signature_residues = []

        for rank, (neighbor_idx, distance) in enumerate(zip(neighbor_indices, neighbor_distances)):
            chain_id = pdb_data.chain_id_list[neighbor_idx]
            residue_number = pdb_data.residue_number_list[neighbor_idx]
            insertion_code = pdb_data.insertion_code_list[neighbor_idx]
            annotation = annotation_store.get(pdb_data.pdb_id, chain_id, residue_number, insertion_code)
            secondary_structure = annotation.secondary_structure
            if secondary_structure is None:
                secondary_structure = pdb_secondary_structure.get((chain_id or "", residue_number, insertion_code), "C")
            if reference_family is None and annotation.family_id:
                reference_family = annotation.family_id

            residue_name = pdb_data.residue_name_list[neighbor_idx] or residues[neighbor_idx]
            rounded_distance = round(float(distance), distance_precision)
            rounded_sasa = None if annotation.sasa is None else round(annotation.sasa, sasa_precision)
            pdu_residues.append(
                {
                    "neighbor_rank": rank,
                    "residue_index": int(neighbor_idx),
                    "chain_id": chain_id,
                    "residue_number": residue_number,
                    "insertion_code": insertion_code,
                    "residue_name": residue_name,
                    "residue_one_letter": residues[neighbor_idx],
                    "distance_angstrom": rounded_distance,
                    "secondary_structure": secondary_structure,
                    "sasa": rounded_sasa,
                }
            )
            signature_residues.append(
                [
                    residues[neighbor_idx],
                    rounded_distance,
                    secondary_structure,
                    rounded_sasa,
                ]
            )

        pdu_hash = hashlib.sha256(json.dumps(signature_residues, separators=(",", ":")).encode("utf-8")).hexdigest()
        pdus.append(
            {
                "pdu_hash": pdu_hash,
                "pdb_id": pdb_data.pdb_id,
                "family_id": reference_family,
                "reference_residue_index": reference_idx,
                "reference_chain_id": pdb_data.chain_id_list[reference_idx],
                "reference_residue_number": pdb_data.residue_number_list[reference_idx],
                "reference_insertion_code": pdb_data.insertion_code_list[reference_idx],
                "reference_residue_name": pdb_data.residue_name_list[reference_idx] or residues[reference_idx],
                "reference_residue_one_letter": residues[reference_idx],
                "radius_angstrom": radius_angstrom,
                "residues": pdu_residues,
            }
        )

    return pdus


def _annotation_key(pdb_id, chain_id, residue_number, insertion_code=""):
    return (str(pdb_id).upper(), str(chain_id or "").strip(), int(residue_number), str(insertion_code or "").strip())


def _first_present(row: dict, *names):
    for name in names:
        if name in row and row[name] not in (None, ""):
            return row[name]
    return None


def _safe_int(value):
    value = value.strip()
    return int(value) if value else None


def _to_float_or_none(value):
    if value in (None, ""):
        return None
    return float(value)


def _mark_range(secondary_structure, chain_id, start, start_icode, end, end_icode, code):
    if start is None or end is None:
        return
    chain_id = chain_id or ""
    start_icode = start_icode or ""
    end_icode = end_icode or ""
    for residue_number in range(start, end + 1):
        insertion_code = start_icode if residue_number == start else end_icode if residue_number == end else ""
        secondary_structure[(chain_id, residue_number, insertion_code)] = code
