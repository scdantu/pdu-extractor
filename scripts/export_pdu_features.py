#!/usr/bin/env python3
import argparse
import logging
import sqlite3
from pathlib import Path

import numpy as np

from kmers.residue_classes import AA_ORDER, RESIDUE_CLASS_ORDER, residue_class
from kmers.logging_utils import add_logging_args, configure_logging

SS_ORDER = ("H", "E", "C")
logger = logging.getLogger("export_pdu_features")


def main():
    parser = argparse.ArgumentParser(description="Export fixed-length PDU feature matrices by reference amino acid.")
    parser.add_argument("--db", default="pdu_output/pdus.sqlite", help="PDU SQLite database.")
    parser.add_argument("--out-dir", default="analysis/features", help="Output directory for .npz files.")
    parser.add_argument("--aa", default=None, help="Optional single reference amino acid to export, e.g. A.")
    parser.add_argument("--bin-width", type=float, default=1.0, help="Radial shell width in Angstroms.")
    parser.add_argument("--radius", type=float, default=15.0, help="Maximum PDU radius in Angstroms.")
    parser.add_argument("--min-pdus", type=int, default=25, help="Skip AA classes with fewer PDUs.")
    parser.add_argument(
        "--residue-encoding",
        choices=("aa", "chemical"),
        default="aa",
        help="Encode neighbors by exact amino-acid identity or broad physicochemical class.",
    )
    add_logging_args(parser)
    args = parser.parse_args()
    configure_logging(args.log_file, args.log_level)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    distance_bins = np.arange(0.0, args.radius + args.bin_width, args.bin_width)
    n_bins = len(distance_bins) - 1
    residue_labels = list(AA_ORDER) if args.residue_encoding == "aa" else list(RESIDUE_CLASS_ORDER)
    feature_names = build_feature_names(distance_bins, residue_labels)

    conn = sqlite3.connect(args.db)
    aa_counts = conn.execute(
        """
        SELECT reference_residue_one_letter, COUNT(*)
        FROM pdu
        GROUP BY reference_residue_one_letter
        """
    ).fetchall()
    if args.aa:
        aa_counts = [(aa, count) for aa, count in aa_counts if aa == args.aa]
        if not aa_counts:
            raise SystemExit(f"No PDUs found for reference amino acid: {args.aa}")

    for aa, count in aa_counts:
        if aa not in AA_ORDER or count < args.min_pdus:
            continue
        pdu_ids = [row[0] for row in conn.execute("SELECT id FROM pdu WHERE reference_residue_one_letter = ? ORDER BY id", (aa,))]
        matrix = np.zeros((len(pdu_ids), len(feature_names)), dtype=np.float32)
        pdu_id_to_row = {pdu_id: idx for idx, pdu_id in enumerate(pdu_ids)}

        placeholders = ",".join("?" for _ in pdu_ids)
        rows = conn.execute(
            f"""
            SELECT pdu_id, residue_one_letter, secondary_structure, distance_angstrom
            FROM pdu_residue
            WHERE pdu_id IN ({placeholders})
            """,
            pdu_ids,
        )
        for pdu_id, residue, secondary_structure, distance in rows:
            distance = float(distance)
            if distance > args.radius:
                continue
            label = residue if args.residue_encoding == "aa" else residue_class(residue)
            try:
                residue_idx = residue_labels.index(label)
            except ValueError:
                continue
            ss_idx = SS_ORDER.index(secondary_structure) if secondary_structure in SS_ORDER else SS_ORDER.index("C")
            bin_idx = min(int(distance // args.bin_width), n_bins - 1)
            col = ((residue_idx * len(SS_ORDER)) + ss_idx) * n_bins + bin_idx
            matrix[pdu_id_to_row[pdu_id], col] += 1.0

        row_sums = matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        matrix = matrix / row_sums

        np.savez_compressed(
            out_dir / f"pdu_features_{aa}.npz",
            X=matrix,
            pdu_ids=np.array(pdu_ids, dtype=np.int64),
            feature_names=np.array(feature_names),
            reference_aa=np.array([aa]),
            distance_bins=distance_bins,
            residue_encoding=np.array([args.residue_encoding]),
            residue_labels=np.array(residue_labels),
            aa_order=np.array(list(AA_ORDER)),
            ss_order=np.array(SS_ORDER),
        )
        logger.info("%s: wrote %s PDUs x %s features", aa, len(pdu_ids), matrix.shape[1])


def build_feature_names(distance_bins, residue_labels):
    names = []
    for aa in residue_labels:
        for ss in SS_ORDER:
            for start, end in zip(distance_bins[:-1], distance_bins[1:]):
                names.append(f"{aa}_{ss}_{start:.1f}-{end:.1f}A")
    return names


if __name__ == "__main__":
    main()
