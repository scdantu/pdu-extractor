#!/usr/bin/env python3
import argparse
import csv
import logging
import math
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

from kmers.residue_classes import AA_ORDER, RESIDUE_CLASS_ORDER, residue_class
from kmers.logging_utils import add_logging_args, configure_logging


SS_ORDER = ("H", "E", "C")
logger = logging.getLogger("analyze_distance_variability")


def main():
    parser = argparse.ArgumentParser(description="Measure residue variability as a function of distance from r_i.")
    parser.add_argument("--db", default="pdu_output/pdus.sqlite", help="PDU SQLite database.")
    parser.add_argument("--out-dir", default="analysis/variability", help="Output directory.")
    parser.add_argument("--bin-width", type=float, default=1.0, help="Distance bin width in Angstroms.")
    parser.add_argument("--radius", type=float, default=15.0, help="Maximum distance in Angstroms.")
    parser.add_argument("--aa", default=None, help="Optional reference amino acid to analyze.")
    parser.add_argument("--include-self", action="store_true", help="Include the reference residue at distance 0.")
    add_logging_args(parser)
    args = parser.parse_args()
    configure_logging(args.log_file, args.log_level)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_bins = int(math.ceil(args.radius / args.bin_width))
    stats = defaultdict(new_bin_stats)
    reference_counts = Counter()

    conn = sqlite3.connect(args.db)
    query = """
        SELECT p.reference_residue_one_letter,
               r.residue_one_letter,
               r.secondary_structure,
               r.distance_angstrom
        FROM pdu p
        JOIN pdu_residue r ON r.pdu_id = p.id
    """
    params = []
    if args.aa:
        query += " WHERE p.reference_residue_one_letter = ?"
        params.append(args.aa)

    for ref_aa, residue, secondary_structure, distance in conn.execute(query, params):
        if ref_aa not in AA_ORDER or residue not in AA_ORDER:
            continue
        distance = float(distance)
        if not args.include_self and distance == 0.0:
            continue
        if distance > args.radius:
            continue
        bin_idx = min(int(distance // args.bin_width), n_bins - 1)
        bin_stats = stats[(ref_aa, bin_idx)]
        bin_stats["neighbor_count"] += 1
        bin_stats["residue_counts"][residue] += 1
        bin_stats["class_counts"][residue_class(residue)] += 1
        bin_stats["ss_counts"][secondary_structure if secondary_structure in SS_ORDER else "C"] += 1

    reference_query = "SELECT reference_residue_one_letter, COUNT(*) FROM pdu"
    reference_params = []
    if args.aa:
        reference_query += " WHERE reference_residue_one_letter = ?"
        reference_params.append(args.aa)
    reference_query += " GROUP BY reference_residue_one_letter"
    for ref_aa, count in conn.execute(reference_query, reference_params):
        reference_counts[ref_aa] = count

    rows = []
    for ref_aa in sorted(reference_counts):
        cumulative_neighbors = 0
        for bin_idx in range(n_bins):
            bin_stats = stats[(ref_aa, bin_idx)]
            neighbor_count = bin_stats["neighbor_count"]
            cumulative_neighbors += neighbor_count
            residue_entropy = entropy(bin_stats["residue_counts"])
            class_entropy = entropy(bin_stats["class_counts"])
            ss_entropy = entropy(bin_stats["ss_counts"])
            rows.append(
                {
                    "reference_aa": ref_aa,
                    "distance_start": bin_idx * args.bin_width,
                    "distance_end": (bin_idx + 1) * args.bin_width,
                    "distance_mid": (bin_idx + 0.5) * args.bin_width,
                    "pdu_count": reference_counts[ref_aa],
                    "neighbor_count": neighbor_count,
                    "mean_neighbors_per_pdu": neighbor_count / reference_counts[ref_aa],
                    "cumulative_mean_neighbors_per_pdu": cumulative_neighbors / reference_counts[ref_aa],
                    "residue_entropy": residue_entropy,
                    "residue_entropy_norm": normalized_entropy(residue_entropy, len(AA_ORDER)),
                    "residue_specificity": 1.0 - normalized_entropy(residue_entropy, len(AA_ORDER)),
                    "class_entropy": class_entropy,
                    "class_entropy_norm": normalized_entropy(class_entropy, len(RESIDUE_CLASS_ORDER)),
                    "class_specificity": 1.0 - normalized_entropy(class_entropy, len(RESIDUE_CLASS_ORDER)),
                    "secondary_structure_entropy": ss_entropy,
                    "secondary_structure_entropy_norm": normalized_entropy(ss_entropy, len(SS_ORDER)),
                    "top_residues": format_counter(bin_stats["residue_counts"], 8),
                    "top_classes": format_counter(bin_stats["class_counts"], 8),
                    "top_secondary_structure": format_counter(bin_stats["ss_counts"], 3),
                }
            )

    output_path = out_dir / "distance_variability.csv"
    write_rows(output_path, rows)
    logger.info("Wrote %s", output_path)


def new_bin_stats():
    return {
        "neighbor_count": 0,
        "residue_counts": Counter(),
        "class_counts": Counter(),
        "ss_counts": Counter(),
    }


def entropy(counter):
    total = sum(counter.values())
    if total == 0:
        return 0.0
    value = 0.0
    for count in counter.values():
        probability = count / total
        value -= probability * math.log2(probability)
    return value


def normalized_entropy(value, n_categories):
    if n_categories <= 1:
        return 0.0
    return value / math.log2(n_categories)


def format_counter(counter, top_n):
    total = sum(counter.values())
    if total == 0:
        return ""
    return ";".join(f"{key}:{value / total:.3f}" for key, value in counter.most_common(top_n))


def write_rows(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
