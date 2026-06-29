#!/usr/bin/env python3
import argparse
import logging
import random
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from kmers.logging_utils import add_logging_args, configure_logging


ENTRIES_INDEX_URL = "https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx"
PDB_ENTRY_URL = "https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/{middle}/pdb{pdb_id}.ent.gz"
logger = logging.getLogger("download_pdb_subset")


def main():
    parser = argparse.ArgumentParser(description="Download a bounded subset of PDB .ent.gz files.")
    parser.add_argument("--count", type=int, default=1000, help="Number of PDB entries to download.")
    parser.add_argument("--output-dir", default="pdb", help="Destination directory used by the pipeline.")
    parser.add_argument("--ids-file", default=None, help="Optional file containing one PDB id per line.")
    parser.add_argument("--seed", type=int, default=1, help="Random seed used when sampling from the PDB index.")
    parser.add_argument("--workers", type=int, default=8, help="Concurrent downloads.")
    add_logging_args(parser)
    args = parser.parse_args()
    configure_logging(args.log_file, args.log_level)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pdb_ids = load_ids(args.ids_file)
    if args.ids_file is None:
        rng = random.Random(args.seed)
        rng.shuffle(pdb_ids)

    pdb_ids = pdb_ids[: args.count]
    if not pdb_ids:
        raise SystemExit("No PDB ids selected.")
    logger.info("Downloading %s PDB entries into %s with %s workers", len(pdb_ids), output_dir, args.workers)

    downloaded = 0
    failed = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(download_entry, pdb_id, output_dir): pdb_id for pdb_id in pdb_ids}
        for idx, future in enumerate(as_completed(futures), start=1):
            pdb_id = futures[future]
            try:
                path = future.result()
                downloaded += 1
                logger.info("[%s/%s] %s -> %s", idx, len(pdb_ids), pdb_id, path)
            except Exception as exc:
                failed.append((pdb_id, str(exc)))
                logger.warning("[%s/%s] %s failed: %s", idx, len(pdb_ids), pdb_id, exc)

    logger.info("Downloaded %s files into %s", downloaded, output_dir)
    if failed:
        print(f"Failed downloads: {len(failed)}", file=sys.stderr)
        for pdb_id, reason in failed[:20]:
            print(f"{pdb_id}: {reason}", file=sys.stderr)
        raise SystemExit(1)


def load_ids(ids_file):
    if ids_file:
        with open(ids_file) as handle:
            return [line.strip().lower() for line in handle if line.strip() and not line.startswith("#")]

    with urllib.request.urlopen(ENTRIES_INDEX_URL, timeout=60) as response:
        lines = response.read().decode("utf-8", errors="replace").splitlines()

    ids = []
    for line in lines:
        if not line or line.startswith("IDCODE") or line.startswith("-"):
            continue
        pdb_id = line.split()[0].lower()
        if len(pdb_id) == 4:
            ids.append(pdb_id)
    return ids


def download_entry(pdb_id, output_dir):
    pdb_id = pdb_id.lower()
    middle = pdb_id[1:3]
    entry_dir = output_dir / middle
    entry_dir.mkdir(parents=True, exist_ok=True)
    output_path = entry_dir / f"pdb{pdb_id}.ent.gz"
    if output_path.exists() and output_path.stat().st_size > 0:
        return output_path

    url = PDB_ENTRY_URL.format(middle=middle, pdb_id=pdb_id)
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            output_path.write_bytes(response.read())
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"{url} returned HTTP {exc.code}") from exc

    if output_path.stat().st_size == 0:
        output_path.unlink(missing_ok=True)
        raise RuntimeError("Downloaded file is empty")
    return output_path


if __name__ == "__main__":
    main()
