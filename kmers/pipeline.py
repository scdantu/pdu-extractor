import argparse
import os
from pathlib import Path

from kmers.pdb_gz_processor import GZProcessor

"""
- 1. OS walk through all PDB gz files 
- 2. Extract the gz files
- 3. extract the PDB coordinates annotated with CA (with C++ script) (PDB_ID.coor)
-     3.1. If unsuccessful (return 1), delete extracted file & created file, next file
-     3.2 If successful, delete extracted file & continue to 4.
- 4. Calculate the Kmers
5. Write them to file (PDB_ID.kmers)
6. Repeat for all files

After batch processing all files:
# 7. Associate each PDB with a protein from the DB
# 8. Select the PDB with the highest resolution
...

9. Read all Kmers into memory (set), with frequency information
    9.1. Truncate to 12 residues
    9.2. Discard ones which are shorter
10. Save natural set where n=12 to disk, along with frequencies
11. Generate synthetic set. (all possible Kmers of n=12)
12. Calculate difference
13. Extract statistics on freq. dist of natural set, fraction appearing in synthetic set, etc.
"""


standard_error_codes = {
    "RESOLUTION_TOO_LOW",
    "MISSING_NON_TERMINAL_RESIDUES",
    "IS_NOT_PROTEIN",
    "NO_ALPHA_CARBON_ATOMS_FOUND",
    "EXCLUDE_SELENOCYSTEINE_AND_PYRROLYSINE",
    "NO_UNIPROT_ID"
}


def check_empty_directory(directory):
    if not os.path.exists(directory):
        return

    if os.listdir(directory):
        raise RuntimeError(f"Error: Output directory '{directory}' is not empty. Aborting.")


def check_file_exists(file_path):
    if not os.path.exists(file_path):
        raise RuntimeError(f"Error: File '{file_path}' does not exist. Aborting.")

    if not os.path.isfile(file_path):
        raise RuntimeError(f"Error: Path '{file_path}' is not a file. Aborting.")


def check_directory_exists(directory_path):
    if not os.path.exists(directory_path):
        raise RuntimeError(f"Error: Directory '{directory_path}' does not exist. Aborting.")

    if not os.path.isdir(directory_path):
        raise RuntimeError(f"Error: Path '{directory_path}' is not a directory. Aborting.")


def check_directory_contains_files(directory_path, extension=""):
    for f in Path(directory_path).rglob(f"*{extension}"):
        if f.is_file():
            return True
    if extension == "":
        raise RuntimeError(f"Error: Directory '{directory_path}' does not contain any files. Aborting.")
    else:
        raise RuntimeError(f"Error: Directory '{directory_path}' does not contain any files with extension "
                           f"'{extension}'. Aborting.")


def parse_bool(value):
    if isinstance(value, bool):
        return value
    normalized = value.lower()
    if normalized in {"true", "1", "yes", "y"}:
        return True
    if normalized in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError("Expected a boolean value")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process PDB files.')
    parser.add_argument('--handle_all_pdbs', required=True, type=parse_bool,
                        help='Set to True to handle all PDBs without checking for uniprot IDs')
    parser.add_argument('--pdu_db', default=None,
                        help='Optional SQLite output path for protein dynamic unit records')
    parser.add_argument('--annotation_csv', default=None,
                        help='Optional mdCATH-style residue annotation CSV with secondary structure, SASA, and family IDs')
    parser.add_argument('--pdu_radius_angstrom', default=15.0, type=float,
                        help='PDU neighborhood radius in Angstroms. 15 Angstroms equals 1.5 nm')
    args = parser.parse_args()

    if args.handle_all_pdbs not in [True, False]:
        raise ValueError("The --handle_all_pdbs argument must be set to either True or False.")

    db_path = os.path.expanduser('uniprotkb/uniprot_sequences.db')

    if not args.handle_all_pdbs:  # Only check for database if we need it
        try:
            check_file_exists(db_path)
        except RuntimeError:
            print("Error: Database of uniprot sequences not found. Run ./prepare.sh after downloading"
                  " uniprot_sprot.fasta.gz and (optionally) uniprot_trembl.fasta.gz into the uniprotkb dir."
                  "Disk space for sprot only is ~400MB, or ~250GB for both.")
            exit(1)

    process_dir = 'pdb'
    try:
        check_directory_exists(process_dir)
        check_directory_contains_files(process_dir, extension='.ent.gz')
    except RuntimeError:
        print("Error: Directory 'pdb' not found or does not contain any .ent.gz files. "
              "It should contain e.g. pdb/a0/1a00.ent.gz, pdb/a0/1a01.ent.gz, etc."
              "Download PDB files from ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/ "
              "and extract them into the pdb directory."
              "A mirror can be found at https://pycom.brunel.ac.uk/misc/ (42GB tar file * 2 = 84GB)")
        exit(1)

    output_dir = 'pdb_output'
    try:
        check_empty_directory(output_dir)
    except RuntimeError:
        print("Error: Output directory 'pdb_output' is not empty. Aborting."
              "If you want to re-run the script, delete the directory first.")
        exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        os.makedirs(os.path.join(output_dir, 'uniprot'))
        os.makedirs(os.path.join(output_dir, 'pdbs'))

    out_uniprot = os.path.join(output_dir, 'uniprot')
    out_pdbs = os.path.join(output_dir, 'pdbs')

    processor = GZProcessor(
        db_path,
        process_dir,
        out_uniprot,
        out_pdbs,
        args.handle_all_pdbs,
        pdu_db_path=args.pdu_db,
        annotation_csv=args.annotation_csv,
        pdu_radius_angstrom=args.pdu_radius_angstrom,
    )
    processor.process_files()
