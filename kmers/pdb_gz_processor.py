import gzip
import sqlite3
import subprocess
import time
from pathlib import Path

from kmers.calculate_kmer import calculate_kmers
from kmers.pdb_data import PDBData
from kmers.pdu import AnnotationStore, PDUWriter, calculate_pdus, parse_pdb_secondary_structure


class GZProcessor:
    def __init__(
        self,
        db_path,
        process_dir,
        out_uniprot_dir,
        out_pdbs_dir,
        handle_all_pdbs,
        pdu_db_path=None,
        annotation_csv=None,
        pdu_radius_angstrom=15.0,
    ):
        self.db_path = db_path
        self.process_dir = process_dir
        self.out_uniprot_dir = out_uniprot_dir
        self.out_pdbs_dir = out_pdbs_dir
        self.handle_all_pdbs = handle_all_pdbs
        self.pdu_radius_angstrom = pdu_radius_angstrom
        self.annotation_store = AnnotationStore(annotation_csv)
        self.pdu_writer = PDUWriter(pdu_db_path) if pdu_db_path else None

        if not self.handle_all_pdbs:
            self.conn = sqlite3.connect(f'file:{self.db_path}?mode=ro', uri=True)

        self.codes = {'SUCCESS': 0}
        self.max_pdb_count = 1  # to avoid division by zero
        self.cur_pdb_count = 0

    def process_gz_file(self, gz_file):

        # 1. extract coordinates
        try:
            parsed_pdb = self.extract_coordinates(gz_file)
            self.codes['SUCCESS'] += 1
        except Exception as e:
            self.codes[str(e)] = self.codes.get(str(e), 0) + 1
            return

        # 2. parse data
        pdb_data = PDBData(parsed_pdb)

        # 3. find matching uniprot entry, reject if not found
        uniprot_id = None
        if not self.handle_all_pdbs:
            uniprot_id = self.get_matching_uniprot_entry(pdb_data)
            if uniprot_id is None:
                self.codes['SUCCESS'] -= 1
                self.codes['NO_UNIPROT_ID'] = self.codes.get('NO_UNIPROT_ID', 0) + 1
                return

        # print(f'{pdb_id} -> {uniprot_id}')

        kmers = calculate_kmers(pdb_data)
        # 4. write data to pdb & uniprot files
        self._write_pdb_file(pdb_data.pdb_id, kmers)
        if not self.handle_all_pdbs:
            self._append_to_uniprot_file(uniprot_id, pdb_data.pdb_id, pdb_data)  # noqa

        if self.pdu_writer:
            secondary_structure = parse_pdb_secondary_structure(gz_file)
            pdus = calculate_pdus(
                pdb_data,
                pdb_secondary_structure=secondary_structure,
                annotation_store=self.annotation_store,
                radius_angstrom=self.pdu_radius_angstrom,
            )
            self.pdu_writer.write_pdus(pdus)

        # 5. pass kmers to natural set parser
        # TODO

    @staticmethod
    def extract_coordinates(gz_file):
        proc = subprocess.Popen(['bin/extract_pdb_coordinates'],
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with gzip.open(gz_file, 'r') as f_in:
            pdb_str = f_in.read()

        parsed_pdb, err = proc.communicate(input=pdb_str)

        if proc.returncode != 0:
            raise Exception(err.decode('utf-8').splitlines()[0])

        return parsed_pdb

    def process_files(self):
        """
        Process all files in the process_dir
        :return:
        """
        self.max_pdb_count = count_files(self.process_dir, '*.ent.gz')
        print(f'Processing {self.max_pdb_count} PDB files files...')

        time_start = time.time()

        for gz_file in Path(self.process_dir).rglob('*.ent.gz'):
            self.process_gz_file(gz_file)
            self.cur_pdb_count += 1

            if self.cur_pdb_count % 100 == 0:
                self.print_progress()

        self.print_progress()
        time_end = time.time()

        self.print_codes()
        print(f'\nCompleted in {time_end - time_start:.2f} seconds')
        if self.pdu_writer:
            self.pdu_writer.close()

    def print_progress(self):
        print(f'\r{self.cur_pdb_count:<{len(str(self.max_pdb_count))}} / {self.max_pdb_count}, '
              f'{self.cur_pdb_count / self.max_pdb_count:.1%}', end='')

    def print_codes(self):
        print()
        for k, v in self.codes.items():
            print(f'{k}: {v}')

    def get_matching_uniprot_entry(self, pdb_data):
        """
        Fetch sequences by uniprot ids first and then check for sequence match.
        """
        sequence = pdb_data.residue_sequence_parsed
        pdb_uniprot_ids = pdb_data.uniprot_ids
        first_residue_number = pdb_data.first_residue_number

        cur = self.conn.cursor()

        # Prepare IDs for the query
        placeholders = ', '.join(['?'] * len(pdb_uniprot_ids))
        query = f'SELECT id, sequence FROM sequences WHERE id IN ({placeholders})'

        cur.execute(query, tuple(pdb_uniprot_ids))
        ids_and_sequences = cur.fetchall()

        if len(ids_and_sequences) == 0:
            return None

        matched_uniprot_ids = []
        # Check for exact match if first_residue_number is not 0
        if first_residue_number != 0:
            # Find sequences that exactly match the input sequence starting and ending at the given points
            matched_uniprot_ids = [entry[0] for entry in ids_and_sequences if sequence ==
                                   entry[1][first_residue_number - 1:first_residue_number + len(sequence) - 1]]

        # If no exact matches or first_residue_number is 0, find sequences that have the input sequence as a substring
        if not matched_uniprot_ids:
            matched_uniprot_ids = [entry[0] for entry in ids_and_sequences if sequence in entry[1]]

        all_matches = [x for x in pdb_uniprot_ids if x in matched_uniprot_ids]

        # if len(all_matches) > 1:
        #     print(f'Found multiple matches for {sequence}: {all_matches}')
        return all_matches[0] if len(all_matches) > 0 else None

    def _write_pdb_file(self, pdb_id: str, kmers: list[str]):
        with open(f'{self.out_pdbs_dir}/{pdb_id}.kmers', 'w') as f_out:
            for kmer in kmers:
                f_out.write(f'{kmer}\n')

    def _append_to_uniprot_file(self, uniprot_id: str, pdb_id: str, pdb_data: PDBData):
        if not Path(f'{self.out_uniprot_dir}/{uniprot_id}.info').exists():
            with open(f'{self.out_uniprot_dir}/{uniprot_id}.info', 'w') as f_out:
                f_out.write(f'>{uniprot_id}\n')

        with open(f'{self.out_uniprot_dir}/{uniprot_id}.info', 'a') as f_out:
            f_out.write(f'{pdb_id} {pdb_data.resolution} {len(pdb_data.residue_sequence_parsed)} '
                        f'{pdb_data.residue_sequence_parsed}\n')


def count_files(directory='.', extension='*'):
    count = 0
    for _ in Path(directory).rglob(extension):
        count += 1
    return count
