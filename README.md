# k-mers
Description and purpose of the project will come at a later point.

## Installation

### Prerequisites
`g++`: through gcc installation

`gzip / gunzip`: usually pre-installed

`python3`: installed and part of path

optional: `sqlite3`: Debian/Ubuntu: `apt install sqlite3`; macOS: `brew install sqlite3`

### Preparation

Clone the repo
```
git clone https://github.com/cemiu/kmers.git && cd kmers
```
Two folders need to be populated:
- `pdb` (required)
- `uniprotkb` (optional)
### pdb
The `pdb` folder has to contain experimental PDB files in the .ent.gz format.

Instructions for downloading can be found here:

https://www.rcsb.org/docs/programmatic-access/file-download-services

https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/

Alternatively a mirror can be found here: https://pycom.brunel.ac.uk/misc/pdb_2023-07-28.tar (42 GB)

Once downloaded they have to be placed in the `pdb` folder **without** being decompressed.
It does not matter whether they are divided; e.g. `pdb/file.ent.gz` or `pdb/<folder>/file.ent.gz`.

### uniprotkb
**Optionally**, the `uniprotkb` folder can be populated with the uniprotkb fasta files.
This is only needed, if the PDBs should be associated to a Protein. If this is not required, **skip this step**.

Files:
- `uniprot_sprot.fasta.gz` (400 MB after processing)
  - Use only Swiss-Prot, has the majority of PDB coverage
- Optionally, `uniprot_trembl.fasta.gz` (250 GB after processing)
  - Use TrEMBL to match more PDBs; might be useful for max. coverage

The latter might result in (slightly) more PDBs which can be associated to a Protein. The difference is expected to be trivial.

Place the files in the `uniprotkb` folder without decompressing them.

Once `run.sh` has been executed and the database `uniprotkb/uniprot_sequences.db` has been created, the files can be deleted.

### Running

To run the script, execute:
```
./run.sh
```

This will
- Compile C++ binaries, if they don't exist
- Ask whether to process the uniprotkb files
  - If yes, create the database `uniprotkb/uniprot_sequences.db`, if it doesn't exist
- Ask for k-mer size (default: k=12)
- Extracts 3d k-mer from the PDBs (into `pdb_output` folder)
- Extracts k-mer of length k into `kmer.txt`, along with frequency

### Testing With A Small PDB Subset

You do not need the full PDB mirror for local testing. Download a bounded random subset:

```
python3 scripts/download_pdb_subset.py --count 1000 --output-dir pdb --workers 8
```

For a deterministic hand-picked set, create a text file with one PDB id per line and run:

```
python3 scripts/download_pdb_subset.py --ids-file pdb_ids.txt --count 1000 --output-dir pdb
```

The script stores files in the same divided `.ent.gz` layout expected by the pipeline, for example
`pdb/hh/pdb4hhb.ent.gz`.

### Protein Dynamic Units

A protein dynamic unit (PDU) can be written for every reference residue by extracting all C-alpha residues within
15 Angstroms / 1.5 nm. The SQLite output stores the reference residue, each neighboring residue, distance from the
reference residue, secondary-structure annotation, SASA value, and a hash that can be counted across families.

Run the pipeline directly to create the PDU database:

```
PYTHONPATH=. python3 kmers/pipeline.py \
  --handle_all_pdbs true \
  --pdu_db pdu_output/pdus.sqlite \
  --annotation_csv path/to/mdcath_residue_annotations.csv
```

Or use the wrapper script:

```
scripts/extract_pdus.sh
```

Useful options:

```
PDU_DB=pdu_output/pdus_10A.sqlite RADIUS_ANGSTROM=10.0 scripts/extract_pdus.sh
RESET_OUTPUT=true scripts/extract_pdus.sh
ANNOTATION_CSV=path/to/annotations.csv scripts/extract_pdus.sh
LOG_DIR=logs/full_pdb scripts/extract_pdus.sh
```

Most Python scripts support:

```
--log-file logs/script_name.log
--log-level INFO
```

Use `DEBUG` for more detail, or `WARNING` to keep server logs quieter.

Generate database statistics after PDU extraction:

```
PYTHONPATH=. python scripts/pdu_db_stats.py \
  --db pdu_output/pdus.sqlite \
  --out-dir analysis/db_stats \
  --log-file logs/db_stats.log
```

For heavier neighbor-level statistics, add:

```
--include-residue-stats
```

This scans `pdu_residue`, so expect it to take longer on full PDB-scale databases.

Export a large SQLite database into per-amino-acid analysis stores:

```
PYTHONPATH=. python scripts/export_analysis_store.py \
  --db /path/to/pdus.sqlite \
  --out-dir analysis/per_aa_sqlite \
  --format split-sqlite \
  --log-file logs/export_split_sqlite.log
```

Or export partitioned Parquet:

```
PYTHONPATH=. python scripts/export_analysis_store.py \
  --db /path/to/pdus.sqlite \
  --out-dir analysis/parquet \
  --format parquet \
  --log-file logs/export_parquet.log
```

For cluster jobs, the safest parallel pattern is one job per amino-acid class:

```
PYTHONPATH=. python scripts/export_analysis_store.py \
  --db /path/to/pdus.sqlite \
  --out-dir analysis/per_aa_sqlite \
  --format split-sqlite \
  --aa A
```

The annotation CSV is optional, but SASA values require it. The loader accepts common column names:

- PDB id: `pdb_id`, `pdb`, `structure_id`, or `entry_id`
- Chain: `chain_id`, `chain`, `auth_chain_id`, or `label_asym_id`
- Residue number: `residue_number`, `resnum`, `residue_id`, `auth_seq_id`, or `label_seq_id`
- Secondary structure: `secondary_structure`, `ss`, `dssp`, or `secstruct`
- SASA: `sasa`, `sasa_value`, `solvent_accessible_surface_area`, or `rsa`
- Family: `family_id`, `cath_id`, `cath_domain`, `superfamily`, or `protein_family`

If no secondary-structure column is provided, HELIX/SHEET records from the PDB file are used (`H`, `E`, otherwise `C`).

Useful SQLite queries:

```
-- Number of unique PDUs
SELECT COUNT(*) FROM unique_pdu_counts;

-- PDUs shared by the most families
SELECT pdu_hash, occurrences, family_count
FROM unique_pdu_counts
ORDER BY family_count DESC, occurrences DESC
LIMIT 20;
```

### Static PDB Clustering

Before processing mdCATH, use the static PDB-derived PDUs to check whether residue environments cluster.
The first pass uses one model per reference amino-acid class.

Create the Python environment:

```
scripts/setup_pypdu.sh
source pypdu/bin/activate
```

Equivalent manual setup:

```
python3 -m venv pypdu
source pypdu/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements.txt
python scripts/check_environment.py
```

On macOS, the normal PyTorch wheel can use Apple Silicon GPU acceleration through MPS when available. On Linux
cluster nodes with NVIDIA GPUs, use the CUDA-specific requirements file that matches the cluster CUDA module, for
example:

```
python -m pip install -r requirements-cuda-cu124.txt
```

Export fixed-length feature vectors from the SQLite database:

```
python3 scripts/export_pdu_features.py \
  --db pdu_output/pdus.sqlite \
  --out-dir analysis/features
```

Each vector is a normalized radial composition of neighboring residue type, secondary structure, and distance shell.
By default this creates 900 features per PDU: 20 amino acids x 3 secondary-structure labels x 15 one-Angstrom shells.

To encode neighboring residues by broad physicochemical class instead of exact amino-acid identity:

```
python3 scripts/export_pdu_features.py \
  --db pdu_output/pdus.sqlite \
  --out-dir analysis/features_chemical \
  --residue-encoding chemical
```

The chemical classes are `hydrophobic`, `aromatic`, `polar`, `positive`, `negative`, and `special` (`G/P`).

Measure how residue variability changes with distance from the reference residue:

```
python3 scripts/analyze_distance_variability.py \
  --db pdu_output/pdus.sqlite \
  --out-dir analysis/variability \
  --bin-width 1.0 \
  --radius 15.0
```

The output reports per-distance-shell entropy and specificity. `residue_specificity` and `class_specificity` are
`1 - normalized_entropy`, so lower values mean more variable and likely less informative shells. The reference residue
itself is excluded by default; pass `--include-self` if you want the distance-zero bin.

Train autoencoders and write 2D coordinates:

```
python3 scripts/train_pdu_autoencoder.py \
  --features-dir analysis/features \
  --out-dir analysis/embeddings \
  --epochs 40
```

To process only one reference amino acid:

```
python3 scripts/train_pdu_autoencoder.py --aa A --epochs 40
```

If `umap-learn` is installed, the CSV output uses UMAP coordinates. If it is not installed, the script still writes
embeddings and uses a PCA fallback for the two plotting columns.

Cluster the embeddings with HDBSCAN:

```
python3 scripts/cluster_pdu_embeddings.py \
  --embeddings-dir analysis/embeddings \
  --out-dir analysis/clusters \
  --space latent \
  --min-cluster-size 25 \
  --min-samples 5
```

The recommended default is to cluster the autoencoder latent vectors (`--space latent`) and use UMAP only for plotting.
To cluster the 2D UMAP coordinates instead:

```
python3 scripts/cluster_pdu_embeddings.py --space umap --min-cluster-size 50 --min-samples 10
```

Outputs:

```
analysis/clusters/pdu_clusters_A_latent.csv
analysis/clusters/cluster_summary_latent.csv
```

Plot clusters on the 2D UMAP coordinates:

```
python3 scripts/plot_pdu_clusters.py \
  --embeddings-dir analysis/embeddings \
  --clusters-dir analysis/clusters \
  --out-dir analysis/plots \
  --aa A \
  --space umap
```

The plot always uses the `x` and `y` embedding columns for display. The `--space` option chooses which cluster labels
to overlay: labels produced from latent vectors or labels produced from UMAP coordinates.

Summarize biological content of clusters:

```
python3 scripts/summarize_pdu_clusters.py \
  --db pdu_output/pdus.sqlite \
  --clusters-dir analysis/clusters_umap \
  --out-dir analysis/cluster_summaries \
  --space umap
```

Outputs include per-cluster size, mean neighbor count, secondary-structure composition, neighbor amino-acid
composition, neighbor physicochemical class composition, and representative PDB residue examples.

Summarize external CATH, EC, GO, or function annotations per cluster:

```
python3 scripts/summarize_cluster_annotations.py \
  --db pdu_output/pdus.sqlite \
  --clusters-dir analysis/clusters_umap \
  --annotations annotations/cath_ec_go.csv \
  --out-dir analysis/annotation_summaries \
  --space umap
```

The annotation CSV can be PDB-level, chain-level, or residue-range-level. Required column: `pdb_id`. Optional
coordinate columns: `chain_id`, `start_residue`, `end_residue`. All other columns are treated as annotation fields,
for example `cath_id`, `ec_number`, `go_terms`, `function`, or `ligand`.
