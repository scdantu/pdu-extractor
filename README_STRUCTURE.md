# PDU Extractor: Project Structure

## Overview

Rational Protein Design framework using PDU (Protein Dynamic Units) clustering. This project validates that introducing residue X at position Y requires specific structural neighborhoods to be "satisfied."

## Directory Structure

```
pdu-extractor/
├── README.md                          # Original project README
├── README_STRUCTURE.md                # This file - directory guide
│
├── docs/                              # All documentation
│   ├── QUICKSTART.md                  # Quick start for ARCHER2
│   ├── WORKFLOW_PYCOM.md              # PyCoM enrichment pipeline
│   ├── WORKFLOW_ENRICHMENT.md         # Enrichment analysis guide
│   ├── WORKFLOW_VALIDATION.md         # Rigorous validation pipeline
│   ├── WORKFLOW_ARCHER2.md            # ARCHER2 distributed training
│   ├── ANALYSIS_RESULTS.md            # Cluster analysis results
│   └── RIGOROUS_VALIDATION_PIPELINE.md
│
├── kmers/                             # Core PDU extraction library
│   ├── __init__.py
│   ├── pdu.py                         # PDU definition & representation
│   ├── neighbors.py                   # Neighboring atom calculation
│   ├── residue_classes.py             # Amino acid categorization
│   ├── pipeline.py                    # Main extraction pipeline
│   └── ...
│
├── bin/                               # Executable scripts (organized by phase)
│
│   ├── extract/                       # Phase 1: PDU Feature Extraction
│   │   ├── extract_pdus.sh            # Extract PDUs from PDB files
│   │   ├── export_pdu_features.py     # Export 900-dim feature vectors
│   │   └── download_pdb_subset.py     # Download PDB structures
│   │
│   ├── train/                         # Phase 2: Representation Learning
│   │   ├── train_pdu_autoencoder.py   # Train 900→16 dimensional AE (GPU)
│   │   ├── train_and_cluster_archer2.py  # Combined AE + HDBSCAN (HPC)
│   │   └── train_pdu_autoencoder_archer2.py
│   │
│   ├── cluster/                       # Phase 3: Local Clustering
│   │   ├── cluster_pdu_embeddings.py  # HDBSCAN clustering (CPU)
│   │   └── cluster_pdu_embeddings_gpu.py  # HDBSCAN clustering (GPU/RAPIDS)
│   │
│   ├── analyze/                       # Phase 4: Cluster Analysis
│   │   ├── analyze_cluster_outputs.py        # Summary stats across AAs
│   │   ├── analyze_clusters_pycom_direct.py  # Biological enrichment
│   │   ├── validate_clusters_vs_central_ss.py # Secondary structure validation
│   │   ├── analyze_cluster_profiles.py       # Cluster composition
│   │   ├── analyze_cath_ec_enrichment.py     # Structural/enzymatic enrichment
│   │   ├── summarize_pdu_clusters.py         # Generate summaries
│   │   └── pdu_db_stats.py                   # Database statistics
│   │
│   ├── validate/                      # Phase 5: Rigorous Validation
│   │   ├── hyperparameter_grid_parallel.py   # Grid search (parallel array job)
│   │   ├── hyperparameter_grid_search.py     # Grid search (sequential)
│   │   ├── validate_stability_parallel.py    # Bootstrap + perturbation
│   │   ├── validate_clusters_stability.py    # Stability testing
│   │   ├── aggregate_grid_results.py         # Combine grid results
│   │   ├── aggregate_stability_results.py    # Combine stability results
│   │   └── feature_ablation_archer2.py       # Feature importance
│   │
│   ├── enrich/                        # Phase 5b: Biological Enrichment
│   │   ├── fetch_pycom_pdb_local.py   # Build PyCoM reference (local DB)
│   │   ├── fetch_pycom_pdb_mapping.py # Build PyCoM reference (REST API)
│   │   └── build_pdb_metadata_via_pycom.py
│   │
│   └── Other utilities
│       ├── inspect_npz.py             # Inspect embedding files
│       └── export_analysis_store.py   # Export analysis data
│
├── pipelines/                         # Pipeline runners
│   ├── run_aa_pipeline.sh             # Single AA processing (skip if done)
│   ├── run_full_pipeline.sh           # All 20 AAs with GPU distribution
│   ├── run_rigorous_validation.py     # Full validation pipeline
│   └── ... (setup scripts)
│
├── submit/                            # ARCHER2 batch job templates
│   ├── submit_archer2.slurm           # Basic training job
│   ├── submit_grid_search_parallel.slurm    # Parallel grid search (array job)
│   ├── submit_stability_parallel.slurm     # Parallel stability (array job)
│   ├── submit_validation_archer2.slurm     # Full validation job
│   └── setup_archer2_venv.sh          # Setup Python env on ARCHER2
│
├── data/                              # Data directories
│   ├── pdb/                           # PDB structure files
│   ├── uniprotkb/                     # UniProtKB metadata
│   └── analysis/                      # Generated analysis outputs
│
├── analysis/                          # Actual analysis outputs (generated)
│   ├── features/                      # 900-dim feature vectors (NPZ)
│   ├── embeddings/                    # 16-dim autoencoder latent space (NPZ)
│   ├── clusters_umap/                 # Cluster assignments + UMAP (CSV)
│   ├── cluster_summaries/             # Cluster statistics
│   ├── cath_ec_enrichment/            # Enrichment results
│   └── plots/                         # Visualizations
│
├── pdb_db/                            # Per-amino acid SQLite databases (generated)
│   ├── pdus_L.sqlite                  # Leucine PDUs
│   ├── pdus_A.sqlite                  # Alanine PDUs
│   └── ... (20 total, one per AA)
│
├── pycom_db/                          # PyCoM reference databases (downloaded)
│   ├── pycom.db                       # Protein metadata (700MB)
│   └── pycom.mat                      # Coevolution matrices (115GB, optional)
│
└── requirements.txt                   # Python dependencies
```

## Workflow: Rational Protein Design

### Step 1: Extract PDU Features
```bash
cd pipelines
./run_aa_pipeline.sh L   # For Leucine
# or
./run_full_pipeline.sh --parallel 3  # For all 20 AAs
```
**Output:** `analysis/features/*.npz` (900-dim vectors per AA)

### Step 2: Train Autoencoders & Cluster
```bash
# Local GPU
python3 bin/train/train_pdu_autoencoder.py --aa L --gpu

# Or on ARCHER2 (HPC)
cd submit
sbatch submit_archer2.slurm
```
**Output:** `analysis/embeddings/*.npz` + `analysis/clusters_umap/*.csv`

### Step 3: Analyze Clusters
```bash
python3 bin/analyze/analyze_cluster_outputs.py
python3 bin/analyze/validate_clusters_vs_central_ss.py --aa L
```
**Output:** Cluster statistics, quality metrics

### Step 4: Biological Enrichment (PyCoM)
**First time only - build mapping:**
```bash
python3 bin/enrich/fetch_pycom_pdb_local.py \
    --db-path pycom_db/pycom.db \
    --out-file pdb_pycom_mapping.csv
```

**Then enrichment analysis:**
```bash
python3 bin/analyze/analyze_clusters_pycom_direct.py \
    --clusters analysis/clusters_umap/pdu_clusters_L_umap.csv \
    --db pdb_db/pdus_L.sqlite \
    --pycom-mapping pdb_pycom_mapping.csv \
    --aa L \
    --out-dir analysis/enrichment
```
**Output:** CATH/EC/GO contingency tables, χ² significance

### Step 5: Rigorous Validation
**If enrichment shows signal:**
```bash
# Local
python3 bin/validate/hyperparameter_grid_search.py --aa L

# Or parallel on ARCHER2 (recommended)
cd submit
sbatch -N 2 --array=0-19 submit_grid_search_parallel.slurm
sbatch -N 1 --array=0-19 submit_stability_parallel.slurm
```
**Output:** Optimal hyperparameters, stability metrics, final report

## Key Files by Use Case

### "I want to run enrichment analysis"
1. Start: `bin/enrich/fetch_pycom_pdb_local.py`
2. Then: `bin/analyze/analyze_clusters_pycom_direct.py`
3. Doc: `docs/WORKFLOW_PYCOM.md`

### "I want to optimize clustering hyperparameters"
1. Start: `bin/validate/hyperparameter_grid_parallel.py`
2. Doc: `docs/WORKFLOW_VALIDATION.md`
3. Submit: `submit/submit_grid_search_parallel.slurm`

### "I want to run on ARCHER2 HPC"
1. Start: `submit/setup_archer2_venv.sh`
2. Doc: `docs/QUICKSTART.md`
3. Submit: `submit/submit_archer2.slurm`

### "I need to extract new PDUs"
1. Start: `pipelines/run_full_pipeline.sh`
2. Doc: `README.md`
3. Output: `analysis/features/*.npz`

## Database Locations

- **PDB structures:** `data/pdb/` (generated)
- **PDU databases:** `pdb_db/pdus_{AA}.sqlite` (generated, 7.4 GB for Leucine)
- **PyCoM reference:** `pycom_db/pycom.db` (downloaded, 700 MB)
- **Analysis results:** `analysis/` (generated)

## Scaling

- **1,000 proteins:** ~10 min extraction + 30 min training + 110 min clustering
- **200,000 proteins:** ~3-4 hours extraction + 90 min training + 4-6 hours clustering
- **Bottleneck:** HDBSCAN clustering (110 min for 861K embeddings on CPU)
- **Optimization:** Use RAPIDS CUML on GPU (10-50x speedup) or ARCHER2 distributed

## Current Status

- ✅ Features extracted for all 20 AAs (200K proteins)
- ✅ Embeddings trained (16-dim latent space)
- ✅ Clustering complete (LEU: 538 clusters, 88.6% noise at min_cluster_size=200)
- ✅ PyCoM enrichment pipeline ready
- ⏳ Enrichment analysis pending (PyCoM validation)
- ⏳ Rigorous validation (grid search + bootstrap)

## Questions?

- **Quick start:** See `docs/QUICKSTART.md`
- **PyCoM enrichment:** See `docs/WORKFLOW_PYCOM.md`
- **Full validation:** See `docs/WORKFLOW_VALIDATION.md`
- **ARCHER2 setup:** See `docs/QUICKSTART.md`
