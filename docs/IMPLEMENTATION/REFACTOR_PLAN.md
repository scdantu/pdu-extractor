# PDU-Extractor Refactoring Plan

## Current State
- 40 standalone Python scripts in `bin/` organized by category (extract, train, cluster, analyze, validate, enrich)
- No shared library structure → code duplication across scripts
- Hard-coded paths, inconsistent argument parsing, mixed responsibilities
- Not easily importable for Jupyter notebooks
- Difficult to add features (functional categories, conservation metrics) consistently across pipeline

## Target State

```
pdu-extractor/
├── pdusearch/                    # Main package (PEP 8 compliant, importable)
│   ├── __init__.py
│   ├── config.py                 # Centralized configuration management
│   ├── logging_utils.py          # Unified logging setup
│   │
│   ├── features/                 # Feature extraction & representation
│   │   ├── __init__.py
│   │   ├── extractor.py          # FeatureExtractor class (load from DB, build vectors)
│   │   ├── encodings.py          # Encoding strategies (20-aa, 5-functional, +conservation)
│   │   └── alignment.py          # Sequence alignment utilities (PyCoM integration)
│   │
│   ├── models/                   # ML models
│   │   ├── __init__.py
│   │   ├── autoencoder.py        # Autoencoder architecture, training logic
│   │   └── training.py           # TrainingConfig, train_autoencoder() function
│   │
│   ├── clustering/               # Clustering & post-processing
│   │   ├── __init__.py
│   │   ├── hdbscan.py            # HDBSCANClusterer class (wrapper with logging)
│   │   ├── metrics.py            # Cluster quality metrics (noise %, confidence)
│   │   └── utils.py              # Cluster I/O, filtering
│   │
│   ├── enrichment/               # Enrichment analysis
│   │   ├── __init__.py
│   │   ├── pycom.py              # PyCoM integration (query alignments, CATH, EC)
│   │   └── statistical.py        # Chi-square tests, enrichment stats
│   │
│   ├── validation/               # Validation pipelines
│   │   ├── __init__.py
│   │   ├── stability.py          # Bootstrap/perturbation validation
│   │   ├── grid_search.py        # Hyperparameter grid search
│   │   └── metrics.py            # Validation quality metrics
│   │
│   ├── io/                       # Data I/O
│   │   ├── __init__.py
│   │   ├── database.py           # PDU database utilities
│   │   ├── npy.py                # NPZ/NPY I/O helpers
│   │   └── csv.py                # CSV I/O, result aggregation
│   │
│   └── utils/                    # Utilities
│       ├── __init__.py
│       ├── paths.py              # Path management (analysis dir, temp dirs)
│       └── types.py              # Type definitions, constants (AA_ORDER, SS_ORDER, etc)
│
├── bin/                          # CLI entry points (thin wrappers)
│   ├── extract_features.py       # → pdusearch.features.FeatureExtractor
│   ├── train_model.py            # → pdusearch.models.train_autoencoder()
│   ├── cluster.py                # → pdusearch.clustering.HDBSCANClusterer
│   ├── analyze_compression.py    # → pdusearch.models, pdusearch.clustering
│   ├── enrich.py                 # → pdusearch.enrichment
│   ├── validate.py               # → pdusearch.validation
│   └── ... (consolidated from 40 scripts)
│
├── notebooks/                    # Jupyter notebooks (use pdusearch library)
│   ├── 01_feature_exploration.ipynb
│   ├── 02_model_training.ipynb
│   ├── 03_clustering_analysis.ipynb
│   ├── 04_enrichment_validation.ipynb
│   └── 05_compression_hypothesis.ipynb
│
├── tests/                        # Unit tests (pytest)
│   ├── test_features.py
│   ├── test_models.py
│   ├── test_clustering.py
│   └── test_enrichment.py
│
├── setup.py                      # Package config (install pdusearch in dev mode)
├── requirements.txt              # Dependencies
├── REFACTOR_PLAN.md             # This file
└── README.md                     # Updated with library usage examples
```

## Key Design Patterns

### 1. Feature Extraction Pipeline
```python
# Current: scattered across 2-3 scripts
# New: unified interface

from pdusearch.features import FeatureExtractor, Encoding
from pdusearch.config import Config

config = Config(db="per_aa_sqlite/pdus_L.sqlite", radius=15.0)
extractor = FeatureExtractor(config, encoding=Encoding.AA_20)
X, pdu_ids = extractor.extract_for_aa("L")

# Also available:
extractor_func = FeatureExtractor(config, encoding=Encoding.FUNCTIONAL_5)
X_func, _ = extractor_func.extract_for_aa("L")

extractor_cons = FeatureExtractor(config, encoding=Encoding.AA_20_WITH_CONSERVATION)
X_cons, _ = extractor_cons.extract_for_aa("L")
```

### 2. Model Training Pipeline
```python
from pdusearch.models import Autoencoder, train_autoencoder
from pdusearch.config import TrainingConfig

config = TrainingConfig(latent_dim=16, epochs=10, batch_size=32, device='cuda')
model, Z, losses = train_autoencoder(X, config)  # Returns model, embeddings, training history
```

### 3. Clustering Pipeline
```python
from pdusearch.clustering import HDBSCANClusterer, ClusterMetrics

clusterer = HDBSCANClusterer(min_cluster_size=200, min_samples=5)
labels, confidences = clusterer.cluster(Z)
metrics = ClusterMetrics(labels, confidences)
print(f"Noise: {metrics.noise_percent:.1f}%, Avg Confidence: {metrics.avg_confidence:.3f}")
```

### 4. Complete End-to-End (Notebook-Friendly)
```python
from pdusearch.features import FeatureExtractor, Encoding
from pdusearch.models import train_autoencoder, TrainingConfig
from pdusearch.clustering import HDBSCANClusterer, ClusterMetrics
import pandas as pd

# Extract
X, pdu_ids = FeatureExtractor(..., encoding=Encoding.AA_20).extract_for_aa("L")

# Train
model, Z, _ = train_autoencoder(X, TrainingConfig(latent_dim=16, epochs=10))

# Cluster
labels, conf = HDBSCANClusterer(min_cluster_size=200).cluster(Z)

# Analyze
metrics = ClusterMetrics(labels, conf)
results = pd.DataFrame({
    'pdu_id': pdu_ids,
    'cluster': labels,
    'confidence': conf
})
```

## CLI Backwards Compatibility

All existing CLI calls should still work via thin bin/ wrappers:

```bash
# Extract features (same interface as before)
python bin/extract_features.py --db per_aa_sqlite/pdus_L.sqlite --aa L --radius 15 --out-dir analysis/features

# Train & cluster (new unified interface)
python bin/train_model.py --features analysis/features/pdu_features_L.npz --latent-dim 16 --epochs 10

# Compress hypothesis test (now parameterizable)
python bin/analyze_compression.py --aa L --distances 10,15 --dims 8,16,32,64,128,256
```

## Refactoring Phases

### Phase 1: Library Structure Setup (Day 1)
- Create `pdusearch/` package with `__init__.py` files
- Create `config.py` with Config, TrainingConfig, ClusteringConfig classes
- Create `logging_utils.py` for centralized logging
- Create `io/database.py` wrapping sqlite3 operations (from export_pdu_features.py)
- Create `utils/types.py` with constants (AA_ORDER, SS_ORDER, etc from kmers/)

### Phase 2: Feature Extraction Module (Day 1-2)
- Create `features/encodings.py` with Encoding enum and strategies for:
  - AA_20 (current 20 amino acids)
  - FUNCTIONAL_5 (hydrophobic, polar, charged+, charged-, special)
  - AA_20_WITH_CONSERVATION (20 AA + Shannon entropy per position)
- Create `features/extractor.py` with FeatureExtractor class
- Create `features/alignment.py` with PyCoM integration
- **First test point**: `python -c "from pdusearch.features import FeatureExtractor; ..."`

### Phase 3: Models Module (Day 1-2)
- Create `models/autoencoder.py` with Autoencoder class (consolidate from train_pdu_autoencoder*.py)
- Create `models/training.py` with train_autoencoder() function and TrainingConfig
- Ensure CUDA/CPU detection works
- **Test point**: Verify models import and train on small dataset

### Phase 4: Clustering Module (Day 1-2)
- Create `clustering/hdbscan.py` with HDBSCANClusterer wrapper
- Create `clustering/metrics.py` with ClusterMetrics class
- **Test point**: Run compression hypothesis test using new library

### Phase 5: Enrichment & Validation (Day 2-3)
- Create `enrichment/pycom.py` with PyCoM query functions
- Create `validation/stability.py` and `validation/grid_search.py`
- Refactor from existing bin/validate*, bin/enrich* scripts

### Phase 6: CLI & Notebooks (Day 3)
- Consolidate 40 bin/*.py scripts into 8-10 focused CLI tools
- Create Jupyter notebooks demonstrating library usage
- Update README with examples

## Implementation Strategy

1. **Don't delete old scripts yet** — create new library in parallel
2. **Test each module independently** before merging phases
3. **Use existing tests** (if any) as regression checks
4. **Keep ARCHER2 jobs running** in background (Task #4) during refactoring
5. **Once refactoring complete** → implement Tasks #1-6 (feature engineering) in clean library

## Expected Benefits

✅ **Code reuse**: 40 scripts → 8-10 focused CLI tools + shared library  
✅ **Notebook compatibility**: Import and use classes directly, not just command-line  
✅ **Feature flexibility**: Add new encodings without touching existing scripts  
✅ **Testing**: Modular code is easier to unit test  
✅ **Maintenance**: Single source of truth for each concept (e.g., one FeatureExtractor, not 5)  
✅ **Future scalability**: GNN, alternative models, new clustering algorithms easier to add  

## Timeline & Effort

- **Refactoring**: ~2-3 days
- **Payoff**: Saves ~1-2 hours per feature iteration (x6 planned features = 6-12 hours saved)
- **ROI**: Positive within first 2-3 feature additions

---

**Start with Phase 1 & 2**: Get FeatureExtractor importable and tested, then expand.
