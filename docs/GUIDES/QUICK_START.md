# Quick Start Guide

## Installation

### Local Setup
```bash
cd pdu-extractor
export PYTHONPATH="$(pwd):$PYTHONPATH"
python -c "from pdusearch.features import FeatureExtractor; print('✓ Ready')"
```

### ARCHER2 Setup
```bash
cd /work/e280/.../pdu-extractor
source /work/.../pypdu/bin/activate
export PYTHONPATH="/work/.../pdu-extractor:$PYTHONPATH"
```

## Run Your First Test

### Quick Local Test (2 minutes)
```bash
python bin/compression_test_unified.py \
    --aa L \
    --encodings functional_5 \
    --dims 16 \
    --epochs 2 \
    --log-level INFO
```

**Output**: `analysis/compression_test/summary_L_functional_5.csv`

### Full ARCHER2 Test (12-20 hours)
```bash
sbatch submit/submit_compression_test_unified.slurm
```

Monitor:
```bash
squeue -j <JOBID>
tail -f logs/compress_unified_*.out
```

Analyze results:
```bash
python bin/analyze_compression_results.py --aa L
python bin/analyze_compression_results.py --aa A
python bin/analyze_compression_results.py --aa G
```

## Basic Usage

### Extract Features
```python
from pdusearch.features import FeatureExtractor, Encoding
from pdusearch.config import Config

config = Config(db="per_aa_sqlite/pdus_L.sqlite", radius=15.0)

# Extract with standard encoding (900D)
extractor = FeatureExtractor(config, encoding=Encoding.AA_20)
X, pdu_ids = extractor.extract_for_aa("L")

# Extract with functional encoding (225D, faster)
extractor_func = FeatureExtractor(config, encoding=Encoding.FUNCTIONAL_5)
X_func, ids_func = extractor_func.extract_for_aa("L")

# Extract with conservation (915D, needs PyCoM)
extractor_cons = FeatureExtractor(config, encoding=Encoding.AA_20_WITH_CONSERVATION)
X_cons, ids_cons = extractor_cons.extract_for_aa("L")
```

### Train Autoencoder
```python
from pdusearch.models import train_autoencoder, TrainingConfig
import numpy as np

# Train autoencoder
config = TrainingConfig(
    latent_dim=16,
    epochs=10,
    batch_size=32,
    device='auto'  # Auto-detects GPU/CPU
)

model, Z, losses = train_autoencoder(X, config)
# Z.shape = (n_samples, 16)
```

### Cluster Embeddings
```python
from pdusearch.clustering import HDBSCANClusterer, ClusterMetrics

# Cluster
clusterer = HDBSCANClusterer(min_cluster_size=200, min_samples=5)
labels, confidences = clusterer.cluster(Z)

# Evaluate
metrics = ClusterMetrics(labels, confidences)
print(f"Clusters: {metrics.n_clusters}")
print(f"Noise: {metrics.noise_percent:.1f}%")
print(f"Quality: {metrics.avg_confidence:.3f}")
```

## File Locations

| Item | Location |
|------|----------|
| Library code | `pdusearch/` |
| CLI scripts | `bin/` |
| Documentation | `docs/` |
| Data output | `analysis/` |
| ARCHER2 batch | `submit/` |

## Documentation

- **[Full Testing Guide](../TESTING/COMPRESSION_TEST_GUIDE.md)** — Detailed test documentation
- **[Architecture](../IMPLEMENTATION/REFACTORING_COMPLETE.md)** — Library design
- **[API Docs](../../pdusearch/)** — Class and function reference

## Troubleshooting

### Import Error
```
ModuleNotFoundError: No module named 'pdusearch'
```
**Fix**: Set PYTHONPATH
```bash
export PYTHONPATH="/path/to/pdu-extractor:$PYTHONPATH"
```

### Database Not Found
```
FileNotFoundError: Database not found
```
**Fix**: Run from pdu-extractor root directory

### Out of Memory
**Fix**: Already handled - autom atically falls back to CPU

## Next Steps

1. **Run quick test** → Verify everything works locally
2. **Submit ARCHER2 job** → Get real results
3. **Analyze results** → Identify optimal configuration
4. **If successful** → Create CLI scripts + notebooks

See [Testing Guide](../TESTING/COMPRESSION_TEST_GUIDE.md) for full details.
