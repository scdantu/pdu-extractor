# PDU Search: Structure-Driven Rational Protein Design

Complete modular framework for extracting, analyzing, and clustering protein dynamic units (PDUs) from protein structures to enable rational protein design.

## Quick Start

### Local Test
```bash
# Test feature extraction + training + clustering
python bin/compression_test_unified.py --aa L --encodings functional_5 --dims 16 --epochs 2
```

### ARCHER2 Full Test
```bash
sbatch submit/submit_compression_test_unified.slurm
python bin/analyze_compression_results.py --aa L
```

### Library Usage
```python
from pdusearch.features import FeatureExtractor, Encoding
from pdusearch.models import train_autoencoder, TrainingConfig
from pdusearch.clustering import HDBSCANClusterer, ClusterMetrics

# Extract features
X, ids = FeatureExtractor(config, encoding=Encoding.FUNCTIONAL_5).extract_for_aa("L")

# Train model
model, Z, losses = train_autoencoder(X, TrainingConfig(latent_dim=16, epochs=10))

# Cluster & evaluate
labels, conf = HDBSCANClusterer(min_cluster_size=200).cluster(Z)
metrics = ClusterMetrics(labels, conf)
print(f"Clusters: {metrics.n_clusters}, Noise: {metrics.noise_percent:.1f}%")
```

## Documentation

See [docs/](docs/) directory for complete documentation:
- **[Architecture Overview](docs/IMPLEMENTATION/REFACTORING_COMPLETE.md)** — Library design
- **[Testing Guide](docs/TESTING/COMPRESSION_TEST_GUIDE.md)** — Run hypothesis test
- **[Implementation Phases](docs/IMPLEMENTATION/)** — Step-by-step details
- **[API Reference](pdusearch/)** — Class & function docs

## Project Structure

```
pdu-extractor/
├── pdusearch/              # Main library (4,057 lines)
├── bin/                    # CLI scripts
├── submit/                 # ARCHER2 batch scripts
├── docs/                   # Documentation
└── notebooks/              # Jupyter examples (coming)
```

## Status

✓ **Library Complete**: Production-ready  
✓ **Tests Ready**: Deployment-ready  
→ **Next**: Execute on ARCHER2

[See docs/ for full documentation →](docs/)
