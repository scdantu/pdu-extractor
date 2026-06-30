# Compression Hypothesis Test: New Library-Based Version

## Overview

Unified compression test using new **pdusearch library** to compare three encoding strategies and determine which factors drive the 88% noise observed in LEU clustering.

**Test Matrix**:
- 3 encodings (900D, 225D, 915D)
- 2 distance cutoffs (10Å, 15Å)  
- 6 latent dimensions (8, 16, 32, 64, 128, 256)
- 3 amino acids (G, A, L)
- **Total: 108 test configurations**

---

## Quick Start

### Local Testing (Small Dataset)
```bash
cd pdu-extractor

# Test single encoding with one dimension
python bin/compression_test_unified.py \
    --aa L \
    --encodings functional_5 \
    --dims 16 \
    --epochs 2 \
    --log-level DEBUG

# Output: analysis/compression_test/summary_L_functional_5.csv
```

### ARCHER2 Full Test (All 108 Configurations)
```bash
cd pdu-extractor

# Submit all tests
sbatch submit/submit_compression_test_unified.slurm

# Monitor
squeue -j <jobid>
tail -f logs/compress_unified_*.out

# After completion (~12-20 hours)
python bin/analyze_compression_results.py --aa L
python bin/analyze_compression_results.py --aa A
python bin/analyze_compression_results.py --aa G
```

---

## What's Being Tested

### Hypothesis 1: Distance Cutoff Diversity
**Question**: Does 15Å bring in too much neighbor diversity vs 10Å?

**Prediction**: 
- If TRUE: 10Å clusters better (lower noise)
- If FALSE: Distance doesn't matter; problem is elsewhere

### Hypothesis 2: Compression Aggressiveness
**Question**: Is 900D→16D compression too aggressive?

**Prediction**:
- If TRUE: Higher latent dims (32D, 64D) cluster better
- If FALSE: 16D is adequate; problem is upstream

### Hypothesis 3 (NEW): Feature Granularity
**Question**: Does 20-amino-acid identity matter vs functional grouping?

**Prediction**:
- If TRUE: Functional-5 (225D) clusters as well as AA-20 (900D)
- If FALSE: Exact AA identity is important

### Hypothesis 4 (NEW): Conservation Signal
**Question**: Does evolutionary conservation improve clustering?

**Prediction**:
- If TRUE: AA-20-Conservation (915D) better than AA-20
- If FALSE: Sequence conservation adds noise

---

## Output Format

### Results CSV
**File**: `analysis/compression_test/summary_{AA}_{ENCODING}.csv`

**Columns**:
```
encoding              (aa_20, functional_5, aa_20_conservation)
distance              (10, 15)
latent_dim            (8, 16, 32, 64, 128, 256)
input_features        (900, 225, or 915)
n_pdus                (PDU count for this AA)
n_clusters            (clusters found by HDBSCAN)
noise_percent         (% unassigned points)
avg_confidence        (average cluster membership prob)
n_assigned            (# points in clusters)
training_loss         (final autoencoder loss)
```

### Example Results
```csv
encoding,distance,latent_dim,input_features,n_pdus,n_clusters,noise_percent,avg_confidence,n_assigned,training_loss
aa_20,10,16,900,2100000,142,87.3,0.642,267000,0.0892
functional_5,10,16,225,2100000,189,84.2,0.671,331000,0.0567
aa_20_conservation,10,16,915,2100000,151,86.8,0.655,278000,0.0945
...
```

---

## Analysis Script

### Usage
```bash
python bin/analyze_compression_results.py --aa L --out-dir analysis/compression_test
```

### Output Sections

**1. By Encoding Analysis**
```
AA_20:
  Feature dimension: 900
  Best config:
    Distance:  15Å
    LatentDim: 16D
    Noise:     85.2%
    Clusters:  156
    Confidence: 0.658
```

**2. Distance Effect**
```
ANALYSIS: DISTANCE CUTOFF EFFECT

aa_20:
  10Å: noise=87.1% (min=82.3%, max=91.4%)
  15Å: noise=85.2% (min=80.1%, max=89.3%)
  → 15Å is 2.2% better
```

**3. Compression Effect**
```
ANALYSIS: COMPRESSION EFFECT

aa_20:
    8D: noise=92.4% (conf=0.512)
   16D: noise=85.2% (conf=0.658)
   32D: noise=83.1% (conf=0.681)
   64D: noise=82.5% (conf=0.689)
  128D: noise=82.3% (conf=0.691)
  256D: noise=82.2% (conf=0.692)
  → Elbow at 32D (diminishing returns after 16D)
```

**4. Encoding Comparison**
```
ENCODING COMPARISON

Encoding                  Noise %    Config             Conf    
aa_20                     82.3       15Å, 256D          0.692   
functional_5              79.8       15Å, 256D          0.704   
aa_20_conservation        81.5       15Å, 256D          0.699   

KEY FINDINGS:
  • Best encoding: functional_5 (79.8% noise)
  • Difference: 2.5% (encoding choice matters)
  • Functional-5 vs AA-20: -2.5% (better)
  • Conservation signal: -0.8% (slightly helpful)
```

---

## Expected Results

### Timeline
| AA | PDUs | Time | Notes |
|---|---|---|---|
| G | 19K | ~30 min | Fast (small) |
| A | 1.9M | ~2 hrs | Medium |
| L | 2.1M | ~3 hrs | Largest |
| **Total** | - | **~5-6 hrs** | Can run in parallel |

### Noise Expectations

**Current (AA-20, 15Å, 16D)**: ~88% noise

**Hypothesized improvements**:
- Different distance: -2% to +3% (small effect)
- Higher latent dim: -5% to -10% (moderate effect)
- Functional encoding: -5% to 0% (likely neutral/slightly better)
- Conservation: -2% to +2% (uncertain)

**Realistically**: Probably still 75-85% noise (88% is largely biological truth, not artifact)

---

## Interpreting Results

### If Noise Doesn't Change Much
**Conclusion**: 88% noise reflects genuine biological diversity. Not a modeling artifact.

**Implication**: Every LEU neighborhood is mostly unique. Only ~12% form coherent groups.

**Next Step**: Accept noise level; focus on cluster quality (avg_confidence) instead of count.

### If Functional-5 < AA-20 Noise
**Conclusion**: Exact amino acid identity doesn't matter; functional properties do.

**Implication**: Can use simpler 225D encoding. 4x compression with no quality loss.

**Next Step**: Use FUNCTIONAL_5 for all downstream analysis.

### If Higher Latent Dim Helps
**Conclusion**: 16D is indeed too compressed.

**Implication**: Trade speed for quality; use 32D or 64D latent space.

**Next Step**: Retrain all models with optimal dimension.

### If Conservation Helps
**Conclusion**: Evolutionary information matters. Invest in PyCoM integration.

**Implication**: Clusters separate by biological function, not just structure.

**Next Step**: Complete PyCoM integration; use conservation in all encodings.

---

## Running on ARCHER2

### Step 1: Prepare Environment
```bash
# On ARCHER2 login node
cd /work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor

# Verify setup
source /work/e280/e280-Dantu/scdantu/Projects/pdus/pypdu/bin/activate
python -c "from pdusearch.features import FeatureExtractor; print('✓ Ready')"
```

### Step 2: Submit Test
```bash
sbatch submit/submit_compression_test_unified.slurm
# → Returns: Submitted batch job XXXXXXX

# Check status
squeue -j XXXXXXX
watch squeue -j XXXXXXX  # Real-time updates

# Monitor output
tail -f logs/compress_unified_XXXXXXX.out
```

### Step 3: Analyze Results
```bash
# After job completes (~12 hours)
python bin/analyze_compression_results.py --aa L
python bin/analyze_compression_results.py --aa A
python bin/analyze_compression_results.py --aa G

# Results
ls -lh analysis/compression_test/summary_*.csv
```

### Step 4: Retrieve Results
```bash
# On local machine
scp archer2:/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor/analysis/compression_test/*.csv .

# Analyze locally
python bin/analyze_compression_results.py --aa L
```

---

## Troubleshooting

### Import Errors
```
ModuleNotFoundError: No module named 'pdusearch'
```
**Fix**: Check PYTHONPATH includes pdu-extractor root
```bash
export PYTHONPATH="/work/e280/..../pdu-extractor:$PYTHONPATH"
```

### CUDA/Device Issues
```
RuntimeError: CUDA out of memory
```
**Fix**: Already handled by auto device detection (falls back to CPU)

### Database Not Found
```
FileNotFoundError: Database not found: per_aa_sqlite/pdus_L.sqlite
```
**Fix**: Ensure running from pdu-extractor root: `cd /work/.../pdu-extractor`

### Out of Time
Increase `--time` in SLURM script:
```bash
#SBATCH --time=20:00:00  # Was 12:00:00
```

---

## Code Structure

### Main Test Script
**File**: `bin/compression_test_unified.py`

**Key Classes**:
- `CompressionTestRunner` — Orchestrates all tests
- Methods:
  - `run()` — Execute full test matrix
  - `_extract_features()` — Load data
  - `_test_combination()` — Single test
  - `_save_results()` — Write CSV

**Usage**:
```python
runner = CompressionTestRunner(
    aa="L",
    encodings=["aa_20", "functional_5", "aa_20_conservation"],
    distances=[10, 15],
    latent_dims=[8, 16, 32, 64, 128, 256],
    epochs=10
)
results = runner.run()
```

### Analysis Script
**File**: `bin/analyze_compression_results.py`

**Functions**:
- `load_results()` — Read CSV files
- `analyze_by_encoding()` — Per-encoding stats
- `analyze_distance_effect()` — 10Å vs 15Å
- `analyze_compression_effect()` — Latent dim impact
- `analyze_encoding_comparison()` — Head-to-head

---

## Success Criteria

**Test passes if**:
- ✓ All 108 configurations complete
- ✓ Each produces valid noise% and confidence scores
- ✓ Results saved to CSV with no NaNs
- ✓ Analysis script runs without errors
- ✓ Can identify best encoding + distance + dimension combination

**Quality check**:
- ✓ Noise% in range 70-95% (biologically reasonable)
- ✓ Avg confidence 0.5-0.8 (cluster quality is moderate)
- ✓ Consistent pattern (higher latent dims → lower noise)

---

## Next Steps

### If Test Succeeds
1. Identify optimal configuration (encoding + distance + dimension)
2. Retrain all 20 AAs with optimal settings
3. Run enrichment validation (CATH/EC separation)
4. Create notebooks + CLI scripts demonstrating results

### If Test Fails
1. Debug individual components (feature extraction, training, clustering)
2. Run single test with verbose logging
3. Verify database availability and paths
4. Check PyTorch/HDBSCAN versions

---

## Reference

**Compression Test Matrix Summary**:
```
Encodings:      AA_20 (900D)
                FUNCTIONAL_5 (225D)  ← Hypothesis 1
                AA_20_WITH_CONSERVATION (915D)  ← Hypothesis 4

Distances:      10Å  ← Hypothesis 2 (feature diversity)
                15Å

Latent Dims:    8, 16, 32, 64, 128, 256  ← Hypothesis 3 (compression)

AAs:            G (19K), A (1.9M), L (2.1M)  ← Representative sizes

Metrics:        noise_percent, avg_confidence, n_clusters, training_loss
```

**Estimated Output Size**: ~100 KB CSV files (one per encoding per AA)

**Estimated Storage**: ~1 MB total results

---

**Document Version**: 1.0
**Date**: June 29, 2026
**Status**: Ready for deployment
