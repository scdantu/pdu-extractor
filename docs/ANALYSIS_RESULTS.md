# Cluster Analysis: All 20 Amino Acids

**Status:** 18/20 AAs clustered (W, C still running)

## Key Findings

### 1. NOISE LEVELS ARE VERY HIGH (88%+)

| AA | PDUs | Clusters | Noise | Noise % |
|----|----|----------|-------|---------|
| **LEU** | 2,086,440 | 538 | 1,848,289 | **88.6%** |
| **ALA** | 1,913,415 | 482 | 1,708,183 | 89.3% |
| **GLY** | 19,435 | 3 | 305 | **1.6%** ← OUTLIER |
| **Average** | - | - | - | **88.3%** |

### 2. CLUSTER DIVERSITY VARIES DRAMATICALLY

```
LEU:    538 clusters / 2.1M PDUs = 1 cluster per 3,873 PDUs
ALA:    482 clusters / 1.9M PDUs = 1 cluster per 3,968 PDUs
GLY:      3 clusters / 19.4K PDUs = 1 cluster per 6,478 PDUs ← Very uniform
PRO:    246 clusters / 1.1M PDUs = 1 cluster per 4,468 PDUs
```

**179.3x difference** between LEU (538 clusters) and GLY (3 clusters)

### 3. CLUSTER QUALITY IS HIGH (for assigned PDUs)

- Mean confidence: 0.887 (good)
- Only GLY has >0.95 confidence (0.987)
- Others range 0.849-0.920

### 4. WHY SO MUCH NOISE?

**Hypothesis:** min_cluster_size=200 is TOO HIGH

Current settings:
- L/A: min_cluster_size=200
- Most others: 50-100

If we use clustering where each cluster has ≥200 PDUs:
- LEU: 2,086,440 / 200 = max ~10,432 clusters (we got 538 + 88% noise)
- This suggests 538 clusters × ~3,873 PDUs/cluster ≈ 2.1M (matches!)

**The noise represents PDUs that don't form cohesive groups of ≥200 at this min_cluster_size.**

## What This Means

✅ **Good news:**
- Clustering is working (high confidence for assigned PDUs)
- 538 LEU clusters is reasonable (~3900 PDUs per cluster)
- Diversity is expected (different structural contexts)

⚠️ **Needs investigation:**
- Is 88% noise acceptable for publication?
- Should we lower min_cluster_size to capture more PDUs?
- Or is this legitimate "no clear neighbors" representation?

## Action Items

### 1. IMMEDIATE: Hyperparameter Grid Search

**This is exactly why we built the validation pipeline!**

Run parallel grid search:
```bash
sbatch -N 2 --array=0-19 submit_grid_search_parallel.slurm
```

**Test different min_cluster_size values:**
- Current: 200 → 88% noise
- Try: 50, 100, 150 → see if noise decreases while maintaining separation

### 2. BEFORE PyCoM ENRICHMENT

Determine: Is 88% noise okay for publication?

**Consider:**
- Rosetta/FoldX typically use all-vs-all interactions (no "noise" concept)
- PDU clusters are MORE conservative (require local coherence)
- Noise might represent biologically valid "non-canonical" neighborhoods

### 3. PINNED TASK FOR FUTURE

After PyCoM validation shows signal:
```
╔════════════════════════════════════════════════════════════╗
║                                                            ║
║  PARALLEL RIGOROUS VALIDATION PIPELINE (PINNED)            ║
║                                                            ║
║  Phase 1: Feature Ablation                                ║
║  Phase 2: Hyperparameter Grid (parallel array)            ║
║  Phase 3: Stability Validation (parallel array)           ║
║  Phase 4: Aggregation & Report                            ║
║                                                            ║
║  Timeline: ~7-8 hours per AA (24h wall time OK)           ║
║                                                            ║
║  Location: pdu-extractor/                                 ║
║    - submit_grid_search_parallel.slurm                    ║
║    - submit_stability_parallel.slurm                      ║
║    - PARALLEL_VALIDATION_WORKFLOW.md                      ║
║                                                            ║
╚════════════════════════════════════════════════════════════╝
```

## Statistics Summary

```
Total PDUs:        21,067,067 (14 AA = 2.1B instances if full 20)
Total Clusters:    5,580
Total Noise:       18,596,597 (88.3%)

Cluster Size Range:
  Min: 3 clusters (GLY - very uniform neighborhoods)
  Max: 538 clusters (LEU - very diverse neighborhoods)
  Avg: 310 clusters

Quality Range:
  Min confidence: 0.849 (V)
  Max confidence: 0.987 (G)
  Avg confidence: 0.887
```

## Next Step

**NOW:** Run hyperparameter grid search on ARCHER2
- 20 grid tasks in parallel → 2-3 hours elapsed
- Find optimal min_cluster_size for each AA
- Decision: Keep current settings or adjust?

**THEN:** If signal in PyCoM enrichment → Full validation pipeline
