# PyCoM-Based Cluster Enrichment Workflow

## Overview
Use Brunel's PyCoM API to perform fast, reliable cluster enrichment analysis instead of hitting external APIs.

**Advantages:**
- ✅ Local Brunel resource (no external API rate limits or timeouts)
- ✅ Precomputed metadata (instant queries)
- ✅ Comprehensive annotation (CATH, EC, GO, Pfam, ligands, cofactors, PTM, disease)
- ✅ Much faster than external APIs

## Two-Step Workflow

### Step 1: Build Unified PDB → PyCoM Metadata Reference (Run Once)

```bash
cd /Users/sarathdantu/Projects/vsprojects/mdpdu/pdu-extractor

python3 scripts/build_pdb_metadata_via_pycom.py \
  --db-dir /Users/sarathdantu/Library/CloudStorage/OneDrive-BrunelUniversityLondon/Work/Research/Manuscript/Me/3D\ Sequence\ Evolution/data/PDUs/pdu_db \
  --out-file pdb_metadata_pycom.json \
  --pdb-uniprot-cache .pdb_uniprot_map.json
```

**What it does:**
1. Extracts unique PDB IDs across all 20 amino acid databases
2. Maps PDB → UniProt (one-time fetch from PDBe API)
3. Queries PyCoM with UniProt IDs → gets CATH, EC, GO, Pfam, ligand, cofactor, PTM, disease
4. Caches locally → reusable for all downstream analysis

**Time:** ~5-10 minutes (mostly waiting for PyCoM API responses)
**Output:** `pdb_metadata_pycom.json` (~2-5 MB)

### Step 2: Analyze Cluster Enrichment (Instant, Run for Each AA)

Once the reference is built, analyze any amino acid instantly:

```bash
# Single AA
python3 scripts/analyze_cluster_enrichment_fast.py \
  --clusters analysis/clusters_umap/pdu_clusters_L_umap.csv \
  --db /Users/sarathdantu/Library/CloudStorage/OneDrive-BrunelUniversityLondon/Work/Research/Manuscript/Me/3D\ Sequence\ Evolution/data/PDUs/pdu_db/pdus_L.sqlite \
  --metadata pdb_metadata_pycom.json \
  --aa L
```

Or batch (all AAs):

```bash
# Batch analysis
for aa in L A E G D V I K R P N S T Q F W Y H M C; do
  echo "Processing AA=$aa..."
  python3 scripts/analyze_cluster_enrichment_fast.py \
    --clusters "analysis/clusters_umap/pdu_clusters_${aa}_umap.csv" \
    --db "/Users/sarathdantu/Library/CloudStorage/OneDrive-BrunelUniversityLondon/Work/Research/Manuscript/Me/3D Sequence Evolution/data/PDUs/pdu_db/pdus_${aa}.sqlite" \
    --metadata pdb_metadata_pycom.json \
    --aa "$aa"
done
```

**Time:** ~30 seconds per AA (no API calls, pure local JSON lookups)
**Output per AA:**
- `contingency_cath_{AA}.csv`: Cluster × CATH superfamily enrichment
- `contingency_ec_{AA}.csv`: Cluster × EC number enrichment
- `clusters_enriched_{AA}.csv`: Full cluster data with all PyCoM annotations

## Outputs Explained

### For LEU (740 clusters):

**`contingency_cath_L.csv`**
- Rows: CATH superfamilies
- Cols: Cluster IDs (0-740)
- Values: PDU counts per (CATH × Cluster) pair
- **Interpretation:** If clusters are CATH-specific → real biological signal

**`contingency_ec_L.csv`**
- Rows: EC numbers
- Cols: Cluster IDs
- **Interpretation:** If certain clusters enrich for specific enzyme types → functional grouping

**`clusters_enriched_L.csv`**
Full data with columns:
- `pdu_id`: PDU ID
- `cluster`: Cluster assignment
- `pdb_id`: PDB ID
- `cath`: CATH superfamily (from PyCoM)
- `enzyme_ec`: EC number (from PyCoM)
- `biological_process`: GO biological process terms
- `molecular_function`: GO molecular function terms
- `cellular_component`: GO cellular component terms
- `domain`: Protein domain families
- `ligand`: Associated ligands
- `cofactor`: Associated cofactors
- `ptm`: Post-translational modifications
- `disease`: Associated diseases

## Interpretation Questions

**For LEU 740 clusters:**

1. **Are clusters meaningful or artifacts?**
   - If `contingency_cath_L.csv` shows cluster-specific CATH families → REAL SIGNAL
   - If all clusters have mixed CATH families → artifacts

2. **Why 740 clusters vs E's 3 clusters?**
   - LEU is abundant (2.1M PDUs) → more diverse structural contexts
   - E is rare (15.6K PDUs) → fewer distinct contexts
   - OR: LEU genuinely has more diverse neighborhoods (richer structure)

3. **Do clusters have distinct "recipes"?**
   - Compare CATH, EC, GO, domain distributions across clusters
   - If cluster 738 is 90% helix-bundle proteins and cluster 687 is 80% β-sandwich → DIFFERENT RECIPES

## Next Steps After Analysis

1. **Visualization:**
   - Heatmap: Cluster × CATH family (show separation)
   - Network: PDB structures colored by cluster assignment
   - Distribution plots: CATH/EC/GO per cluster

2. **Cluster Profiles:**
   - For each cluster, extract "design recipe" (avg CATH, EC, GO, domain composition)
   - Compare 740 LEU recipes: how many distinct patterns?

3. **Cross-AA Comparison:**
   - Same analysis for A, E, G
   - Why such different clustering counts?
   - Are universal principles or AA-specific?

## Files in This Workflow

- `build_pdb_metadata_via_pycom.py`: Build PyCoM metadata reference (run once)
- `analyze_cluster_enrichment_fast.py`: Fast enrichment analysis (run per AA)
- `pdb_metadata_pycom.json`: Reusable reference (commit to git)
- `PYCOM_WORKFLOW.md`: This file
