# Cluster Enrichment Analysis Workflow

## Problem
Fetching CATH/EC/GO metadata for each amino acid separately is inefficient — PDB IDs overlap across all 20 AAs.

## Solution
Two-step workflow:

### Step 1: Build Unified PDB Metadata Reference (Run Once)
Extract unique PDB IDs across all 20 amino acid databases, fetch metadata once, save to a single reference file.

```bash
cd /Users/sarathdantu/Projects/vsprojects/mdpdu/pdu-extractor

python3 scripts/build_pdb_metadata_reference.py \
  --db-dir /Users/sarathdantu/Library/CloudStorage/OneDrive-BrunelUniversityLondon/Work/Research/Manuscript/Me/3D\ Sequence\ Evolution/data/PDUs/pdu_db \
  --out-file pdb_metadata_reference.json \
  --api
```

**What it does:**
1. Queries all 20 per-AA databases
2. Extracts unique PDB IDs (~1000-5000 structures)
3. Fetches from PDBe API:
   - **CATH**: Structural classification (superfamily)
   - **EC**: Enzyme classification
   - **GO**: Gene Ontology (biological process, molecular function, cellular component)
   - **Pfam**: Protein domain families
   - **UniProt**: Accession IDs
4. Saves to JSON with checkpoints (resume-friendly)

**Time:** ~30-60 minutes (includes API rate limiting)
**Output:** `pdb_metadata_reference.json` (~5-10 MB)

### Step 2: Run Fast Enrichment for Each AA
Once reference exists, analyze any amino acid instantly (no API calls).

```bash
# Single AA
python3 scripts/analyze_cluster_enrichment_fast.py \
  --clusters analysis/clusters_umap/pdu_clusters_L_umap.csv \
  --db /path/to/pdus_L.sqlite \
  --metadata pdb_metadata_reference.json \
  --aa L

# Batch (all AAs)
for aa in L A E G D V I K R P N S T Q F W Y H M C; do
  python3 scripts/analyze_cluster_enrichment_fast.py \
    --clusters "analysis/clusters_umap/pdu_clusters_${aa}_umap.csv" \
    --db "/path/to/pdus_${aa}.sqlite" \
    --metadata pdb_metadata_reference.json \
    --aa $aa
done
```

**What it does:**
1. Loads cluster assignments for AA
2. Maps PDU IDs → PDB IDs (from database)
3. Looks up metadata in reference file (instant)
4. Analyzes enrichment:
   - **CATH**: Do clusters separate by structural family?
   - **EC**: Do clusters correspond to enzyme types?
   - **GO**: Do clusters enrich for specific biological processes?
   - **Pfam**: Do clusters contain specific protein domains?
5. Outputs:
   - `contingency_cath_{aa}.csv`: Cluster × CATH family crosstab
   - `contingency_ec_{aa}.csv`: Cluster × EC number crosstab
   - `clusters_enriched_{aa}.csv`: Full cluster data with metadata

**Time:** ~1-2 minutes per AA
**Output:** 3 CSVs per AA

## Interpretation

### For LEU (740 clusters):
- **CATH enrichment**: If clusters cluster by structural family → real biological signal
- **EC enrichment**: If certain clusters are enzyme-rich → functional grouping
- **GO enrichment**: If clusters correlate with specific processes → biological meaning

### Cross-AA comparison:
- **Why does L have 740 clusters but E has 3?**
  - L is large (2.1M PDUs) → more diverse contexts
  - E is small (15.6K PDUs) → few distinct contexts
  - OR: L contexts are genuinely more diverse (structurally richer)

## Next: Visualization
Once enrichment is complete:
- Heatmap: Cluster × CATH family (show if clusters are CATH-specific)
- Network: PDB structures colored by cluster assignment
- Table: Per-cluster CATH/EC/GO profiles (design recipes)

## Files Created
- `build_pdb_metadata_reference.py`: Build unified metadata
- `analyze_cluster_enrichment_fast.py`: Fast enrichment analysis
- `pdb_metadata_reference.json`: Reusable reference (checked into git)
