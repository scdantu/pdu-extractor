# PyCoM Direct Enrichment (No UniProt Mapping)

**Problem:** PDB → UniProt mapping yielded 0 results. UniProt lookup too brittle.

**Solution:** Query PyCoM API directly for PDB entries, build local reference, skip UniProt entirely.

## Workflow

### Step 1: Fetch PyCoM PDB Mapping (One-time, ~10 min)

```bash
python3 pdu-extractor/scripts/fetch_pycom_pdb_mapping.py \
    --out-file pdb_pycom_mapping.csv \
    --batch-size 100
```

**What it does:**
1. Query PyCoM API with `has_pdb=true`
2. Extract for each PDB entry:
   - PDB ID
   - CATH classification
   - EC number (enzyme)
   - GO terms (biological process, molecular function, cellular component)
   - Pfam domains
   - Disease associations
   - Organism, length

3. Save to CSV: `pdb_pycom_mapping.csv`

**Output columns:**
```
pdb_id | cath | enzyme_ec | biological_process | molecular_function | cellular_component | pfam | disease | organism | length
```

**Size estimate:** 200K+ PDB entries = ~100 MB CSV (reusable for all AAs)

### Step 2: Run Enrichment for Each AA

```bash
# LEU
python3 pdu-extractor/scripts/analyze_clusters_pycom_direct.py \
    --clusters analysis/clusters_umap/pdu_clusters_L_umap.csv \
    --db pdus_L.sqlite \
    --pycom-mapping pdb_pycom_mapping.csv \
    --aa L \
    --out-dir analysis/pycom_enrichment

# All AAs (batch)
for aa in A C D E F G H I K L M N P Q R S T V W Y; do
  python3 pdu-extractor/scripts/analyze_clusters_pycom_direct.py \
    --clusters "analysis/clusters_umap/pdu_clusters_${aa}_umap.csv" \
    --db "pdus_${aa}.sqlite" \
    --pycom-mapping pdb_pycom_mapping.csv \
    --aa "$aa" \
    --out-dir analysis/pycom_enrichment
done
```

**Mapping chain:**
```
PDU ID → [database] → PDB ID → [PyCoM CSV] → CATH/EC/GO
```

**Runtime:** ~2 min per AA

### Step 3: Check Results

```bash
# LEU enrichment
ls -lh analysis/pycom_enrichment/

# View enrichment
head contingency_cath_L.csv
head contingency_ec_L.csv
```

## Expected Outputs

```
analysis/pycom_enrichment/
├── pdb_pycom_mapping.csv          ← One-time reference
├── contingency_cath_L.csv         ← Cluster × CATH family
├── contingency_ec_L.csv           ← Cluster × EC number
├── clusters_enriched_L.csv        ← Full cluster data + metadata
├── contingency_cath_A.csv
├── contingency_ec_A.csv
├── clusters_enriched_A.csv
└── ... (repeat for all AAs)
```

## Key Advantage

✅ **Direct API querying** — no intermediate UniProt mapping
✅ **Single reference file** — reusable across all 20 AAs
✅ **Complete metadata** — CATH, EC, GO, Pfam, disease all in one place
✅ **Reproducible** — same CSV source for all analyses

## Interpretation

**For LEU clusters:**

**CATH enrichment example:**
```
χ² = 2,847, p = 1.2e-45 ← HIGHLY SIGNIFICANT

Contingency: Cluster × CATH family
  If cluster_738 is 87% "all-α fold"
  While cluster_687 is 63% "all-α fold"
  → Clusters do separate by structural family
```

**EC enrichment example:**
```
χ² = 1,203, p = 3.4e-22 ← SIGNIFICANT

  If cluster_123 enriches for "1.1.1.*" (oxidoreductases)
  While cluster_456 is mostly "unannotated"
  → Cluster structure correlates with enzymatic function
```

**GO enrichment example:**
```
If cluster_687 PDUs come from proteins annotating to:
  - "ATP binding"
  - "Protein kinase activity"
  
While cluster_738 PDUs annotate to:
  - "Hydrolase activity"
  - "Catalytic activity"
  
→ Clusters are functionally distinct
```

## If Results Show Signal

✅ CATH/EC/GO separate by cluster → Clusters have biological meaning
✅ χ² p-value < 0.05 → Statistically significant
✅ Large effect size → Worth publishing

→ **Trigger full validation pipeline** (hyperparameter grid + bootstrap)

## If Results Show No Signal

❌ No CATH/EC/GO separation → Clusters are noise
❌ p-value > 0.05 → Not significant
❌ Small effect size → Hard to publish

→ **Need to rethink approach:**
- Lower min_cluster_size? (too many noise points now)
- Different features? (current 900-dim not sufficient)
- Different embedding? (16-dim latent space not separable)
