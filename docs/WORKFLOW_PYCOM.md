# PyCoM Direct Enrichment (Local Database)

**Problem:** REST API and UniProt mapping both failed to return results.

**Solution:** Use LOCAL PyCoM Python library with downloaded database files.

## Setup: Download PyCoM Database

1. Visit: https://pycom.brunel.ac.uk/database.html
2. Download both files:
   - `pycom.db` (main database)
   - `pycom.mat` (coevolution matrices)
3. Place them in a known directory (e.g., `/Volumes/local_drive/pycom/`)

## Workflow

### Step 1: Build PyCoM PDB Mapping (One-time, ~20 min)

```bash
python3 pdu-extractor/scripts/fetch_pycom_pdb_local.py \
    --db-path /path/to/pycom.db \
    --mat-path /path/to/pycom.mat \
    --out-file pdb_pycom_mapping.csv
```

**What it does:**
1. Initialize local PyCoM database (Python library)
2. Query for all entries with `has_pdb=true`
3. Use `DataLoader` to add biological features:
   - CATH classification
   - EC number (enzyme)
   - PDB IDs
   - Pfam domains
   - Disease associations
4. Save to CSV: `pdb_pycom_mapping.csv`

**Output columns:**
```
uniprot_id | pdb_id | cath_class | enzyme_commission | domain | disease_name | organism_id | sequence_length | helix_frac | turn_frac | strand_frac
```

**Size estimate:** 450K+ UniProt entries with PDB = ~100-200 MB CSV (reusable for all AAs)

### Step 2: Run Enrichment for Each AA

Once `pdb_pycom_mapping.csv` is created, run enrichment for each amino acid:

```bash
# LEU (example)
python3 pdu-extractor/scripts/analyze_clusters_pycom_direct.py \
    --clusters analysis/clusters_umap/pdu_clusters_L_umap.csv \
    --db pdb_db/pdus_L.sqlite \
    --pycom-mapping pdb_pycom_mapping.csv \
    --aa L \
    --out-dir analysis/pycom_enrichment

# All AAs (parallel batch)
for aa in A C D E F G H I K L M N P Q R S T V W Y; do
  python3 pdu-extractor/scripts/analyze_clusters_pycom_direct.py \
    --clusters "analysis/clusters_umap/pdu_clusters_${aa}_umap.csv" \
    --db "pdb_db/pdus_${aa}.sqlite" \
    --pycom-mapping pdb_pycom_mapping.csv \
    --aa "$aa" \
    --out-dir analysis/pycom_enrichment &
done
wait
```

**Mapping chain:**
```
PDU ID → [SQLite database] → PDB ID → [pdb_pycom_mapping.csv] → CATH/EC/GO
```

**Runtime:** ~1-2 min per AA (parallel: ~5 min total for all 20)

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
