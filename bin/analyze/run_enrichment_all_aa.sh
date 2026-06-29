#!/bin/bash
# Run PyCoM enrichment for all 20 amino acids

set -e

cd "$(dirname "$0")/../../"

PYCOM_MAPPING="pycom/pycom_pdb_mapping_fixed.csv"
OUT_DIR="analysis/enrichment"

# Check if mapping exists
if [ ! -f "$PYCOM_MAPPING" ]; then
    echo "Error: PyCoM mapping not found: $PYCOM_MAPPING"
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "========================================================================"
echo "PyCoM Enrichment Analysis: All 20 Amino Acids"
echo "========================================================================"
echo ""

AAS=(A C D E F G H I K L M N P Q R S T V W Y)

for aa in "${AAS[@]}"; do
    CLUSTERS="analysis/clusters_umap/pdu_clusters_${aa}_umap.csv"
    DB="pdb_db/pdus_${aa}.sqlite"

    if [ ! -f "$CLUSTERS" ]; then
        echo "⊘ Skipping $aa: clusters not found ($CLUSTERS)"
        continue
    fi

    if [ ! -f "$DB" ]; then
        echo "⊘ Skipping $aa: database not found ($DB)"
        continue
    fi

    echo "▶ Running enrichment for $aa..."
    python3 bin/analyze/analyze_clusters_pycom_direct.py \
        --clusters "$CLUSTERS" \
        --db "$DB" \
        --pycom-mapping "$PYCOM_MAPPING" \
        --aa "$aa" \
        --out-dir "$OUT_DIR" 2>&1 | grep -E "(PDUs with|χ²|✓)"

    echo ""
done

echo "========================================================================"
echo "✓ All enrichment analyses complete"
echo "========================================================================"
echo ""
echo "Results saved to: $OUT_DIR"
echo ""
echo "Next: Run validation pipeline"
echo "  python3 bin/validate/hyperparameter_grid_parallel.py"
echo "  python3 bin/validate/validate_stability_parallel.py"
