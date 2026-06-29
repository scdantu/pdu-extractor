#!/bin/bash
# Batch CATH/EC enrichment analysis for multiple amino acids

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

DB_DIR="${DB_DIR:-/Users/sarathdantu/Library/CloudStorage/OneDrive-BrunelUniversityLondon/Work/Research/Manuscript/Me/3D\ Sequence\ Evolution/data/PDUs/pdu_db}"
CLUSTERS_DIR="${CLUSTERS_DIR:-analysis/clusters_umap}"
OUT_DIR="${OUT_DIR:-analysis/cath_ec_enrichment}"

mkdir -p "${OUT_DIR}"

echo "=========================================="
echo "BATCH CATH/EC ENRICHMENT ANALYSIS"
echo "=========================================="
echo ""

AAS=(L A E G)
CACHE_DIR="${OUT_DIR}/.cache"
mkdir -p "${CACHE_DIR}"

for AA in "${AAS[@]}"; do
    DB="${DB_DIR}/pdus_${AA}.sqlite"
    CLUSTERS="${CLUSTERS_DIR}/pdu_clusters_${AA}_umap.csv"
    CATH_CACHE="${CACHE_DIR}/cath_cache_${AA}.json"
    EC_CACHE="${CACHE_DIR}/ec_cache_${AA}.json"

    if [[ ! -f "${DB}" ]]; then
        echo "⚠️  Database not found: ${DB}"
        continue
    fi

    if [[ ! -f "${CLUSTERS}" ]]; then
        echo "⚠️  Clusters not found: ${CLUSTERS}"
        continue
    fi

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "AA=${AA}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    python3 scripts/analyze_cath_ec_enrichment.py \
        --clusters "${CLUSTERS}" \
        --db "${DB}" \
        --aa "${AA}" \
        --out-dir "${OUT_DIR}" \
        --cath-cache "${CATH_CACHE}" \
        --ec-cache "${EC_CACHE}"

    echo ""
done

echo "=========================================="
echo "✓ Analysis complete"
echo "Results in: ${OUT_DIR}/"
echo "=========================================="
