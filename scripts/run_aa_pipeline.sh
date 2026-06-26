#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

AA="${1:-${AA:-}}"
if [[ -z "${AA}" ]]; then
    echo "Usage: scripts/run_aa_pipeline.sh <AA>"
    echo "Example: scripts/run_aa_pipeline.sh A"
    exit 1
fi

# Use per-AA database (hardcoded path)
DB_DIR="${DB_DIR:-/mnt/disk_b/pdu-extractor/per_aa_sqlite}"
DB="${DB_DIR}/pdus_${AA}.sqlite"

# Verify database exists
if [[ ! -f "${DB}" ]]; then
    echo "ERROR: Database not found at ${DB}"
    echo "DB_DIR=${DB_DIR}"
    ls -la "${DB_DIR}"/ | head -10
    exit 1
fi
FEATURES_DIR="${FEATURES_DIR:-analysis/features}"
EMBEDDINGS_DIR="${EMBEDDINGS_DIR:-analysis/embeddings}"
CLUSTERS_DIR="${CLUSTERS_DIR:-analysis/clusters_umap}"
LOG_DIR="${LOG_DIR:-logs/aa_pipeline}"
RADIUS="${RADIUS:-15.0}"
BIN_WIDTH="${BIN_WIDTH:-1.0}"
RESIDUE_ENCODING="${RESIDUE_ENCODING:-aa}"
EPOCHS="${EPOCHS:-15}"
BATCH_SIZE="${BATCH_SIZE:-256}"
LATENT_DIM="${LATENT_DIM:-16}"
LEARNING_RATE="${LEARNING_RATE:-0.001}"
VALIDATION_FRACTION="${VALIDATION_FRACTION:-0.1}"
PATIENCE="${PATIENCE:-0}"
SPACE="${SPACE:-umap}"

# Tune min_cluster_size per amino acid based on PDU count (from db_stats)
case "${AA}" in
    L|A) MIN_CLUSTER_SIZE="${MIN_CLUSTER_SIZE:-200}" ;;  # Large classes (1.9–2.1M PDUs)
    G|V|E|S|D|T|K|I) MIN_CLUSTER_SIZE="${MIN_CLUSTER_SIZE:-80}" ;;  # Medium classes (1.2–1.8M PDUs)
    R|P|N|F|Q|Y) MIN_CLUSTER_SIZE="${MIN_CLUSTER_SIZE:-50}" ;;  # Small-medium (860k–1.1M PDUs)
    H|M) MIN_CLUSTER_SIZE="${MIN_CLUSTER_SIZE:-30}" ;;  # Small classes (515k–575k PDUs)
    W|C) MIN_CLUSTER_SIZE="${MIN_CLUSTER_SIZE:-15}" ;;  # Tiny classes (360k–363k PDUs)
    *) MIN_CLUSTER_SIZE="${MIN_CLUSTER_SIZE:-50}" ;;
esac

MIN_SAMPLES="${MIN_SAMPLES:-10}"
STANDARDIZE="${STANDARDIZE:-false}"
LOG_LEVEL="${LOG_LEVEL:-INFO}"

mkdir -p "${LOG_DIR}"
PIPELINE_LOG="${LOG_DIR}/${AA}_pipeline.log"
exec > >(tee -a "${PIPELINE_LOG}") 2>&1

echo "$(date '+%Y-%m-%d %H:%M:%S') INFO run_aa_pipeline: AA=${AA}"
echo "$(date '+%Y-%m-%d %H:%M:%S') INFO run_aa_pipeline: DB=${DB}"

# Log GPU info if using multiple GPUs
if [[ -n "${CUDA_VISIBLE_DEVICES:-}" ]]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') INFO run_aa_pipeline: CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
fi

export PYTHONPATH="${REPO_DIR}:${PYTHONPATH:-}"

python scripts/export_pdu_features.py \
    --db "${DB}" \
    --out-dir "${FEATURES_DIR}" \
    --aa "${AA}" \
    --radius "${RADIUS}" \
    --bin-width "${BIN_WIDTH}" \
    --residue-encoding "${RESIDUE_ENCODING}" \
    --log-file "${LOG_DIR}/${AA}_features.log" \
    --log-level "${LOG_LEVEL}"

python scripts/train_pdu_autoencoder.py \
    --features-dir "${FEATURES_DIR}" \
    --out-dir "${EMBEDDINGS_DIR}" \
    --aa "${AA}" \
    --latent-dim "${LATENT_DIM}" \
    --epochs "${EPOCHS}" \
    --batch-size "${BATCH_SIZE}" \
    --learning-rate "${LEARNING_RATE}" \
    --validation-fraction "${VALIDATION_FRACTION}" \
    --patience "${PATIENCE}" \
    --log-file "${LOG_DIR}/${AA}_training.log" \
    --log-level "${LOG_LEVEL}"

cluster_cmd=(
    python scripts/cluster_pdu_embeddings.py
    --embeddings-dir "${EMBEDDINGS_DIR}"
    --out-dir "${CLUSTERS_DIR}"
    --aa "${AA}"
    --space "${SPACE}"
    --min-cluster-size "${MIN_CLUSTER_SIZE}"
    --min-samples "${MIN_SAMPLES}"
    --log-file "${LOG_DIR}/${AA}_clustering.log"
    --log-level "${LOG_LEVEL}"
)

if [[ "${STANDARDIZE}" == "true" ]]; then
    cluster_cmd+=(--standardize)
fi

"${cluster_cmd[@]}"

echo "$(date '+%Y-%m-%d %H:%M:%S') INFO run_aa_pipeline: completed AA=${AA}"
