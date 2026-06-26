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

DB="${DB:-pdu_output/pdus.sqlite}"
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
SPACE="${SPACE:-umap}"
MIN_CLUSTER_SIZE="${MIN_CLUSTER_SIZE:-50}"
MIN_SAMPLES="${MIN_SAMPLES:-10}"
STANDARDIZE="${STANDARDIZE:-false}"
LOG_LEVEL="${LOG_LEVEL:-INFO}"

mkdir -p "${LOG_DIR}"
PIPELINE_LOG="${LOG_DIR}/${AA}_pipeline.log"
exec > >(tee -a "${PIPELINE_LOG}") 2>&1

echo "$(date '+%Y-%m-%d %H:%M:%S') INFO run_aa_pipeline: AA=${AA}"
echo "$(date '+%Y-%m-%d %H:%M:%S') INFO run_aa_pipeline: DB=${DB}"

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
