#!/usr/bin/env bash
# Run full PDU clustering pipeline across all amino acids
# Usage: scripts/run_full_pipeline.sh [--parallel N] [--db-dir path]

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

PARALLEL=1
DB_DIR="pdu_db"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --parallel)
            PARALLEL="$2"
            shift 2
            ;;
        --db-dir)
            DB_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# All amino acids (order: from largest to smallest for better pipelining)
AAS=(L A G V E S D T K I R P N F Q Y H M W C)

echo "$(date '+%Y-%m-%d %H:%M:%S') Starting full PDU pipeline (${#AAS[@]} amino acids, ${PARALLEL} parallel)"
echo "Database directory: ${DB_DIR}"

# Create log directory
mkdir -p logs/aa_pipeline

# Run with GNU parallel if available, otherwise serial
if command -v parallel &>/dev/null && [[ ${PARALLEL} -gt 1 ]]; then
    echo "Using GNU parallel with ${PARALLEL} jobs"
    export DB_DIR
    export REPO_DIR
    export PYTHONPATH="${REPO_DIR}:${PYTHONPATH:-}"
    # Distribute GPUs round-robin (requires nvidia-smi, optional)
    if command -v nvidia-smi &>/dev/null; then
        N_GPUS=$(nvidia-smi --list-gpus | wc -l)
        echo "Detected ${N_GPUS} GPU(s) - distributing across parallel jobs"
        printf '%s\n' "${AAS[@]}" | parallel -j "${PARALLEL}" \
            "CUDA_VISIBLE_DEVICES={%${N_GPUS}} scripts/run_aa_pipeline.sh {}"
    else
        printf '%s\n' "${AAS[@]}" | parallel -j "${PARALLEL}" \
            "scripts/run_aa_pipeline.sh {}"
    fi
elif [[ ${PARALLEL} -gt 1 ]]; then
    echo "GNU parallel not found, using xargs with ${PARALLEL} jobs"
    export DB_DIR
    export REPO_DIR
    export PYTHONPATH="${REPO_DIR}:${PYTHONPATH:-}"
    printf '%s\n' "${AAS[@]}" | xargs -P "${PARALLEL}" -I {} bash -c 'scripts/run_aa_pipeline.sh "$@"' _ {}
else
    echo "Running serially (1 AA at a time)"
    for AA in "${AAS[@]}"; do
        echo "$(date '+%Y-%m-%d %H:%M:%S') Running AA=${AA}"
        DB_DIR="${DB_DIR}" bash scripts/run_aa_pipeline.sh "${AA}" || {
            echo "ERROR: Pipeline failed for AA=${AA}"
            exit 1
        }
    done
fi

echo "$(date '+%Y-%m-%d %H:%M:%S') Pipeline complete"
echo "Results:"
echo "  - Features: analysis/features/"
echo "  - Embeddings: analysis/embeddings/"
echo "  - Clusters: analysis/clusters_umap/"
echo "  - Logs: logs/aa_pipeline/"
