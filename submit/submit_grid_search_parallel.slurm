#!/bin/bash
#SBATCH --job-name=PDU_GRID_L
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

#SBATCH --account=e280-dantu
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --output=logs/grid_%a.out
#SBATCH --error=logs/grid_%a.err

#SBATCH --array=0-19

# ARCHER2 Parallel Hyperparameter Grid Search
#
# Runs 20 grid search tasks in parallel (4 at a time on 2 nodes each).
# Each task: different (latent_dim, min_cluster_size) combination.
#
# Grid:
#   latent_dim: [8, 16, 32, 64]
#   min_cluster_size: [50, 100, 200, 300, 500]
#   Total: 4 × 5 = 20 combinations
#
# Usage:
#   Single AA (L):
#     sbatch submit_grid_search_parallel.slurm
#
#   Different AA (A):
#     sbatch -J PDU_GRID_A --export=PDU_AA=A submit_grid_search_parallel.slurm
#
# Monitor:
#   squeue --array
#   sacct -j <jobid> --format=JobID,JobName,State,Elapsed

# Setup environment
source /work/e280/e280-Dantu/scdantu/Projects/pdus/pypdu/bin/activate

# Set Python path to include pdu-extractor scripts
export PYTHONPATH="$PYTHONPATH:/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor/bin"

# Working directory
export WORK_DIR="/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor"
export OUT_DIR="$WORK_DIR/analysis"
export FEATURES_DIR="$OUT_DIR/features"

# Thread settings for CPU efficiency
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled

# Change to work directory
cd "$WORK_DIR" || exit 1

# Parameters
AA="${PDU_AA:-L}"

# Create output directories
mkdir -p "$OUT_DIR/embeddings" "$OUT_DIR/clusters_umap" "$OUT_DIR/grid_search_parallel" logs

# Task ID from array
TASK_ID=$SLURM_ARRAY_TASK_ID

echo "=========================================="
echo "ARCHER2: Parallel Grid Search"
echo "=========================================="
echo "AA: $AA"
echo "Task ID: $TASK_ID / 19"
echo "Nodes: $SLURM_NNODES"
echo "CPUs per node: 128"
echo "=========================================="
echo ""

# Run grid search for this task
python3 bin/validate/hyperparameter_grid_parallel.py \
    --task-id $TASK_ID \
    --features "${FEATURES_DIR}/pdu_features_${AA}.npz" \
    --aa "$AA" \
    --out-dir "$OUT_DIR/grid_search_parallel" \
    --epochs 5

echo ""
echo "=========================================="
echo "✓ Task $TASK_ID complete"
echo "=========================================="
