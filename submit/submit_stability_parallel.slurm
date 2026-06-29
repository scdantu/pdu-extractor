#!/bin/bash
#SBATCH --job-name=PDU_STABILITY_L
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

#SBATCH --account=e280-dantu
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --output=logs/stability_%a.out
#SBATCH --error=logs/stability_%a.err

#SBATCH --array=0-19

# ARCHER2 Parallel Stability Validation
#
# Runs 20 validation tasks in parallel:
#   Task 0-9: Bootstrap iterations (10 total)
#   Task 10-19: Perturbation iterations (10 total)
#
# Each task runs independently, results aggregated at the end.
#
# Usage:
#   sbatch submit_stability_parallel.slurm
#   sbatch -J PDU_STABILITY_A --export=PDU_AA=A submit_stability_parallel.slurm
#
# Monitor:
#   squeue --array
#   tail -f logs/stability_0.out

module load cray-python

export VENV_DIR="venv_archer2"
source "$VENV_DIR/bin/activate"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled

# Parameters
AA="${PDU_AA:-L}"
EMBEDDINGS_DIR="pdu-extractor/analysis/embeddings"
OUT_DIR="pdu-extractor/analysis/stability_parallel"

mkdir -p "$OUT_DIR" logs

# Task ID from array
TASK_ID=$SLURM_ARRAY_TASK_ID

# Map task ID to iteration type
if [ $TASK_ID -lt 10 ]; then
    ITER_TYPE="bootstrap"
    ITER_NUM=$TASK_ID
else
    ITER_TYPE="perturbation"
    ITER_NUM=$((TASK_ID - 10))
fi

echo "=========================================="
echo "ARCHER2: Parallel Stability Validation"
echo "=========================================="
echo "AA: $AA"
echo "Task ID: $TASK_ID / 19"
echo "Type: $ITER_TYPE (iteration $ITER_NUM)"
echo "CPUs: 128"
echo "=========================================="
echo ""

# Run stability validation for this task
python3 pdu-extractor/scripts/validate_stability_parallel.py \
    --task-id $TASK_ID \
    --embeddings "${EMBEDDINGS_DIR}/pdu_embedding_${AA}.npz" \
    --aa "$AA" \
    --out-dir "$OUT_DIR" \
    --min-cluster-size 200

echo ""
echo "=========================================="
echo "✓ Task $TASK_ID complete"
echo "=========================================="
