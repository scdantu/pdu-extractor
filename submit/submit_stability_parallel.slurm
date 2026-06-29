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

# Setup environment
source /work/e280/e280-Dantu/scdantu/Projects/pdus/pypdu/bin/activate

# Set Python path to include pdu-extractor root (for kmers module) + bin
export PYTHONPATH="/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor:/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor/bin:$PYTHONPATH"

# Working directory
export WORK_DIR="/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor"
export OUT_DIR="$WORK_DIR/analysis"

# Thread settings for CPU efficiency
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled

# Change to work directory
cd "$WORK_DIR" || exit 1

# Parameters
AA="${PDU_AA:-L}"
EMBEDDINGS_DIR="$OUT_DIR/embeddings"

# Create output directories
mkdir -p "$OUT_DIR/stability_parallel" logs

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
python3 bin/validate/validate_stability_parallel.py \
    --task-id $TASK_ID \
    --embeddings "${EMBEDDINGS_DIR}/pdu_embedding_${AA}.npz" \
    --aa "$AA" \
    --out-dir "$OUT_DIR/stability_parallel" \
    --min-cluster-size 200

echo ""
echo "=========================================="
echo "✓ Task $TASK_ID complete"
echo "=========================================="
