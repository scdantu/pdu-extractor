#!/bin/bash
#SBATCH --job-name=PDU_VALIDATION_L
#SBATCH --nodes=4
#SBATCH --ntasks=512
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00

#SBATCH --account=e280-dantu
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --output=logs/validation_%j.out
#SBATCH --error=logs/validation_%j.err

# Rigorous cluster validation pipeline on ARCHER2
#
# Tests:
#   1. Feature ablation (which features matter?)
#   2. Hyperparameter grid search (optimal latent_dim + min_cluster_size?)
#   3. Stability validation (bootstrap + perturbation)
#
# Usage:
#   sbatch submit_validation_archer2.slurm
#   sbatch -J PDU_VALIDATION_A --export=PDU_AA=A submit_validation_archer2.slurm

module load cray-python

export VENV_DIR="venv_archer2"
if [ ! -d "$VENV_DIR" ]; then
    echo "Error: Virtual environment not found at $VENV_DIR"
    exit 1
fi

source "$VENV_DIR/bin/activate"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled

# Parameters
AA="${PDU_AA:-L}"
DB_DIR="${DB_DIR:-.}"
FEATURES_DIR="pdu-extractor/analysis/features"
OUT_DIR="pdu-extractor/analysis/rigorous_validation"

mkdir -p "$OUT_DIR" logs

echo "=========================================="
echo "ARCHER2: Rigorous Cluster Validation"
echo "=========================================="
echo "AA: $AA"
echo "Nodes: $SLURM_NNODES"
echo "Tasks: $SLURM_NTASKS"
echo "Time limit: $SLURM_TIMELIMIT"
echo "=========================================="
echo ""

# Check inputs
DB_FILE="${DB_DIR}/pdus_${AA}.sqlite"
if [ ! -f "$DB_FILE" ]; then
    echo "Error: Database not found: $DB_FILE"
    exit 1
fi

FEATURES_FILE="${FEATURES_DIR}/pdu_features_${AA}.npz"
if [ ! -f "$FEATURES_FILE" ]; then
    echo "Error: Features not found: $FEATURES_FILE"
    exit 1
fi

echo "Using:"
echo "  DB: $DB_FILE"
echo "  Features: $FEATURES_FILE"
echo ""

# Run validation pipeline
python3 pdu-extractor/scripts/run_rigorous_validation.py \
    --db "$DB_FILE" \
    --features "$FEATURES_FILE" \
    --aa "$AA" \
    --out-dir "$OUT_DIR"

echo ""
echo "=========================================="
echo "✓ Validation pipeline complete"
echo "=========================================="
echo ""
echo "Results in: $OUT_DIR"
echo ""
