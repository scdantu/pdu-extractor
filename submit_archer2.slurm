#!/bin/bash
#SBATCH --job-name=PDU_TRAIN_L
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00

#SBATCH --account=e280-dantu
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --output=logs/pdu_train_%j.out
#SBATCH --error=logs/pdu_train_%j.err

# PDU Autoencoder Training on ARCHER2
#
# Usage:
#   Single node, single AA (L):
#     sbatch submit_archer2.slurm
#
#   Single node, different AA:
#     sbatch -J PDU_TRAIN_A --export=PDU_AA=A submit_archer2.slurm
#
#   Multiple nodes (4):
#     sbatch -N 4 submit_archer2.slurm
#
#   Monitor:
#     squeue -u $USER
#     sacct -j <jobid> --format=JobID,JobName,State,Elapsed,CPUTime

# Load modules
module load cray-python

# Setup virtual environment
export VENV_DIR="venv_archer2"
if [ ! -d "$VENV_DIR" ]; then
    echo "Error: Virtual environment not found at $VENV_DIR"
    echo "Run: bash pdu-extractor/setup_archer2_venv.sh"
    exit 1
fi

source "$VENV_DIR/bin/activate"

# Threading (ARCHER2 specific)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled

# Amino acid to train (default L, override with --export=PDU_AA=X)
AA="${PDU_AA:-L}"

# Features and output directories
FEATURES_DIR="pdu-extractor/analysis/features"
OUT_DIR="pdu-extractor/analysis/embeddings"

# Create output directories
mkdir -p "$OUT_DIR" logs

echo "=========================================="
echo "ARCHER2 PDU Autoencoder Training"
echo "=========================================="
echo "AA: $AA"
echo "Nodes: $SLURM_NNODES"
echo "Tasks: $SLURM_NTASKS"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Time limit: $SLURM_TIMELIMIT"
echo "Account: $SLURM_ACCOUNT"
echo "=========================================="

# Run training via srun (distributed across all allocated cores)
srun python3 pdu-extractor/scripts/train_pdu_autoencoder_archer2.py \
    --features-dir "$FEATURES_DIR" \
    --out-dir "$OUT_DIR" \
    --aa "$AA" \
    --batch-size 1024 \
    --epochs 15 \
    --learning-rate 0.001

echo "=========================================="
echo "✓ Training complete"
echo "=========================================="
