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

# PDU Autoencoder Training + HDBSCAN Clustering on ARCHER2
# Assumes features already on ARCHER2
#
# Usage:
#   Single node, single AA (L):
#     sbatch submit_archer2.slurm
#
#   Different AA:
#     sbatch -J PDU_TRAIN_A --export=PDU_AA=A submit_archer2.slurm
#
#   Multiple nodes (for multi-threaded clustering):
#     sbatch -N 2 submit_archer2.slurm
#
#   Monitor:
#     squeue -u $USER
#     sacct -j <jobid> --format=JobID,JobName,State,Elapsed

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

# Threading (ARCHER2 CPU-friendly)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled

# Amino acid to process (default L)
AA="${PDU_AA:-L}"

# Directories (features must already be on ARCHER2)
FEATURES_DIR="pdu-extractor/analysis/features"
OUT_DIR="pdu-extractor/analysis"

# Create output directories
mkdir -p "$OUT_DIR/embeddings" "$OUT_DIR/clusters_umap" logs

echo "=========================================="
echo "ARCHER2: Train Autoencoder + Cluster"
echo "=========================================="
echo "AA: $AA"
echo "Nodes: $SLURM_NNODES"
echo "Tasks: $SLURM_NTASKS"
echo "Time limit: $SLURM_TIMELIMIT"
echo "=========================================="
echo ""

# Run combined training + clustering
python3 pdu-extractor/scripts/train_and_cluster_archer2.py \
    --aa "$AA" \
    --features-dir "$FEATURES_DIR" \
    --out-dir "$OUT_DIR" \
    --latent-dim 16 \
    --batch-size 1024 \
    --epochs 15 \
    --learning-rate 0.001 \
    --min-cluster-size 200 \
    --min-samples 10

# Move clusters to proper location
if [ -f "$OUT_DIR/pdu_clusters_${AA}_umap.csv" ]; then
    mv "$OUT_DIR/pdu_clusters_${AA}_umap.csv" "$OUT_DIR/clusters_umap/"
fi

if [ -f "$OUT_DIR/pdu_embedding_${AA}.npz" ]; then
    mv "$OUT_DIR/pdu_embedding_${AA}.npz" "$OUT_DIR/embeddings/"
fi

echo ""
echo "=========================================="
echo "✓ Complete: Training + Clustering done"
echo "=========================================="
