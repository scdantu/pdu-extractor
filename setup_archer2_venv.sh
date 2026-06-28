#!/bin/bash
# Setup Python virtual environment on ARCHER2
# Run on login node:
#   bash pdu-extractor/setup_archer2_venv.sh

set -e

echo "=========================================="
echo "ARCHER2 Virtual Environment Setup"
echo "=========================================="
echo ""

# Check current directory
if [ ! -f "pdu-extractor/requirements-archer2.txt" ]; then
    echo "Error: Run this script from the parent directory containing pdu-extractor/"
    echo "  cd /path/to/mdpdu"
    echo "  bash pdu-extractor/setup_archer2_venv.sh"
    exit 1
fi

# Load Python module
echo "[1/5] Loading Python module..."
module load cray-python
python3 --version

# Create venv
VENV_DIR="venv_archer2"
echo "[2/5] Creating virtual environment: $VENV_DIR"
python3 -m venv "$VENV_DIR"

# Activate venv
echo "[3/5] Activating virtual environment..."
source "$VENV_DIR/bin/activate"

# Upgrade pip
echo "[4/5] Upgrading pip..."
pip install --upgrade pip setuptools wheel

# Install requirements
echo "[5/5] Installing ARCHER2-optimized requirements..."
echo "      This may take 5-10 minutes..."
pip install -r pdu-extractor/requirements-archer2.txt

echo ""
echo "=========================================="
echo "✓ Virtual environment ready"
echo "=========================================="
echo ""
echo "To activate on ARCHER2:"
echo "  source venv_archer2/bin/activate"
echo ""
echo "To use in Slurm jobs, load modules then activate:"
echo "  module load cray-python"
echo "  source venv_archer2/bin/activate"
echo ""
