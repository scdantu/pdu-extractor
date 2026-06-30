#!/bin/bash
#SBATCH --job-name=COMPRESS_DDP
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --account=e280-dantu
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --output=logs/compress_ddp_%j.out
#SBATCH --error=logs/compress_ddp_%j.err

# ARCHER2: Compression Test with Distributed Data Parallel (DDP)
#
# Uses torchrun to parallelize training across all 128 CPUs
# Tests all 3 AAs × 3 encodings × 2 distances × 6 dims = 108 configs
# Expected time: ~30-45 minutes (vs 4+ hours without DDP)
#
# DDP provides: 128x CPU utilization + automatic gradient synchronization

# Setup environment
source /work/e280/e280-Dantu/scdantu/Projects/pdus/pypdu/bin/activate

# Set Python path
export PYTHONPATH="/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor:$PYTHONPATH"

# Working directory
export WORK_DIR="/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor"
export OUT_DIR="$WORK_DIR/analysis/compression_test"

# Thread settings
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled

cd "$WORK_DIR" || exit 1

echo "=========================================================================="
echo "ARCHER2: Compression Test with DDP"
echo "=========================================================================="
echo "Node: $(hostname)"
echo "CPUs: 128"
echo "Time: 2 hours"
echo "Method: Distributed Data Parallel (torchrun)"
echo "Expected duration: 30-45 minutes"
echo "=========================================================================="
echo ""

# Create output directory
mkdir -p "$OUT_DIR" logs

# Test 3 representative AAs
AAS=(G A L)
ENCODINGS=(aa_20 functional_5 aa_20_conservation)

echo "▶ Launching DDP training: 3 AAs × 3 encodings = 9 sequential jobs"
echo "   Each job runs with all 128 CPUs in parallel"
echo ""

for AA in "${AAS[@]}"; do
    for ENC in "${ENCODINGS[@]}"; do
        echo "Launching: AA=$AA, encoding=$ENC"

        # Use torchrun to distribute across 128 CPUs
        # Each job trains on full dataset split across 128 processes
        torchrun --nproc_per_node=128 \
            bin/compression_test_ddp.py \
            --aa "$AA" \
            --encodings "$ENC" \
            --distances 10 15 \
            --dims 8 16 32 64 128 256 \
            --epochs 10 \
            --out-dir "$OUT_DIR" \
            --log-level INFO \
            --log-file "$OUT_DIR/test_${AA}_${ENC}_ddp.log"

        if [ $? -eq 0 ]; then
            echo "✓ $AA-$ENC complete"
        else
            echo "✗ $AA-$ENC failed"
        fi

        echo ""
    done
done

echo "=========================================================================="
echo "✓ All DDP training complete"
echo "=========================================================================="
echo ""
echo "Output files:"
ls -lh "$OUT_DIR"/summary_*.csv 2>/dev/null || echo "  (No results)"

echo ""
echo "Analyzing results:"
for AA in "${AAS[@]}"; do
    if [ -f "$OUT_DIR/summary_${AA}_aa_20.csv" ]; then
        echo "  ▶ $AA analysis:"
        python3 bin/analyze_compression_results.py --aa "$AA" 2>&1 | head -20
        echo ""
    fi
done

echo "Full analysis commands:"
echo "  python3 bin/analyze_compression_results.py --aa G"
echo "  python3 bin/analyze_compression_results.py --aa A"
echo "  python3 bin/analyze_compression_results.py --aa L"
echo ""
echo "Results saved to: $OUT_DIR/"
