#!/bin/bash
#SBATCH --job-name=COMPRESS_UNIFIED
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --account=e280-dantu
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --output=logs/compress_unified_%j.out
#SBATCH --error=logs/compress_unified_%j.err

# ARCHER2: Unified Compression Hypothesis Test (PARALLEL VERSION)
# Tests encoding granularity, distance cutoff, and compression jointly
#
# PARALLELIZATION: 3 AAs × 3 encodings = 9 parallel jobs
# Each job tests: 2 distances × 6 latent dims = 12 tests
# Total: 108 tests
#
# Sequential: ~24 hours
# Parallel (9 jobs): ~2-3 hours (108 ÷ 9 = 12 tests per process)

# Setup environment
source /work/e280/e280-Dantu/scdantu/Projects/pdus/pypdu/bin/activate

# Set Python path for pdusearch library
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
echo "ARCHER2: Unified Compression Hypothesis Test"
echo "=========================================================================="
echo "Node: $(hostname)"
echo "CPUs: 128"
echo "Time: 12 hours"
echo "Repository: pdusearch library (new modular design)"
echo "=========================================================================="
echo ""

# Test 3 representative AAs
AAS=(G A L)
ENCODINGS=(aa_20 functional_5 aa_20_conservation)

# Create output directory
mkdir -p "$OUT_DIR" logs

echo "▶ Launching parallel jobs: 3 AAs × 3 encodings = 9 processes"
echo ""

# Launch all AA × encoding combinations in parallel
declare -a PIDS
declare -a JOB_LABELS

job_count=0
for AA in "${AAS[@]}"; do
    for ENC in "${ENCODINGS[@]}"; do
        job_count=$((job_count + 1))

        echo "[$job_count/9] Launching: AA=$AA, encoding=$ENC"

        # Run test in background
        python3 bin/compression_test_unified.py \
            --aa "$AA" \
            --encodings "$ENC" \
            --distances 10 15 \
            --dims 8 16 32 64 128 256 \
            --epochs 10 \
            --out-dir "$OUT_DIR" \
            --log-level INFO \
            --log-file "$OUT_DIR/test_${AA}_${ENC}.log" &

        PIDS+=($!)
        JOB_LABELS+=("$AA-$ENC")
    done
done

echo ""
echo "Waiting for all 9 jobs to complete..."
echo "Monitor with: tail -f logs/compress_unified_*.out"
echo ""

# Wait for all background jobs and track completion
completed=0
for i in "${!PIDS[@]}"; do
    pid=${PIDS[$i]}
    label=${JOB_LABELS[$i]}

    if wait $pid; then
        echo "✓ [$((i+1))/9] $label complete"
        completed=$((completed + 1))
    else
        echo "✗ [$((i+1))/9] $label failed"
    fi
done

echo ""
echo "=========================================================================="
echo "✓ All tests complete: $completed/9 jobs successful"
echo "=========================================================================="
echo ""
echo "Output files:"
ls -lh "$OUT_DIR"/summary_*.csv 2>/dev/null || echo "  (No results yet)"

echo ""
echo "Analyzing results for each AA:"
for AA in "${AAS[@]}"; do
    if [ -f "$OUT_DIR/summary_${AA}_aa_20.csv" ]; then
        echo "  ▶ Analyzing $AA..."
        python3 bin/analyze_compression_results.py --aa "$AA" 2>&1 | grep -E "^(Best|Noise|Clusters|→)" || true
    fi
done

echo ""
echo "Full analysis:"
echo "  python3 bin/analyze_compression_results.py --aa G"
echo "  python3 bin/analyze_compression_results.py --aa A"
echo "  python3 bin/analyze_compression_results.py --aa L"
echo ""
echo "Results saved to: $OUT_DIR/"
