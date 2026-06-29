#!/bin/bash
# Run compression hypothesis test for a few amino acids

# Setup environment
source /work/e280/e280-Dantu/scdantu/Projects/pdus/pypdu/bin/activate
export PYTHONPATH="$PYTHONPATH:/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor/bin"
export WORK_DIR="/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor"
cd "$WORK_DIR" || exit 1

echo "=========================================================================="
echo "Running Compression Hypothesis Test"
echo "=========================================================================="
echo ""

# Test these 3 representative AAs
AAS=(L A G)  # Leucine (large), Alanine (small), Glycine (tiny)

for AA in "${AAS[@]}"; do
    echo "▶ Testing $AA..."
    ./pipelines/test_compression_hypothesis.sh "$AA" 2>&1 | tee "logs/compression_test_${AA}.log"
    echo ""
done

echo "=========================================================================="
echo "✓ All tests complete"
echo "=========================================================================="
echo ""
echo "Analyzing results..."
echo ""

for AA in "${AAS[@]}"; do
    echo "━━━ $AA ━━━"
    python3 bin/analyze/analyze_compression_test.py "$AA"
    echo ""
done

echo "Summary files:"
ls -lh analysis/compression_test/summary_*.csv

echo ""
echo "Plots:"
ls -lh analysis/compression_test/compression_analysis_*.png
