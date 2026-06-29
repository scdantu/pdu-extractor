#!/bin/bash
# Test compression hypothesis: 15Å features × varying latent dimensions
#
# Hypothesis: 15Å might be fine, but 900D → 16D is TOO aggressive
# Test: 10Å vs 15Å × 8D, 16D, 32D, 64D, 128D, 256D
#
# Setup environment
source /work/e280/e280-Dantu/scdantu/Projects/pdus/pypdu/bin/activate
export PYTHONPATH="$PYTHONPATH:/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor/bin"
export WORK_DIR="/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor"
cd "$WORK_DIR" || exit 1

# Parameters
AA="${1:-L}"  # Default: Leucine
DISTANCES=(10 15)
DIMS=(8 16 32 64 128 256)
EPOCHS=10

OUT_DIR="$WORK_DIR/analysis/compression_test"
mkdir -p "$OUT_DIR" logs

echo "========================================================================"
echo "Testing Compression Hypothesis"
echo "========================================================================"
echo "AA: $AA"
echo "Distances: ${DISTANCES[@]} Å"
echo "Latent dims: ${DIMS[@]}"
echo "Epochs: $EPOCHS"
echo "========================================================================"
echo ""

# Results summary file
SUMMARY="$OUT_DIR/summary_${AA}.csv"
echo "Distance,LatentDim,Features,Embeddings,Clusters,AvgConf,NoisePercent,CATHchi2" > "$SUMMARY"

# Test each combination
for DIST in "${DISTANCES[@]}"; do
    echo "======== Distance: ${DIST}Å ========"

    # Extract features (one-time per distance)
    FEATURES_FILE="$OUT_DIR/features_${AA}_${DIST}a.npz"
    if [ ! -f "$FEATURES_FILE" ]; then
        echo "→ Extracting features (${DIST}Å)..."
        python3 bin/extract/export_pdu_features.py \
            --db pdb_db/pdus_${AA}.sqlite \
            --out-file "$FEATURES_FILE" \
            --distance-cutoff $DIST \
            --aa "$AA" \
            2>&1 | tee "logs/extract_${AA}_${DIST}a.log"
    else
        echo "→ Features exist: $FEATURES_FILE"
    fi

    # Test each dimension
    for DIM in "${DIMS[@]}"; do
        echo "  Testing: ${DIST}Å × ${DIM}D"

        # Train autoencoder
        EMBEDDING_FILE="$OUT_DIR/embedding_${AA}_${DIST}a_${DIM}d.npz"
        if [ ! -f "$EMBEDDING_FILE" ]; then
            echo "    → Training autoencoder (${DIM}D, $EPOCHS epochs)..."
            python3 bin/train/train_pdu_autoencoder.py \
                --features "$FEATURES_FILE" \
                --latent-dim "$DIM" \
                --epochs "$EPOCHS" \
                --batch-size 32 \
                --out-file "$EMBEDDING_FILE" \
                --aa "$AA" \
                2>&1 | tee "logs/train_${AA}_${DIST}a_${DIM}d.log"
        fi

        # Cluster
        CLUSTER_FILE="$OUT_DIR/clusters_${AA}_${DIST}a_${DIM}d.csv"
        if [ ! -f "$CLUSTER_FILE" ]; then
            echo "    → Clustering..."
            python3 bin/cluster/cluster_pdu_embeddings_gpu.py \
                --embeddings "$EMBEDDING_FILE" \
                --out-file "$CLUSTER_FILE" \
                --min-cluster-size 200 \
                --aa "$AA" \
                2>&1 | tee "logs/cluster_${AA}_${DIST}a_${DIM}d.log"
        fi

        # Validate (quick stats)
        echo "    → Computing stats..."
        STATS=$(python3 << PYEOF
import pandas as pd
import numpy as np
df = pd.read_csv('$CLUSTER_FILE')
assigned = df[df['cluster'] != -1]
noise_pct = 100 * (df['cluster'] == -1).sum() / len(df)
avg_conf = assigned['cluster_probability'].mean() if len(assigned) > 0 else 0
print(f"{noise_pct:.1f},{avg_conf:.3f}")
PYEOF
)

        NOISE=$(echo "$STATS" | cut -d, -f1)
        CONF=$(echo "$STATS" | cut -d, -f2)

        echo "$DIST,$DIM,$FEATURES_FILE,$EMBEDDING_FILE,$CLUSTER_FILE,$CONF,$NOISE,?" >> "$SUMMARY"
        echo "    ✓ Noise: ${NOISE}%, Conf: ${CONF}"
    done

    echo ""
done

echo "========================================================================"
echo "✓ Testing complete"
echo "========================================================================"
echo "Summary: $SUMMARY"
echo ""
echo "Results:"
cat "$SUMMARY"
echo ""
echo "Next: Analyze results to find optimal distance × dimension combination"
echo "========================================================================"
