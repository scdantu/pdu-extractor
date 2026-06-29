#!/bin/bash
#SBATCH --job-name=PDU_COMPRESS
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00

#SBATCH --account=e280-dantu
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --output=logs/compress_%j.out
#SBATCH --error=logs/compress_%j.err

# ARCHER2: Compression Hypothesis Test
# Tests: 10Å vs 15Å × 8D to 256D latent dimensions
# AAs: L, A, G (representative subset)

# Setup environment
source /work/e280/e280-Dantu/scdantu/Projects/pdus/pypdu/bin/activate

# Set Python path to include pdu-extractor root (for kmers module) + bin
export PYTHONPATH="/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor:/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor/bin:$PYTHONPATH"

# Working directory
export WORK_DIR="/work/e280/e280-Dantu/scdantu/Projects/pdus/pdu-extractor"
export OUT_DIR="$WORK_DIR/analysis"

# Thread settings
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled

cd "$WORK_DIR" || exit 1

echo "=========================================================================="
echo "ARCHER2: Compression Hypothesis Test"
echo "=========================================================================="
echo "Node: $(hostname)"
echo "CPUs: 128"
echo "Time: 8 hours"
echo "=========================================================================="
echo ""

# Test 3 representative AAs
AAS=(G A L)  # Order: fastest to slowest
DISTANCES=(10 15)
DIMS=(8 16 32 64 128 256)
EPOCHS=10

mkdir -p "$OUT_DIR/compression_test" logs

for AA in "${AAS[@]}"; do
    echo "▶ Testing AA=$AA (distances: ${DISTANCES[@]}Å, dims: ${DIMS[@]}D)"
    echo ""

    OUT_SUBDIR="$OUT_DIR/compression_test"
    SUMMARY="$OUT_SUBDIR/summary_${AA}.csv"
    echo "Distance,LatentDim,Features,Embeddings,Clusters,AvgConf,NoisePercent,CATHchi2" > "$SUMMARY"

    for DIST in "${DISTANCES[@]}"; do
        echo "  Distance: ${DIST}Å"

        # Extract features (one-time per distance)
        FEATURES_FILE="$OUT_SUBDIR/features_${AA}_${DIST}a.npz"
        if [ ! -f "$FEATURES_FILE" ]; then
            echo "    → Extracting features..."
            python3 bin/extract/export_pdu_features.py \
                --db pdb_db/pdus_${AA}.sqlite \
                --out-file "$FEATURES_FILE" \
                --distance-cutoff $DIST \
                --aa "$AA" \
                2>&1 | grep -E "(INFO|ERROR|✓)" || true
        else
            echo "    → Features exist (skipping)"
        fi

        # Test each dimension
        for DIM in "${DIMS[@]}"; do
            echo "    Testing: ${DIM}D"

            # Train autoencoder
            EMBEDDING_FILE="$OUT_SUBDIR/embedding_${AA}_${DIST}a_${DIM}d.npz"
            if [ ! -f "$EMBEDDING_FILE" ]; then
                python3 bin/train/train_pdu_autoencoder.py \
                    --features "$FEATURES_FILE" \
                    --latent-dim "$DIM" \
                    --epochs "$EPOCHS" \
                    --batch-size 32 \
                    --out-file "$EMBEDDING_FILE" \
                    --aa "$AA" \
                    2>&1 | grep -E "(Epoch|Loss|✓)" | tail -3 || true
            fi

            # Cluster
            CLUSTER_FILE="$OUT_SUBDIR/clusters_${AA}_${DIST}a_${DIM}d.csv"
            if [ ! -f "$CLUSTER_FILE" ]; then
                python3 bin/cluster/cluster_pdu_embeddings_gpu.py \
                    --embeddings "$EMBEDDING_FILE" \
                    --out-file "$CLUSTER_FILE" \
                    --min-cluster-size 200 \
                    --aa "$AA" \
                    2>&1 | grep -E "(clusters|noise|✓)" || true
            fi

            # Get stats
            STATS=$(python3 << PYEOF
import pandas as pd
df = pd.read_csv('$CLUSTER_FILE')
assigned = df[df['cluster'] != -1]
noise_pct = 100 * (df['cluster'] == -1).sum() / len(df)
avg_conf = assigned['cluster_probability'].mean() if len(assigned) > 0 else 0
print(f"{noise_pct:.1f},{avg_conf:.3f}")
PYEOF
)
            echo "$DIST,$DIM,$FEATURES_FILE,$EMBEDDING_FILE,$CLUSTER_FILE,$(echo $STATS),?" >> "$SUMMARY"
        done
    done

    echo ""
    echo "  Results for $AA:"
    cat "$SUMMARY" | tail -7 | column -t -s','
    echo ""
done

echo "=========================================================================="
echo "✓ Testing complete"
echo "=========================================================================="
echo ""
echo "Analyzing results..."
for AA in "${AAS[@]}"; do
    python3 bin/analyze/analyze_compression_test.py "$AA" 2>&1 | grep -E "(FINDING|CONCLUSION|→)" || true
done

echo ""
echo "Summary files: $OUT_DIR/compression_test/summary_*.csv"
echo "Plot files: $OUT_DIR/compression_test/compression_analysis_*.png"
