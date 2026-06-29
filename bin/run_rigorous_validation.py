#!/usr/bin/env python3
"""
Master orchestrator: Run full rigorous validation pipeline.

Pipeline:
  1. Feature ablation: Test AA-only, SS-only, distance-only vs full
  2. Hyperparameter grid search: Find optimal latent_dim + min_cluster_size
  3. Stability validation: Bootstrap + Perturbation + Silhouette

Usage:
  python3 run_rigorous_validation.py \
    --db /path/to/pdus_L.sqlite \
    --features analysis/features/pdu_features_L.npz \
    --aa L \
    --out-dir analysis/rigorous_validation
"""

import subprocess
import argparse
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def run_command(cmd, description):
    """Run shell command and log."""
    logger.info(f"\n{'='*80}")
    logger.info(f"STEP: {description}")
    logger.info(f"{'='*80}")
    logger.info(f"Running: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        logger.error(f"Failed: {description}")
        return False

    logger.info(f"✓ Complete: {description}")
    return True


def main():
    parser = argparse.ArgumentParser(description="Rigorous validation pipeline")
    parser.add_argument("--db", required=True, help="SQLite database")
    parser.add_argument("--features", required=True, help="Features NPZ file")
    parser.add_argument("--aa", required=True, help="Amino acid")
    parser.add_argument("--out-dir", default="analysis/rigorous_validation", help="Output directory")
    parser.add_argument("--skip-ablation", action="store_true", help="Skip feature ablation")
    parser.add_argument("--skip-grid", action="store_true", help="Skip hyperparameter grid")
    parser.add_argument("--skip-stability", action="store_true", help="Skip stability validation")

    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info("RIGOROUS VALIDATION PIPELINE")
    logger.info("=" * 80)
    logger.info(f"AA: {args.aa}")
    logger.info(f"Output: {args.out_dir}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    scripts_dir = Path(__file__).parent

    # Step 1: Feature ablation
    if not args.skip_ablation:
        logger.info("\n" + "=" * 80)
        logger.info("PHASE 1: FEATURE ABLATION")
        logger.info("=" * 80)

        for subset in ["full", "aa_only", "ss_only", "dist_only"]:
            cmd = [
                "python3", str(scripts_dir / "feature_ablation_archer2.py"),
                "--db", args.db,
                "--aa", args.aa,
                "--subset", subset,
                "--out-dir", str(out_dir / "ablation"),
            ]
            if not run_command(cmd, f"Feature ablation: {subset}"):
                return

    # Step 2: Hyperparameter grid search
    if not args.skip_grid:
        logger.info("\n" + "=" * 80)
        logger.info("PHASE 2: HYPERPARAMETER GRID SEARCH")
        logger.info("=" * 80)

        cmd = [
            "python3", str(scripts_dir / "hyperparameter_grid_search.py"),
            "--features", args.features,
            "--aa", args.aa,
            "--out-dir", str(out_dir / "grid_search"),
            "--latent-dims", "8", "16", "32", "64",
            "--min-cluster-sizes", "50", "100", "200", "300", "500",
            "--epochs", "5",  # Quick search with 5 epochs
        ]
        if not run_command(cmd, "Hyperparameter grid search"):
            return

    # Step 3: Stability validation (with optimal hyperparameters)
    if not args.skip_stability:
        logger.info("\n" + "=" * 80)
        logger.info("PHASE 3: STABILITY VALIDATION")
        logger.info("=" * 80)

        cmd = [
            "python3", str(scripts_dir / "validate_clusters_stability.py"),
            "--embeddings", args.features.replace("pdu_features", "pdu_embedding").replace(".npz", "_full.npz"),
            "--aa", args.aa,
            "--n-bootstrap", "20",
            "--n-perturbations", "10",
            "--min-cluster-size", "200",
            "--out-dir", str(out_dir / "stability"),
        ]
        if not run_command(cmd, "Stability validation"):
            return

    logger.info("\n" + "=" * 80)
    logger.info("✓ RIGOROUS VALIDATION COMPLETE")
    logger.info("=" * 80)
    logger.info(f"\nResults in: {out_dir}")
    logger.info("\nNext steps:")
    logger.info("  1. Review ablation results → which features matter?")
    logger.info("  2. Review grid search → optimal hyperparameters?")
    logger.info("  3. Review stability → are clusters robust?")
    logger.info("  4. If all checks pass → ready for publication")
    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
