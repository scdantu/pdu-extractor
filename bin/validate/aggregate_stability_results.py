#!/usr/bin/env python3
"""
Aggregate parallel stability validation results.

Reads all task_XX_*.csv files and creates summary with statistics.

Usage:
  python3 aggregate_stability_results.py \
    --stability-dir analysis/stability_parallel \
    --aa L
"""

import argparse
from pathlib import Path
import logging

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Aggregate stability results")
    parser.add_argument("--stability-dir", required=True, help="Stability results directory")
    parser.add_argument("--aa", required=True, help="Amino acid")
    parser.add_argument("--out-dir", default=None, help="Output directory")

    args = parser.parse_args()

    stability_dir = Path(args.stability_dir)
    out_dir = Path(args.out_dir) if args.out_dir else stability_dir

    logger.info("=" * 80)
    logger.info("Aggregating Stability Validation Results")
    logger.info("=" * 80)
    logger.info(f"Input: {stability_dir}")

    # Find all task files
    task_files = sorted(stability_dir.glob("task_*.csv"))
    logger.info(f"Found {len(task_files)} task results")

    if not task_files:
        logger.error(f"No task files found in {stability_dir}")
        return

    # Aggregate results
    all_results = []
    for task_file in task_files:
        df = pd.read_csv(task_file)
        all_results.append(df)

    results_df = pd.concat(all_results, ignore_index=True)

    # Separate bootstrap and perturbation
    bootstrap_df = results_df[results_df['type'] == 'bootstrap']
    perturbation_df = results_df[results_df['type'] == 'perturbation']

    logger.info(f"Bootstrap: {len(bootstrap_df)} iterations")
    logger.info(f"Perturbation: {len(perturbation_df)} iterations")

    # Statistics
    bootstrap_ari = bootstrap_df['ari'].values
    perturbation_ari = perturbation_df['ari'].values

    summary = {
        'aa': args.aa,
        'bootstrap_mean_ari': np.mean(bootstrap_ari),
        'bootstrap_std_ari': np.std(bootstrap_ari),
        'bootstrap_min_ari': np.min(bootstrap_ari),
        'bootstrap_max_ari': np.max(bootstrap_ari),
        'perturbation_mean_ari': np.mean(perturbation_ari),
        'perturbation_std_ari': np.std(perturbation_ari),
        'perturbation_min_ari': np.min(perturbation_ari),
        'perturbation_max_ari': np.max(perturbation_ari),
    }

    # Save results
    out_dir.mkdir(parents=True, exist_ok=True)

    # Individual results
    agg_file = out_dir / f"stability_results_{args.aa}_aggregated.csv"
    results_df.to_csv(agg_file, index=False)

    # Summary stats
    summary_file = out_dir / f"stability_summary_{args.aa}.csv"
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(summary_file, index=False)

    logger.info(f"\n✓ Aggregated {len(results_df)} results")
    logger.info(f"Saved: {agg_file}")
    logger.info(f"Summary: {summary_file}")

    # Display summary
    logger.info("\n" + "=" * 80)
    logger.info("STABILITY VALIDATION SUMMARY")
    logger.info("=" * 80)

    logger.info(f"\nBootstrap Resampling ({len(bootstrap_ari)} iterations):")
    logger.info(f"  Mean ARI: {summary['bootstrap_mean_ari']:.3f} ± {summary['bootstrap_std_ari']:.3f}")
    logger.info(f"  Range: {summary['bootstrap_min_ari']:.3f} - {summary['bootstrap_max_ari']:.3f}")

    if summary['bootstrap_mean_ari'] > 0.85:
        logger.info(f"  ✓ STABLE (ARI >0.85)")
    else:
        logger.warning(f"  ⚠ UNSTABLE (ARI <0.85)")

    logger.info(f"\nPerturbation Analysis ({len(perturbation_ari)} iterations, noise=0.01):")
    logger.info(f"  Mean ARI: {summary['perturbation_mean_ari']:.3f} ± {summary['perturbation_std_ari']:.3f}")
    logger.info(f"  Range: {summary['perturbation_min_ari']:.3f} - {summary['perturbation_max_ari']:.3f}")

    if summary['perturbation_mean_ari'] > 0.85:
        logger.info(f"  ✓ ROBUST (ARI >0.85)")
    else:
        logger.warning(f"  ⚠ SENSITIVE (ARI <0.85)")

    logger.info("\n" + "=" * 80)

    if summary['bootstrap_mean_ari'] > 0.85 and summary['perturbation_mean_ari'] > 0.85:
        logger.info("✓ CLUSTERS ARE STABLE AND ROBUST - READY FOR PUBLICATION")
    else:
        logger.warning("⚠ Clusters may have stability issues - review results")

    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
