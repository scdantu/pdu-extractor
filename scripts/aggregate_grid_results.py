#!/usr/bin/env python3
"""
Aggregate parallel grid search results into single report.

Reads all task_XX_*.csv files and creates summary.

Usage:
  python3 aggregate_grid_results.py \
    --grid-dir analysis/grid_search_parallel \
    --aa L
"""

import argparse
from pathlib import Path
import logging

import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Aggregate grid search results")
    parser.add_argument("--grid-dir", required=True, help="Grid search results directory")
    parser.add_argument("--aa", required=True, help="Amino acid")
    parser.add_argument("--out-dir", default=None, help="Output directory (default: same as grid-dir)")

    args = parser.parse_args()

    grid_dir = Path(args.grid_dir)
    out_dir = Path(args.out_dir) if args.out_dir else grid_dir

    logger.info("=" * 80)
    logger.info("Aggregating Grid Search Results")
    logger.info("=" * 80)
    logger.info(f"Input: {grid_dir}")
    logger.info(f"Output: {out_dir}")

    # Find all task files
    task_files = sorted(grid_dir.glob("task_*.csv"))
    logger.info(f"Found {len(task_files)} task results")

    if not task_files:
        logger.error(f"No task files found in {grid_dir}")
        return

    # Aggregate results
    all_results = []
    for task_file in task_files:
        df = pd.read_csv(task_file)
        all_results.append(df)

    results_df = pd.concat(all_results, ignore_index=True)

    # Sort by silhouette score (descending)
    results_df = results_df.sort_values('silhouette', ascending=False)

    # Save aggregated results
    out_dir.mkdir(parents=True, exist_ok=True)
    agg_file = out_dir / f"grid_search_results_{args.aa}_aggregated.csv"
    results_df.to_csv(agg_file, index=False)

    logger.info(f"\n✓ Aggregated {len(results_df)} results")
    logger.info(f"Saved: {agg_file}")

    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("TOP 10 CONFIGURATIONS")
    logger.info("=" * 80)
    logger.info(f"{'Rank':<5} {'latent_dim':<12} {'min_cluster_size':<18} {'silhouette':<12} {'n_clusters':<12} {'noise_pct':<10}")
    logger.info("-" * 80)

    for idx, (_, row) in enumerate(results_df.head(10).iterrows(), 1):
        logger.info(
            f"{idx:<5} {int(row['latent_dim']):<12} {int(row['min_cluster_size']):<18} "
            f"{row['silhouette']:<12.3f} {int(row['n_clusters']):<12} {row['noise_pct']:<10.1f}%"
        )

    logger.info("=" * 80)

    # Check if current settings are in top 5
    top_5 = results_df.head(5)
    current_config = top_5[
        (top_5['latent_dim'] == 16) & (top_5['min_cluster_size'] == 200)
    ]

    if len(current_config) > 0:
        rank = list(results_df.index).index(current_config.index[0]) + 1
        logger.info(f"\n✓ Current settings (latent_dim=16, min_cluster_size=200)")
        logger.info(f"  Rank: {rank} / {len(results_df)}")
        logger.info(f"  Silhouette: {current_config.iloc[0]['silhouette']:.3f}")
    else:
        logger.warning(f"\n⚠ Current settings not in results")

    logger.info("=" * 80 + "\n")


if __name__ == "__main__":
    main()
