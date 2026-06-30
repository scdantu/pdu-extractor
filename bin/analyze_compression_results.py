#!/usr/bin/env python3
"""Analyze compression hypothesis test results.

Reads results from compression test CSV files and determines:
  1. Which encoding gives best clustering (lowest noise)?
  2. Which distance cutoff matters (10Å vs 15Å)?
  3. Which latent dimension suffices?
  4. Overall conclusions for optimal pipeline

Usage:
    python bin/analyze_compression_results.py --aa L --out-dir analysis/compression_test
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List
import numpy as np
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from pdusearch.logging_utils import configure_logging, log_separator

logger = logging.getLogger(__name__)


def load_results(aa: str, out_dir: Path) -> Dict[str, any]:
    """Load compression test results for given AA.

    Args:
        aa: Amino acid
        out_dir: Output directory

    Returns:
        Dictionary mapping encoding_name -> dataframe
    """
    try:
        import pandas as pd
    except ImportError:
        logger.error("pandas required for analysis")
        return {}

    results = {}
    for csv_file in out_dir.glob(f"summary_{aa}_*.csv"):
        encoding = csv_file.stem.replace(f"summary_{aa}_", "")
        df = pd.read_csv(csv_file)
        results[encoding] = df
        logger.info(f"Loaded: {csv_file.name} ({len(df)} tests)")

    return results


def analyze_by_encoding(results: Dict) -> None:
    """Analyze results grouped by encoding.

    Args:
        results: Dictionary mapping encoding_name -> dataframe
    """
    log_separator(logger)
    logger.info("ANALYSIS BY ENCODING")
    log_separator(logger)

    for encoding, df in sorted(results.items()):
        logger.info(f"\n{encoding.upper()}:")
        logger.info(f"  Feature dimension: {df['input_features'].iloc[0]}")

        # Best configuration (lowest noise)
        best_idx = df["noise_percent"].idxmin()
        best = df.loc[best_idx]

        logger.info(f"  Best config:")
        logger.info(f"    Distance:  {best['distance']}Å")
        logger.info(f"    LatentDim: {best['latent_dim']}D")
        logger.info(f"    Noise:     {best['noise_percent']:.1f}%")
        logger.info(f"    Clusters:  {best['n_clusters']:.0f}")
        logger.info(f"    Confidence: {best['avg_confidence']:.3f}")

        # Statistics
        logger.info(f"  Noise statistics:")
        logger.info(f"    Min:  {df['noise_percent'].min():.1f}%")
        logger.info(f"    Max:  {df['noise_percent'].max():.1f}%")
        logger.info(f"    Mean: {df['noise_percent'].mean():.1f}%")
        logger.info(f"    Std:  {df['noise_percent'].std():.1f}%")


def analyze_distance_effect(results: Dict) -> None:
    """Analyze effect of distance cutoff.

    Args:
        results: Dictionary mapping encoding_name -> dataframe
    """
    log_separator(logger)
    logger.info("ANALYSIS: DISTANCE CUTOFF EFFECT")
    log_separator(logger)

    for encoding, df in sorted(results.items()):
        logger.info(f"\n{encoding}:")

        for dist in sorted(df["distance"].unique()):
            dist_data = df[df["distance"] == dist]
            logger.info(
                f"  {dist}Å: noise={dist_data['noise_percent'].mean():.1f}% "
                f"(min={dist_data['noise_percent'].min():.1f}%, "
                f"max={dist_data['noise_percent'].max():.1f}%)"
            )

        # Compare distances
        noise_10 = df[df["distance"] == 10]["noise_percent"].mean()
        noise_15 = df[df["distance"] == 15]["noise_percent"].mean()
        diff_pct = 100 * (noise_10 - noise_15) / noise_15

        if abs(diff_pct) > 5:
            better = "10Å" if noise_10 < noise_15 else "15Å"
            logger.info(f"  → {better} is {abs(diff_pct):.1f}% better")
        else:
            logger.info(f"  → Distance has minimal effect ({diff_pct:.1f}% difference)")


def analyze_compression_effect(results: Dict) -> None:
    """Analyze effect of latent dimension (compression).

    Args:
        results: Dictionary mapping encoding_name -> dataframe
    """
    log_separator(logger)
    logger.info("ANALYSIS: COMPRESSION EFFECT (Latent Dimension)")
    log_separator(logger)

    for encoding, df in sorted(results.items()):
        logger.info(f"\n{encoding}:")

        # Group by latent dimension
        for dim in sorted(df["latent_dim"].unique()):
            dim_data = df[df["latent_dim"] == dim]
            logger.info(
                f"  {dim:3d}D: noise={dim_data['noise_percent'].mean():.1f}% "
                f"(conf={dim_data['avg_confidence'].mean():.3f})"
            )

        # Find elbow: where does increasing dimension help diminish?
        dims = sorted(df["latent_dim"].unique())
        noise_by_dim = [df[df["latent_dim"] == d]["noise_percent"].mean() for d in dims]

        improvements = [noise_by_dim[0] - noise_by_dim[i] for i in range(len(noise_by_dim))]
        try:
            # Find where improvement < 5% relative
            for i in range(1, len(improvements)):
                if improvements[i] < improvements[i - 1] * 0.5:
                    logger.info(
                        f"  → Elbow at {dims[i]}D (diminishing returns after {dims[i-1]}D)"
                    )
                    break
        except:
            pass


def analyze_encoding_comparison(results: Dict) -> None:
    """Compare encodings directly.

    Args:
        results: Dictionary mapping encoding_name -> dataframe
    """
    log_separator(logger)
    logger.info("ENCODING COMPARISON")
    log_separator(logger)

    # Find best config for each encoding
    best_by_encoding = {}
    for encoding, df in results.items():
        best_idx = df["noise_percent"].idxmin()
        best_by_encoding[encoding] = df.loc[best_idx]

    # Table
    logger.info(f"\n{'Encoding':<25} {'Noise %':<10} {'Config':<30} {'Conf':<8}")
    logger.info("-" * 75)

    for encoding in sorted(best_by_encoding.keys()):
        best = best_by_encoding[encoding]
        config = f"{best['distance']:.0f}Å, {best['latent_dim']:.0f}D"
        logger.info(
            f"{encoding:<25} {best['noise_percent']:<10.1f} {config:<30} "
            f"{best['avg_confidence']:<8.3f}"
        )

    # Conclusions
    logger.info("\nKEY FINDINGS:")

    # Noise comparison
    noise_values = {enc: best_by_encoding[enc]["noise_percent"] for enc in best_by_encoding}
    best_enc = min(noise_values, key=noise_values.get)
    worst_enc = max(noise_values, key=noise_values.get)
    improvement = noise_values[worst_enc] - noise_values[best_enc]

    logger.info(f"  • Best encoding: {best_enc} ({noise_values[best_enc]:.1f}% noise)")
    logger.info(f"  • Difference: {improvement:.1f}% (encoding choice matters: {improvement > 5})")

    # Functional vs AA20
    if "functional_5" in noise_values and "aa_20" in noise_values:
        func_vs_aa = noise_values["functional_5"] - noise_values["aa_20"]
        logger.info(
            f"  • Functional-5 vs AA-20: {func_vs_aa:+.1f}% "
            f"({'worse' if func_vs_aa > 0 else 'better'})"
        )

    # Conservation signal
    if "aa_20_conservation" in noise_values:
        cons_vs_aa = noise_values["aa_20_conservation"] - noise_values["aa_20"]
        logger.info(
            f"  • Conservation signal: {cons_vs_aa:+.1f}% "
            f"({'not helpful' if cons_vs_aa > 0 else 'helpful'})"
        )


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Analyze compression hypothesis test results")
    parser.add_argument("--aa", required=True, help="Amino acid tested (L, A, G)")
    parser.add_argument("--out-dir", default="analysis/compression_test", help="Results directory")
    parser.add_argument("--log-level", default="INFO", help="Logging level")

    args = parser.parse_args()

    # Setup logging
    configure_logging(args.log_level)

    out_dir = Path(args.out_dir)
    if not out_dir.exists():
        logger.error(f"Output directory not found: {out_dir}")
        sys.exit(1)

    # Load results
    results = load_results(args.aa, out_dir)
    if not results:
        logger.error("No results found")
        sys.exit(1)

    # Analyze
    log_separator(logger)
    logger.info(f"COMPRESSION HYPOTHESIS TEST ANALYSIS: {args.aa}")
    log_separator(logger)

    analyze_by_encoding(results)
    analyze_distance_effect(results)
    analyze_compression_effect(results)
    analyze_encoding_comparison(results)

    log_separator(logger)
    logger.info("✓ Analysis complete")


if __name__ == "__main__":
    main()
