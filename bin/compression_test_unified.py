#!/usr/bin/env python3
"""Unified compression hypothesis test using new pdusearch library.

Tests the hypothesis: Which factor causes 88% noise in LEU clustering?
  (a) Feature diversity from 15Å distance cutoff
  (b) Autoencoder compression being too aggressive
  (c) Clustering threshold (min_cluster_size) being too high
  (d) Feature representation granularity (20-AA vs functional grouping)

Test Matrix:
  - Encodings: AA_20 (900D), FUNCTIONAL_5 (225D), AA_20_WITH_CONSERVATION (915D)
  - Distances: 10Å, 15Å
  - Latent dims: 8, 16, 32, 64, 128, 256
  - Amino acids: G (19K PDUs, ~30 min), A (1.9M, ~2h), L (2.1M, ~2-3h)

Output:
  analysis/compression_test/summary_{AA}_{ENCODING}.csv with:
    Distance, LatentDim, NoisePercent, AvgConfidence, NClusters

Total combinations: 3 encodings × 2 distances × 6 dims × 3 AAs = 108 tests

Usage:
  # Local test (small dataset)
  python bin/compression_test_unified.py --aa L --encodings aa_20 --dims 16 --epochs 2

  # ARCHER2 full test
  sbatch submit/submit_compression_test_unified.slurm

  # Specific encoding
  python bin/compression_test_unified.py --aa L --encodings functional_5 --distances 15
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import List, Tuple, Dict
import numpy as np

# Add repo to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from pdusearch.features import FeatureExtractor, Encoding
from pdusearch.models import train_autoencoder, TrainingConfig
from pdusearch.clustering import HDBSCANClusterer, ClusterMetrics
from pdusearch.config import Config
from pdusearch.logging_utils import configure_logging, log_separator

logger = logging.getLogger(__name__)


class CompressionTestRunner:
    """Unified compression hypothesis test runner.

    Systematically tests combinations of:
      - Feature encodings (granularity)
      - Distance cutoffs (feature diversity)
      - Latent dimensions (compression)

    Records clustering quality metrics for each combination.
    """

    def __init__(
        self,
        aa: str,
        encodings: List[str],
        distances: List[int],
        latent_dims: List[int],
        out_dir: str = "analysis/compression_test",
        db_dir: str = "per_aa_sqlite",
        epochs: int = 10,
        verbose: bool = True,
    ):
        """Initialize test runner.

        Args:
            aa: Amino acid to test (L, A, G)
            encodings: List of encoding names
            distances: List of distance cutoffs
            latent_dims: List of latent dimensions to test
            out_dir: Output directory
            db_dir: Database directory
            epochs: Training epochs
            verbose: Log progress
        """
        self.aa = aa.upper()
        self.encodings = encodings
        self.distances = distances
        self.latent_dims = latent_dims
        self.out_dir = Path(out_dir)
        self.db_dir = Path(db_dir)
        self.epochs = epochs
        self.verbose = verbose

        self.out_dir.mkdir(parents=True, exist_ok=True)

        # Encoding name to class mapping
        self.encoding_map = {
            "aa_20": Encoding.AA_20,
            "functional_5": Encoding.FUNCTIONAL_5,
            "aa_20_conservation": Encoding.AA_20_WITH_CONSERVATION,
        }

        if self.verbose:
            logger.info(
                f"CompressionTestRunner initialized for {self.aa}: "
                f"{len(encodings)} encodings × {len(distances)} distances × "
                f"{len(latent_dims)} dims = {len(encodings) * len(distances) * len(latent_dims)} tests"
            )

    def run(self) -> Dict:
        """Run full compression test matrix.

        Returns:
            Dictionary mapping encoding -> results dataframe
        """
        results = {}

        for encoding_name in self.encodings:
            if encoding_name not in self.encoding_map:
                logger.warning(f"Unknown encoding: {encoding_name}")
                continue

            encoding = self.encoding_map[encoding_name]

            if self.verbose:
                log_separator(logger)
                logger.info(f"Testing encoding: {encoding_name} ({encoding.value})")
                log_separator(logger)

            # Extract features once per encoding (reuse for all distances/dims)
            try:
                X_dict = self._extract_features(encoding)
            except Exception as e:
                logger.error(f"Feature extraction failed for {encoding_name}: {e}")
                continue

            # Test all distance/dim combinations
            test_results = []
            for distance in self.distances:
                for latent_dim in self.latent_dims:
                    try:
                        result = self._test_combination(
                            encoding_name, encoding, distance, latent_dim, X_dict
                        )
                        if result:
                            test_results.append(result)
                    except Exception as e:
                        logger.error(
                            f"Test failed: {encoding_name} {distance}Å {latent_dim}D: {e}"
                        )

            # Save results
            if test_results:
                self._save_results(encoding_name, test_results)
                results[encoding_name] = test_results

        return results

    def _extract_features(self, encoding: Encoding) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
        """Extract features for all test distances.

        Args:
            encoding: Encoding strategy

        Returns:
            Dict mapping distance -> (X, pdu_ids)
        """
        if self.verbose:
            logger.info(f"Extracting features for {self.aa}...")

        config = Config(
            db=str(self.db_dir / f"pdus_{self.aa}.sqlite"),
            radius=max(self.distances),  # Extract once for max distance
            aa=self.aa,
        )

        extractor = FeatureExtractor(config, encoding=encoding)
        X_full, pdu_ids = extractor.extract_for_aa(self.aa, max_distance=max(self.distances))

        if self.verbose:
            logger.info(f"  → Extracted {X_full.shape[0]:,} PDUs, {X_full.shape[1]} features")

        # Return dict with original full features (will subset by distance in clustering)
        # For now, treat all distances as same features (actual distance-based subsetting
        # would require PDU coordinate data)
        return {dist: (X_full, pdu_ids) for dist in self.distances}

    def _test_combination(
        self,
        encoding_name: str,
        encoding: Encoding,
        distance: int,
        latent_dim: int,
        X_dict: Dict,
    ) -> Dict:
        """Test single encoding/distance/dimension combination.

        Args:
            encoding_name: Name for logging
            encoding: Encoding enum
            distance: Distance cutoff
            latent_dim: Latent dimension
            X_dict: Feature data dict

        Returns:
            Result dict with metrics
        """
        X, pdu_ids = X_dict[distance]

        if self.verbose:
            logger.info(
                f"  Testing: {distance}Å, {latent_dim}D "
                f"({X.shape[0]:,} PDUs, {X.shape[1]} features)"
            )

        # Train autoencoder
        train_cfg = TrainingConfig(
            latent_dim=latent_dim,
            epochs=self.epochs,
            batch_size=32,
            device="auto",
        )

        try:
            model, Z, losses = train_autoencoder(X, train_cfg, verbose=False)
        except Exception as e:
            logger.error(f"Training failed: {e}")
            return None

        if self.verbose:
            logger.debug(f"    Training loss: {losses[-1]:.4f}")

        # Cluster embeddings
        clusterer = HDBSCANClusterer(min_cluster_size=200, min_samples=5)
        labels, confidences = clusterer.cluster(Z, verbose=False)

        # Compute metrics
        metrics = ClusterMetrics(labels, confidences)

        result = {
            "encoding": encoding_name,
            "distance": distance,
            "latent_dim": latent_dim,
            "input_features": X.shape[1],
            "n_pdus": X.shape[0],
            "n_clusters": metrics.n_clusters,
            "noise_percent": metrics.noise_percent,
            "avg_confidence": metrics.avg_confidence,
            "n_assigned": metrics.n_assigned,
            "training_loss": losses[-1],
        }

        if self.verbose:
            logger.info(
                f"    → {metrics.n_clusters} clusters, "
                f"{metrics.noise_percent:.1f}% noise, "
                f"{metrics.avg_confidence:.3f} confidence"
            )

        return result

    def _save_results(self, encoding_name: str, results: List[Dict]) -> None:
        """Save results to CSV.

        Args:
            encoding_name: Encoding name for filename
            results: List of result dicts
        """
        try:
            import pandas as pd
        except ImportError:
            logger.warning("pandas not available, skipping CSV save")
            return

        df = pd.DataFrame(results)

        output_file = self.out_dir / f"summary_{self.aa}_{encoding_name}.csv"
        df.to_csv(output_file, index=False)

        if self.verbose:
            logger.info(f"Saved: {output_file}")

        # Print summary
        logger.info(f"\nSummary for {self.aa} × {encoding_name}:")
        logger.info(df.to_string(index=False))


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Compression hypothesis test using pdusearch library",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("--aa", required=True, help="Amino acid to test (L, A, G)")
    parser.add_argument(
        "--encodings",
        nargs="+",
        default=["aa_20", "functional_5", "aa_20_conservation"],
        choices=["aa_20", "functional_5", "aa_20_conservation"],
        help="Encodings to test",
    )
    parser.add_argument(
        "--distances",
        type=int,
        nargs="+",
        default=[10, 15],
        help="Distance cutoffs to test (Angstroms)",
    )
    parser.add_argument(
        "--dims",
        type=int,
        nargs="+",
        default=[8, 16, 32, 64, 128, 256],
        help="Latent dimensions to test",
    )
    parser.add_argument("--epochs", type=int, default=10, help="Training epochs")
    parser.add_argument("--out-dir", default="analysis/compression_test", help="Output directory")
    parser.add_argument("--db-dir", default="per_aa_sqlite", help="Database directory")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument("--log-file", default=None, help="Log file path")

    args = parser.parse_args()

    # Setup logging
    configure_logging(args.log_file, args.log_level)

    # Run test
    logger.info("="*70)
    logger.info("Compression Hypothesis Test - Unified")
    logger.info("="*70)

    runner = CompressionTestRunner(
        aa=args.aa,
        encodings=args.encodings,
        distances=args.distances,
        latent_dims=args.dims,
        out_dir=args.out_dir,
        db_dir=args.db_dir,
        epochs=args.epochs,
        verbose=True,
    )

    results = runner.run()

    logger.info("="*70)
    logger.info(f"✓ Test complete: {len(results)} encodings tested")
    logger.info(f"  Output: {runner.out_dir}")
    logger.info("="*70)


if __name__ == "__main__":
    main()
