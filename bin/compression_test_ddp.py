#!/usr/bin/env python3
"""Compression test with Distributed Data Parallel (DDP) training.

This script uses PyTorch's DDP to parallelize training across multiple CPUs/GPUs.
Run with: torchrun --nproc_per_node=<N> bin/compression_test_ddp.py [args]

Example on ARCHER2 (128 CPUs):
    torchrun --nproc_per_node=128 \\
        bin/compression_test_ddp.py \\
        --aa L --encodings aa_20 --distances 15 --dims 16 --epochs 10

Expected speedup: ~128x (full utilization of all CPUs)
Expected time for all configs: 30-45 minutes (vs 4+ hours without DDP)
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import List, Tuple, Dict
import numpy as np
import os

# Add repo to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from pdusearch.features import FeatureExtractor, Encoding
from pdusearch.models import train_autoencoder_ddp, TrainingConfig
from pdusearch.clustering import HDBSCANClusterer, ClusterMetrics
from pdusearch.config import Config
try:
    from pdusearch.logging_utils import configure_logging, log_separator
except ImportError:
    # Fallback if logging_utils doesn't exist
    import logging
    def configure_logging(log_file=None, log_level="INFO"):
        logging.basicConfig(
            level=getattr(logging, log_level),
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            filename=log_file if log_file else None
        )
    def log_separator(logger):
        logger.info("=" * 70)

logger = logging.getLogger(__name__)


class CompressionTestDDP:
    """Compression test runner with DDP support."""

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
        """Initialize DDP test runner."""
        self.aa = aa.upper()
        self.encodings = encodings
        self.distances = distances
        self.latent_dims = latent_dims
        self.out_dir = Path(out_dir)
        self.db_dir = Path(db_dir)
        self.epochs = epochs
        self.verbose = verbose

        self.out_dir.mkdir(parents=True, exist_ok=True)

        self.encoding_map = {
            "aa_20": Encoding.AA_20,
            "functional_5": Encoding.FUNCTIONAL_5,
            "aa_20_conservation": Encoding.AA_20_WITH_CONSERVATION,
        }

        # Only rank 0 logs
        rank = int(os.environ.get("RANK", 0))
        if verbose and rank == 0:
            logger.info(
                f"CompressionTestDDP initialized for {self.aa}: "
                f"{len(encodings)} encodings × {len(distances)} distances × "
                f"{len(latent_dims)} dims"
            )

    def run(self) -> Dict:
        """Run compression test with DDP."""
        rank = int(os.environ.get("RANK", 0))
        world_size = int(os.environ.get("WORLD_SIZE", 1))

        results = {}

        for encoding_name in self.encodings:
            if encoding_name not in self.encoding_map:
                if rank == 0:
                    logger.warning(f"Unknown encoding: {encoding_name}")
                continue

            encoding = self.encoding_map[encoding_name]

            if rank == 0 and self.verbose:
                log_separator(logger)
                logger.info(f"Testing encoding: {encoding_name} ({encoding.value})")
                log_separator(logger)

            # Extract features once (only on rank 0)
            if rank == 0:
                try:
                    X_dict = self._extract_features(encoding)
                except Exception as e:
                    logger.error(f"Feature extraction failed: {e}")
                    continue
            else:
                X_dict = {}  # Non-rank-0 don't extract

            # Broadcast X_dict to all ranks if needed (for now, extract on each)
            # For simplicity, each rank extracts independently
            if rank == 0:
                X_dict = self._extract_features(encoding)

            test_results = []
            for distance in self.distances:
                for latent_dim in self.latent_dims:
                    try:
                        result = self._test_combination_ddp(
                            encoding_name, encoding, distance, latent_dim, X_dict
                        )
                        if result and rank == 0:
                            test_results.append(result)
                    except Exception as e:
                        if rank == 0:
                            logger.error(
                                f"Test failed: {encoding_name} {distance}Å {latent_dim}D: {e}"
                            )

            # Save results (rank 0 only)
            if test_results and rank == 0:
                self._save_results(encoding_name, test_results)
                results[encoding_name] = test_results

        return results

    def _extract_features(self, encoding: Encoding) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
        """Extract features for all test distances."""
        rank = int(os.environ.get("RANK", 0))

        if rank == 0 and self.verbose:
            logger.info(f"Extracting features for {self.aa}...")

        config = Config(
            db=str(self.db_dir / f"pdus_{self.aa}.sqlite"),
            radius=max(self.distances),
            aa=self.aa,
        )

        extractor = FeatureExtractor(config, encoding=encoding)
        X_full, pdu_ids = extractor.extract_for_aa(self.aa, max_distance=max(self.distances))

        if rank == 0 and self.verbose:
            logger.info(f"  → Extracted {X_full.shape[0]:,} PDUs, {X_full.shape[1]} features")

        return {dist: (X_full, pdu_ids) for dist in self.distances}

    def _test_combination_ddp(
        self,
        encoding_name: str,
        encoding: Encoding,
        distance: int,
        latent_dim: int,
        X_dict: Dict,
    ) -> Dict:
        """Test single combination using DDP."""
        rank = int(os.environ.get("RANK", 0))

        X, pdu_ids = X_dict[distance]

        if rank == 0 and self.verbose:
            logger.info(f"  Testing: {distance}Å, {latent_dim}D (DDP on all CPUs)")

        # Train with DDP
        train_cfg = TrainingConfig(
            latent_dim=latent_dim,
            epochs=self.epochs,
            batch_size=32,
            device="auto",
        )

        try:
            model, Z, losses = train_autoencoder_ddp(X, train_cfg, verbose=False)
        except Exception as e:
            if rank == 0:
                logger.error(f"Training failed: {e}")
            return None

        # Only rank 0 continues
        if rank != 0 or model is None:
            return None

        if self.verbose:
            logger.debug(f"    Training loss: {losses[-1]:.4f}")

        # Cluster
        clusterer = HDBSCANClusterer(min_cluster_size=200, min_samples=5)
        labels, confidences = clusterer.cluster(Z, verbose=False)

        # Metrics
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
        """Save results to CSV (rank 0 only)."""
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


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Compression hypothesis test with DDP",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("--aa", required=True, help="Amino acid to test")
    parser.add_argument(
        "--encodings",
        nargs="+",
        default=["aa_20"],
        help="Encodings to test",
    )
    parser.add_argument(
        "--distances",
        type=int,
        nargs="+",
        default=[10, 15],
        help="Distance cutoffs",
    )
    parser.add_argument(
        "--dims",
        type=int,
        nargs="+",
        default=[8, 16, 32, 64, 128, 256],
        help="Latent dimensions",
    )
    parser.add_argument("--epochs", type=int, default=10, help="Training epochs")
    parser.add_argument("--out-dir", default="analysis/compression_test", help="Output dir")
    parser.add_argument("--db-dir", default="per_aa_sqlite", help="Database dir")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument("--log-file", default=None, help="Log file")

    args = parser.parse_args()

    # Setup logging (rank 0 only)
    rank = int(os.environ.get("RANK", 0))
    if rank == 0:
        configure_logging(args.log_file, args.log_level)
        logger.info("=" * 70)
        logger.info("Compression Hypothesis Test - DDP")
        logger.info("=" * 70)

    # Run test
    runner = CompressionTestDDP(
        aa=args.aa,
        encodings=args.encodings,
        distances=args.distances,
        latent_dims=args.dims,
        out_dir=args.out_dir,
        db_dir=args.db_dir,
        epochs=args.epochs,
        verbose=(rank == 0),
    )

    results = runner.run()

    if rank == 0:
        logger.info("=" * 70)
        logger.info(f"✓ Test complete: {len(results)} encodings tested")
        logger.info(f"  Output: {runner.out_dir}")
        logger.info("=" * 70)


if __name__ == "__main__":
    main()
