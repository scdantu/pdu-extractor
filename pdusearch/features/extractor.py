"""PDU feature extraction with configurable encoding strategies.

Main module for extracting fixed-length feature vectors from PDU database.
Supports multiple encoding strategies (AA_20, FUNCTIONAL_5, AA_20_WITH_CONSERVATION)
and provides both single-AA and batch extraction interfaces.

Example:
    >>> from pdusearch.features import FeatureExtractor, Encoding
    >>> from pdusearch.config import Config
    >>>
    >>> config = Config(db="per_aa_sqlite/pdus_L.sqlite", radius=15.0)
    >>> extractor = FeatureExtractor(config, encoding=Encoding.AA_20)
    >>> X, pdu_ids = extractor.extract_for_aa("L")
    >>> print(X.shape)  # (n_pdus, 900)
    >>>
    >>> # Extract with functional encoding (4x fewer features)
    >>> extractor_func = FeatureExtractor(config, encoding=Encoding.FUNCTIONAL_5)
    >>> X_func, ids_func = extractor_func.extract_for_aa("L")
    >>> print(X_func.shape)  # (n_pdus, 225)
"""

from pathlib import Path
from typing import Optional, List, Tuple, Dict
import numpy as np
import logging

from pdusearch.config import Config
from pdusearch.io import PDUDatabase
from pdusearch.logging_utils import get_logger
from pdusearch.utils import DEFAULT_RADIUS

from .encodings import Encoding, get_encoding_strategy, EncodingStrategy
from .alignment import AlignmentProvider, NullAlignmentProvider

logger = logging.getLogger(__name__)


class FeatureExtractor:
    """Extract PDU features with configurable encoding.

    Loads PDU neighbor data from SQLite database and encodes it as fixed-length
    feature vectors using the specified encoding strategy. Supports efficient
    batch extraction and multiple output formats.

    Attributes:
        config: Configuration object with database and radius settings
        encoding: Feature encoding strategy (AA_20, FUNCTIONAL_5, etc.)
        alignment_provider: Optional provider for conservation scores
    """

    def __init__(
        self,
        config: Config,
        encoding: Encoding = Encoding.AA_20,
        alignment_provider: Optional[AlignmentProvider] = None,
    ):
        """Initialize feature extractor.

        Args:
            config: Config object specifying database path, radius, etc.
            encoding: Encoding strategy to use (default: AA_20)
            alignment_provider: Optional provider for conservation metrics
                               (only used with AA_20_WITH_CONSERVATION encoding)

        Raises:
            FileNotFoundError: If database does not exist
            ValueError: If encoding is not recognized
        """
        self.config = config
        self.encoding = encoding
        self.alignment_provider = alignment_provider or NullAlignmentProvider()

        # Validate database exists
        db_path = Path(config.db)
        if not db_path.exists():
            raise FileNotFoundError(f"Database not found: {config.db}")

        # Initialize database connection (lazy loaded on first use)
        self.db = PDUDatabase(config.db)

        # Get encoding strategy
        try:
            self.strategy = get_encoding_strategy(encoding)
        except ValueError as e:
            raise ValueError(f"Invalid encoding: {e}")

        self.logger = get_logger(__name__)
        self.logger.debug(
            f"FeatureExtractor initialized: encoding={encoding.value}, "
            f"radius={config.radius}Å, db={config.db}"
        )

    def extract_for_aa(
        self,
        aa: str,
        max_distance: Optional[float] = None,
        save_path: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Extract features for all PDUs of a specific amino acid.

        Retrieves all PDUs with the specified central amino acid from the
        database, encodes their neighborhood composition as feature vectors,
        and optionally saves to disk.

        Args:
            aa: Amino acid one-letter code (e.g., 'L', 'A', 'G')
            max_distance: Override config radius for this extraction (Angstroms)
            save_path: Optional path to save features as NPZ file

        Returns:
            Tuple of (X, pdu_ids) where:
                - X: Feature matrix of shape (n_pdus, n_features) as float32
                - pdu_ids: Array of PDU IDs of shape (n_pdus,)

        Raises:
            ValueError: If amino acid is invalid or no PDUs found

        Example:
            >>> X, ids = extractor.extract_for_aa("L", save_path="features_L.npz")
            >>> print(f"Extracted {X.shape[0]} PDUs, {X.shape[1]} features")
        """
        # Normalize input
        aa = aa.upper()
        if len(aa) != 1:
            raise ValueError(f"Invalid amino acid: {aa}. Expected single character.")

        # Use config radius or override
        distance = max_distance if max_distance is not None else self.config.radius

        self.logger.info(
            f"Extracting {aa} features: {distance}Å radius, {self.strategy.get_encoding_name()} encoding"
        )

        # Load from database
        self.db.connect()
        try:
            # Get PDU IDs for this amino acid
            pdu_ids = self.db.get_pdu_ids_for_aa(aa)
            if not pdu_ids:
                self.logger.warning(f"No PDUs found for amino acid: {aa}")
                return np.array([]), np.array([])

            self.logger.info(f"  → Found {len(pdu_ids):,} PDUs")

            # Get neighbor data (residue type, SS, distance) for all PDUs of this AA
            neighbors_data = self.db.get_pdu_neighbors(aa, max_distance=distance)
            self.logger.info(f"  → Found {len(neighbors_data):,} neighbor relationships")

            # Organize neighbors by PDU ID for efficient lookup
            neighbors_by_pdu: Dict[int, List[Tuple]] = {pdu_id: [] for pdu_id in pdu_ids}
            for neighbor in neighbors_data:
                pdu_id = neighbor[0]
                if pdu_id in neighbors_by_pdu:
                    neighbors_by_pdu[pdu_id].append(neighbor[1:])

            # Calculate feature dimension
            n_distance_bins = self.config.n_distance_bins
            n_features = self.strategy.get_n_features(n_distance_bins)

            # Allocate feature matrix
            X = np.zeros((len(pdu_ids), n_features), dtype=np.float32)

            # Encode each PDU
            self.logger.info(f"  → Encoding {len(pdu_ids):,} PDUs...")
            for idx, pdu_id in enumerate(pdu_ids):
                neighbors = neighbors_by_pdu.get(pdu_id, [])

                # Encode neighbors (residue_type, ss, distance)
                # For conservation encoding, also pass conservation scores
                if self.encoding == Encoding.AA_20_WITH_CONSERVATION:
                    # Extract conservation scores if available
                    conservation_scores = self._get_conservation_scores(
                        aa, pdu_id, neighbors
                    )
                    X[idx] = self.strategy.encode_neighbors(
                        neighbors,
                        n_distance_bins,
                        radius=distance,
                        bin_width=self.config.bin_width,
                        conservation_scores=conservation_scores,
                    )
                else:
                    X[idx] = self.strategy.encode_neighbors(
                        neighbors,
                        n_distance_bins,
                        radius=distance,
                        bin_width=self.config.bin_width,
                    )

                # Log progress
                if (idx + 1) % max(1, len(pdu_ids) // 10) == 0:
                    self.logger.debug(f"    Encoded {idx + 1:,}/{len(pdu_ids):,} PDUs")

            # Convert to numpy array
            pdu_ids_array = np.array(pdu_ids, dtype=np.int64)

            # Save if requested
            if save_path:
                self.save_features(X, pdu_ids_array, save_path, aa=aa)

            self.logger.info(
                f"✓ Extracted {X.shape[0]:,} PDUs × {X.shape[1]} features "
                f"for {aa} ({self.strategy.get_encoding_name()})"
            )

            return X, pdu_ids_array

        finally:
            self.db.close()

    def extract_batch(
        self,
        aa_list: Optional[List[str]] = None,
        save_dir: Optional[str] = None,
    ) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """Extract features for multiple amino acids.

        Convenience method for batch extraction. Extracts and optionally saves
        features for each amino acid in the list.

        Args:
            aa_list: List of amino acids to extract. If None, uses all 20 standard AAs.
            save_dir: Optional directory to save NPZ files per AA

        Returns:
            Dictionary mapping amino acid -> (X, pdu_ids) tuples

        Example:
            >>> results = extractor.extract_batch(["L", "A", "G"])
            >>> for aa, (X, ids) in results.items():
            ...     print(f"{aa}: {X.shape}")
        """
        if aa_list is None:
            aa_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                      "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

        results = {}
        for aa in aa_list:
            try:
                save_path = None
                if save_dir:
                    Path(save_dir).mkdir(parents=True, exist_ok=True)
                    save_path = f"{save_dir}/pdu_features_{aa}.npz"

                X, ids = self.extract_for_aa(aa, save_path=save_path)
                if len(ids) > 0:
                    results[aa] = (X, ids)
                else:
                    self.logger.warning(f"No features extracted for {aa}")

            except Exception as e:
                self.logger.error(f"Failed to extract {aa}: {e}")

        self.logger.info(f"✓ Batch extraction complete: {len(results)} amino acids")
        return results

    def save_features(
        self,
        X: np.ndarray,
        pdu_ids: np.ndarray,
        output_path: str,
        aa: Optional[str] = None,
    ) -> None:
        """Save extracted features to NPZ file.

        Args:
            X: Feature matrix of shape (n_pdus, n_features)
            pdu_ids: Array of PDU IDs
            output_path: Path to output NPZ file
            aa: Optional amino acid label for metadata

        Example:
            >>> extractor.save_features(X, ids, "features_L.npz", aa="L")
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save with metadata
        save_dict = {
            "X": X,
            "pdu_ids": pdu_ids,
            "encoding": self.encoding.value,
            "radius": self.config.radius,
            "bin_width": self.config.bin_width,
        }

        if aa is not None:
            save_dict["amino_acid"] = aa

        np.savez(output_path, **save_dict)
        self.logger.info(f"Saved features: {output_path} ({X.nbytes / 1e6:.1f} MB)")

    @staticmethod
    def load_features(npz_path: str) -> Tuple[np.ndarray, np.ndarray, Dict]:
        """Load features from NPZ file.

        Args:
            npz_path: Path to NPZ file saved by save_features()

        Returns:
            Tuple of (X, pdu_ids, metadata)
                - X: Feature matrix
                - pdu_ids: PDU identifiers
                - metadata: Dictionary with encoding, radius, etc.

        Example:
            >>> X, ids, meta = FeatureExtractor.load_features("features_L.npz")
            >>> print(f"Encoding: {meta['encoding']}, Radius: {meta['radius']}")
        """
        data = np.load(npz_path)

        X = data["X"]
        pdu_ids = data["pdu_ids"]

        # Extract metadata
        metadata = {}
        for key in ["encoding", "radius", "bin_width", "amino_acid"]:
            if key in data:
                val = data[key]
                # Handle numpy scalar types
                if hasattr(val, "item"):
                    metadata[key] = val.item()
                else:
                    metadata[key] = val

        logger.info(f"Loaded {X.shape[0]:,} PDUs × {X.shape[1]} features from {npz_path}")

        return X, pdu_ids, metadata

    def _get_conservation_scores(
        self,
        aa: str,
        pdu_id: int,
        neighbors: List[Tuple[str, str, float]],
    ) -> np.ndarray:
        """Get conservation scores for neighbors (placeholder).

        This is a placeholder that returns zeros. Will be populated by PyCoM
        integration once alignment provider is connected to PDB metadata.

        Args:
            aa: Central amino acid
            pdu_id: PDU identifier
            neighbors: List of neighbor tuples (residue_type, ss, distance)

        Returns:
            Array of conservation scores, one per neighbor
        """
        # Placeholder implementation
        if self.alignment_provider is None:
            return np.zeros(len(neighbors), dtype=np.float32)

        # This would be populated with actual PyCoM queries in future
        conservation_scores = np.zeros(len(neighbors), dtype=np.float32)
        return conservation_scores

    def get_feature_names(self) -> List[str]:
        """Get human-readable feature names.

        Returns:
            List of feature names for interpretation and visualization
        """
        n_distance_bins = self.config.n_distance_bins
        bin_width = self.config.bin_width

        feature_names = []

        if self.encoding == Encoding.AA_20:
            # 20 AA types × 3 SS types × distance bins
            aa_names = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                       "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
            ss_names = ["H", "E", "C"]

            for aa in aa_names:
                for ss in ss_names:
                    for bin_idx in range(n_distance_bins):
                        dist_start = bin_idx * bin_width
                        dist_end = (bin_idx + 1) * bin_width
                        feature_names.append(f"{aa}_{ss}_{dist_start:.1f}-{dist_end:.1f}Å")

        elif self.encoding == Encoding.FUNCTIONAL_5:
            # 5 functional classes × 3 SS types × distance bins
            func_names = ["HY", "PO", "PO+", "PO-", "SP"]
            ss_names = ["H", "E", "C"]

            for func in func_names:
                for ss in ss_names:
                    for bin_idx in range(n_distance_bins):
                        dist_start = bin_idx * bin_width
                        dist_end = (bin_idx + 1) * bin_width
                        feature_names.append(f"{func}_{ss}_{dist_start:.1f}-{dist_end:.1f}Å")

        elif self.encoding == Encoding.AA_20_WITH_CONSERVATION:
            # AA20 features + entropy per bin
            # (same as AA_20 for first 900 features)
            aa_names = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                       "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
            ss_names = ["H", "E", "C"]

            for aa in aa_names:
                for ss in ss_names:
                    for bin_idx in range(n_distance_bins):
                        dist_start = bin_idx * bin_width
                        dist_end = (bin_idx + 1) * bin_width
                        feature_names.append(f"{aa}_{ss}_{dist_start:.1f}-{dist_end:.1f}Å")

            # Conservation scores per bin
            for bin_idx in range(n_distance_bins):
                dist_start = bin_idx * bin_width
                dist_end = (bin_idx + 1) * bin_width
                feature_names.append(f"conservation_{dist_start:.1f}-{dist_end:.1f}Å")

        return feature_names

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"FeatureExtractor("
            f"db={self.config.db}, "
            f"encoding={self.encoding.value}, "
            f"radius={self.config.radius}Å"
            f")"
        )
