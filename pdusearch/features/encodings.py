"""Feature encoding strategies for PDU representation.

This module provides multiple ways to encode PDU neighbor information into fixed-length
feature vectors. Each encoding strategy represents neighborhood composition differently,
enabling hypothesis testing on feature representation choices.

Encoding Types:
    - AA_20: Standard 20 amino acid types × 3 secondary structures × distance bins
    - FUNCTIONAL_5: Functional residue classes (hydrophobic, polar, charged, special)
    - AA_20_WITH_CONSERVATION: 20 AA types + Shannon entropy conservation scores

Example:
    >>> from pdusearch.features.encodings import Encoding, AA20Encoding
    >>> strategy = AA20Encoding()
    >>> n_features = strategy.get_n_features(n_distance_bins=15)
    >>> print(n_features)  # 900
"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import List, Tuple, Optional, Dict
import numpy as np
import logging

from pdusearch.utils import (
    AA_ORDER,
    SS_ORDER,
    AA_TO_FUNCTIONAL,
    FunctionalClass,
    FUNCTIONAL_CLASS_ORDER,
)

logger = logging.getLogger(__name__)


class Encoding(Enum):
    """Feature encoding strategies for PDUs.

    Attributes:
        AA_20: Standard 20 amino acid types (current method)
        FUNCTIONAL_5: 5 functional residue classes (hydrophobic, polar, charged+/-, special)
        AA_20_WITH_CONSERVATION: 20 AA types + Shannon entropy per distance shell
    """
    AA_20 = "aa_20"
    FUNCTIONAL_5 = "functional_5"
    AA_20_WITH_CONSERVATION = "aa_20_conservation"


class EncodingStrategy(ABC):
    """Abstract base class for PDU feature encoding strategies.

    Defines the interface that all encoding strategies must implement.
    Each strategy converts a list of neighbor residues into a fixed-length
    feature vector representing the neighborhood composition.
    """

    @abstractmethod
    def get_n_features(self, n_distance_bins: int) -> int:
        """Calculate total feature vector dimension.

        Args:
            n_distance_bins: Number of distance shells (radius / bin_width)

        Returns:
            Total number of features in encoded vector
        """
        pass

    @abstractmethod
    def encode_neighbors(
        self,
        neighbors: List[Tuple[str, str, float]],
        n_distance_bins: int,
    ) -> np.ndarray:
        """Encode neighbor residues into feature vector.

        Args:
            neighbors: List of (residue_type, secondary_structure, distance_angstrom) tuples
            n_distance_bins: Number of distance bins for quantization

        Returns:
            Feature vector of shape (n_features,) as float32
        """
        pass

    @abstractmethod
    def get_encoding_name(self) -> str:
        """Get human-readable name of encoding strategy.

        Returns:
            Name suitable for logging and filenames
        """
        pass


class AA20Encoding(EncodingStrategy):
    """Standard 20-amino-acid encoding.

    Encodes neighborhood as: 20 amino acids × 3 secondary structures × distance bins
    Total features: 20 * 3 * n_distance_bins = 900 features (for 15 distance bins)

    This is the baseline encoding used in the current PDU extraction pipeline.
    Each dimension counts how many residues of type (AA, SS) are in each distance shell.

    Example:
        >>> encoding = AA20Encoding()
        >>> neighbors = [("L", "H", 3.5), ("A", "H", 5.2), ("V", "E", 8.1)]
        >>> features = encoding.encode_neighbors(neighbors, n_distance_bins=15)
        >>> features.shape
        (900,)
    """

    def __init__(self):
        """Initialize AA20 encoding strategy."""
        self.aa_list = list(AA_ORDER)
        self.ss_list = list(SS_ORDER)
        self.logger = logging.getLogger(__name__)

    def get_n_features(self, n_distance_bins: int) -> int:
        """Calculate feature dimension.

        Args:
            n_distance_bins: Number of distance shells

        Returns:
            20 * 3 * n_distance_bins features
        """
        return len(self.aa_list) * len(self.ss_list) * n_distance_bins

    def encode_neighbors(
        self,
        neighbors: List[Tuple[str, str, float]],
        n_distance_bins: int,
        radius: float = 15.0,
        bin_width: float = 1.0,
    ) -> np.ndarray:
        """Encode neighbors using amino acid identity.

        Creates a vector where each dimension represents the count of residues
        with specific (amino_acid_type, secondary_structure, distance_shell) properties.

        Args:
            neighbors: List of (residue_type, ss, distance) tuples
            n_distance_bins: Total number of distance bins
            radius: Maximum distance cutoff (Angstroms)
            bin_width: Width of each distance bin (Angstroms)

        Returns:
            Feature vector of shape (n_features,) as float32

        Details:
            - Distance bins are 1Å intervals: [0-1], [1-2], ..., [14-15]
            - Missing neighbors contribute 0 to feature vector
            - Neighbors beyond radius are ignored
        """
        n_features = self.get_n_features(n_distance_bins)
        features = np.zeros(n_features, dtype=np.float32)

        # Build lookup tables for fast indexing
        # Feature index = aa_idx * 3 * n_distance_bins + ss_idx * n_distance_bins + bin_idx
        aa_to_idx = {aa: i for i, aa in enumerate(self.aa_list)}
        ss_to_idx = {ss: i for i, ss in enumerate(self.ss_list)}

        # Process each neighbor
        for residue_type, secondary_structure, distance in neighbors:
            # Validate inputs
            if residue_type not in aa_to_idx:
                # Skip non-standard residues (X, etc.)
                continue
            if secondary_structure not in ss_to_idx:
                # Skip if secondary structure is not H/E/C
                continue
            if distance > radius:
                # Skip neighbors beyond radius
                continue

            # Calculate distance bin
            distance_bin = int(distance / bin_width)
            if distance_bin >= n_distance_bins:
                distance_bin = n_distance_bins - 1  # Clamp to last bin

            # Calculate feature index
            aa_idx = aa_to_idx[residue_type]
            ss_idx = ss_to_idx[secondary_structure]
            feature_idx = aa_idx * len(self.ss_list) * n_distance_bins + ss_idx * n_distance_bins + distance_bin

            # Increment feature (count residues in each category)
            features[feature_idx] += 1.0

        return features

    def get_encoding_name(self) -> str:
        """Get name of this encoding strategy."""
        return "AA_20"


class Functional5Encoding(EncodingStrategy):
    """Functional residue class encoding.

    Groups 20 amino acids into 5 functional classes:
        - HYDROPHOBIC: A, V, L, I, M, F, W, P
        - POLAR: S, T, N, Q, C
        - POSITIVE: K, R, H
        - NEGATIVE: D, E
        - SPECIAL: G, Y

    Encodes as: 5 functional classes × 3 secondary structures × distance bins
    Total features: 5 * 3 * n_distance_bins = 225 features (for 15 distance bins)

    This encoding tests the hypothesis that amino acid *functional properties*
    matter more than exact identity. Reduces feature dimensionality 4x while
    potentially capturing essential neighborhood properties.

    Example:
        >>> encoding = Functional5Encoding()
        >>> neighbors = [("L", "H", 3.5), ("A", "H", 5.2), ("D", "E", 8.1)]
        >>> features = encoding.encode_neighbors(neighbors, n_distance_bins=15)
        >>> features.shape
        (225,)
    """

    def __init__(self):
        """Initialize functional encoding strategy."""
        self.functional_classes = [cls for cls in FunctionalClass]
        self.ss_list = list(SS_ORDER)
        self.logger = logging.getLogger(__name__)

    def get_n_features(self, n_distance_bins: int) -> int:
        """Calculate feature dimension.

        Args:
            n_distance_bins: Number of distance shells

        Returns:
            5 * 3 * n_distance_bins features
        """
        return len(self.functional_classes) * len(self.ss_list) * n_distance_bins

    def encode_neighbors(
        self,
        neighbors: List[Tuple[str, str, float]],
        n_distance_bins: int,
        radius: float = 15.0,
        bin_width: float = 1.0,
    ) -> np.ndarray:
        """Encode neighbors using functional residue classes.

        Groups amino acids by physicochemical properties and encodes composition.

        Args:
            neighbors: List of (residue_type, ss, distance) tuples
            n_distance_bins: Total number of distance bins
            radius: Maximum distance cutoff (Angstroms)
            bin_width: Width of each distance bin (Angstroms)

        Returns:
            Feature vector of shape (n_features,) as float32

        Details:
            - Same distance binning as AA20Encoding
            - Non-standard residues (X, etc.) contribute 0
            - Reduces dimensionality: 900D → 225D (4x compression)
            - Interpretation: "What functional types surround this residue?"
        """
        n_features = self.get_n_features(n_distance_bins)
        features = np.zeros(n_features, dtype=np.float32)

        # Build lookup tables
        func_to_idx = {fc: i for i, fc in enumerate(self.functional_classes)}
        ss_to_idx = {ss: i for i, ss in enumerate(self.ss_list)}

        # Process each neighbor
        for residue_type, secondary_structure, distance in neighbors:
            # Validate inputs
            if residue_type not in AA_TO_FUNCTIONAL:
                # Skip non-standard residues
                continue
            if secondary_structure not in ss_to_idx:
                continue
            if distance > radius:
                continue

            # Map amino acid to functional class
            functional_class = AA_TO_FUNCTIONAL[residue_type]
            if functional_class not in func_to_idx:
                continue

            # Calculate distance bin
            distance_bin = int(distance / bin_width)
            if distance_bin >= n_distance_bins:
                distance_bin = n_distance_bins - 1

            # Calculate feature index
            func_idx = func_to_idx[functional_class]
            ss_idx = ss_to_idx[secondary_structure]
            feature_idx = func_idx * len(self.ss_list) * n_distance_bins + ss_idx * n_distance_bins + distance_bin

            # Increment feature
            features[feature_idx] += 1.0

        return features

    def get_encoding_name(self) -> str:
        """Get name of this encoding strategy."""
        return "FUNCTIONAL_5"


class AA20ConservationEncoding(EncodingStrategy):
    """20-amino-acid encoding with conservation metrics.

    Combines 20-AA identity with Shannon entropy conservation scores.

    Encodes as: 20 AA × 3 SS × distance bins + entropy per distance bin
    Total features: 900 + 15 = 915 features (for 15 distance bins)

    This encoding tests whether evolutionary conservation at neighbor positions
    improves cluster coherence. Conservation scores come from sequence alignments
    (e.g., HHBLITS via PyCoM).

    Example:
        >>> encoding = AA20ConservationEncoding()
        >>> neighbors = [("L", "H", 3.5), ("A", "H", 5.2), ("V", "E", 8.1)]
        >>> conservation = [2.1, 1.8, 3.5]  # Shannon entropy scores
        >>> features = encoding.encode_neighbors(neighbors, n_distance_bins=15, conservation_scores=conservation)
        >>> features.shape
        (915,)
    """

    def __init__(self):
        """Initialize conservation encoding strategy."""
        self.aa_list = list(AA_ORDER)
        self.ss_list = list(SS_ORDER)
        self.logger = logging.getLogger(__name__)

    def get_n_features(self, n_distance_bins: int) -> int:
        """Calculate feature dimension.

        Args:
            n_distance_bins: Number of distance shells

        Returns:
            20 * 3 * n_distance_bins + n_distance_bins features
        """
        # Base AA20 features + entropy per distance bin
        return len(self.aa_list) * len(self.ss_list) * n_distance_bins + n_distance_bins

    def encode_neighbors(
        self,
        neighbors: List[Tuple[str, str, float]],
        n_distance_bins: int,
        radius: float = 15.0,
        bin_width: float = 1.0,
        conservation_scores: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Encode neighbors with amino acid identity + conservation.

        Args:
            neighbors: List of (residue_type, ss, distance) tuples
            n_distance_bins: Total number of distance bins
            radius: Maximum distance cutoff (Angstroms)
            bin_width: Width of each distance bin (Angstroms)
            conservation_scores: Optional array of conservation scores per neighbor

        Returns:
            Feature vector of shape (n_features,) as float32

        Details:
            - First 900D: Same as AA20Encoding
            - Last 15D: Average conservation score per distance bin
            - If conservation_scores not provided, defaults to zeros
            - Conservation score interpretation: Shannon entropy (0=conserved, ~3.8=random)
        """
        n_features = self.get_n_features(n_distance_bins)
        features = np.zeros(n_features, dtype=np.float32)

        # Build lookup tables
        aa_to_idx = {aa: i for i, aa in enumerate(self.aa_list)}
        ss_to_idx = {ss: i for i, ss in enumerate(self.ss_list)}

        # Track conservation scores per distance bin for averaging
        conservation_per_bin = np.zeros(n_distance_bins, dtype=np.float32)
        counts_per_bin = np.zeros(n_distance_bins, dtype=np.float32)

        # Process each neighbor
        for neighbor_idx, (residue_type, secondary_structure, distance) in enumerate(neighbors):
            # Validate inputs
            if residue_type not in aa_to_idx:
                continue
            if secondary_structure not in ss_to_idx:
                continue
            if distance > radius:
                continue

            # Calculate distance bin
            distance_bin = int(distance / bin_width)
            if distance_bin >= n_distance_bins:
                distance_bin = n_distance_bins - 1

            # Update AA20 features (same as AA20Encoding)
            aa_idx = aa_to_idx[residue_type]
            ss_idx = ss_to_idx[secondary_structure]
            feature_idx = aa_idx * len(self.ss_list) * n_distance_bins + ss_idx * n_distance_bins + distance_bin
            features[feature_idx] += 1.0

            # Accumulate conservation score for this distance bin
            if conservation_scores is not None and neighbor_idx < len(conservation_scores):
                conservation_per_bin[distance_bin] += conservation_scores[neighbor_idx]
                counts_per_bin[distance_bin] += 1.0

        # Compute average conservation per distance bin and append to features
        base_features = len(self.aa_list) * len(self.ss_list) * n_distance_bins
        for bin_idx in range(n_distance_bins):
            if counts_per_bin[bin_idx] > 0:
                avg_conservation = conservation_per_bin[bin_idx] / counts_per_bin[bin_idx]
            else:
                avg_conservation = 0.0  # No neighbors in this bin
            features[base_features + bin_idx] = avg_conservation

        return features

    def get_encoding_name(self) -> str:
        """Get name of this encoding strategy."""
        return "AA_20_CONSERVATION"


# Factory function for creating encoding strategies
def get_encoding_strategy(encoding: Encoding) -> EncodingStrategy:
    """Create encoding strategy from enum.

    Args:
        encoding: Encoding enum value

    Returns:
        EncodingStrategy instance

    Raises:
        ValueError: If encoding type is not recognized

    Example:
        >>> strategy = get_encoding_strategy(Encoding.FUNCTIONAL_5)
        >>> features = strategy.encode_neighbors(neighbors, n_distance_bins=15)
    """
    strategies: Dict[Encoding, EncodingStrategy] = {
        Encoding.AA_20: AA20Encoding(),
        Encoding.FUNCTIONAL_5: Functional5Encoding(),
        Encoding.AA_20_WITH_CONSERVATION: AA20ConservationEncoding(),
    }

    if encoding not in strategies:
        raise ValueError(f"Unknown encoding: {encoding}. Available: {list(strategies.keys())}")

    return strategies[encoding]


def get_encoding_dimension_info() -> Dict[str, int]:
    """Get feature dimension for each encoding at default radius (15Å).

    Returns:
        Dictionary mapping encoding name to feature dimension

    Example:
        >>> dims = get_encoding_dimension_info()
        >>> dims["AA_20"]
        900
        >>> dims["FUNCTIONAL_5"]
        225
    """
    n_distance_bins = 15  # Default 15Å radius with 1Å bin width

    return {
        "AA_20": 20 * 3 * n_distance_bins,  # 900
        "FUNCTIONAL_5": 5 * 3 * n_distance_bins,  # 225
        "AA_20_CONSERVATION": 20 * 3 * n_distance_bins + n_distance_bins,  # 915
    }
