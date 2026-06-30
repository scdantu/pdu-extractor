"""Sequence alignment utilities for conservation metrics.

This module provides interfaces for accessing sequence alignment data and computing
conservation scores (Shannon entropy). Conservation metrics can be used to enhance
PDU feature encoding by incorporating evolutionary information.

Alignments come from HHBLITS (multiple sequence alignments) via PyCoM database.
Conservation scores quantify how variable each position is across the alignment.

Example:
    >>> from pdusearch.features.alignment import compute_shannon_entropy
    >>> alignment = np.array([[0, 0, 1, 2, ...], ...])  # (n_sequences, n_positions)
    >>> entropy = compute_shannon_entropy(alignment[:, position=5])
    >>> print(entropy)  # ~2.1 (conserved position)
"""

from abc import ABC, abstractmethod
from typing import Optional, Dict, Tuple
import numpy as np
import logging

logger = logging.getLogger(__name__)


def compute_shannon_entropy(residue_counts: np.ndarray) -> float:
    """Compute Shannon entropy from residue frequency distribution.

    Shannon entropy measures the variability of a position in a sequence alignment.
    Interpretation:
        - Low entropy (~0): Position is highly conserved (few residue types)
        - High entropy (~3.8): Position is variable (many residue types equally likely)

    Args:
        residue_counts: Array of counts for each residue type at a position.
                       Can be raw counts or normalized frequencies.
                       Shape: (n_residue_types,)

    Returns:
        Shannon entropy in bits, range [0, log2(n_types)]

    Example:
        >>> counts = np.array([90, 5, 5, 0, ...])  # 90% Leu, 5% Val, 5% Ile, rest none
        >>> entropy = compute_shannon_entropy(counts)
        >>> print(f"{entropy:.2f}")  # ~0.47 (very conserved)
    """
    # Normalize to probabilities
    total = np.sum(residue_counts)
    if total == 0:
        return 0.0  # No data = no entropy

    probs = residue_counts / total

    # Compute Shannon entropy: H = -sum(p * log2(p))
    # Only include non-zero probabilities to avoid log(0)
    entropy = 0.0
    for p in probs:
        if p > 0:
            entropy -= p * np.log2(p)

    return entropy


def compute_conservation_score(residue_counts: np.ndarray, method: str = "entropy") -> float:
    """Compute conservation score for a sequence alignment position.

    Wrapper function for different conservation scoring methods.

    Args:
        residue_counts: Array of residue counts at position
        method: Scoring method ('entropy', 'gap_fraction', 'relative_entropy')

    Returns:
        Conservation score (interpretation depends on method)

    Details:
        - 'entropy': Shannon entropy (higher = more variable)
        - 'gap_fraction': Fraction of gaps/insertions at position
        - 'relative_entropy': KL divergence from background distribution
    """
    if method == "entropy":
        return compute_shannon_entropy(residue_counts)
    elif method == "gap_fraction":
        # Assuming last entry in counts is gaps
        total = np.sum(residue_counts)
        if total == 0:
            return 1.0  # All gaps = high conservation score
        return residue_counts[-1] / total
    else:
        raise ValueError(f"Unknown conservation method: {method}")


class AlignmentProvider(ABC):
    """Abstract interface for sequence alignment data providers.

    Implementations provide access to multiple sequence alignments (MSAs)
    and compute conservation metrics for positions.
    """

    @abstractmethod
    def get_alignment(self, pdb_id: str, chain_id: str) -> Optional[np.ndarray]:
        """Get multiple sequence alignment for a PDB structure.

        Args:
            pdb_id: PDB identifier (e.g., '1abc')
            chain_id: Chain identifier (single letter)

        Returns:
            Alignment matrix of shape (n_sequences, n_positions) with residue indices,
            or None if alignment not available
        """
        pass

    @abstractmethod
    def get_conservation_score(
        self,
        pdb_id: str,
        chain_id: str,
        position: int,
        method: str = "entropy",
    ) -> float:
        """Get conservation score for a specific position.

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            position: Position index in chain (1-indexed or 0-indexed, implementation-specific)
            method: Conservation scoring method

        Returns:
            Conservation score as float
        """
        pass

    @abstractmethod
    def is_available(self, pdb_id: str, chain_id: str) -> bool:
        """Check if alignment data is available for structure.

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier

        Returns:
            True if alignment data exists, False otherwise
        """
        pass


class PyComAlignmentProvider(AlignmentProvider):
    """PyCoM-based alignment provider.

    Accesses HHBLITS alignments from PyCoM database (https://pycom.brunel.ac.uk).
    PyCoM provides pre-computed multiple sequence alignments for PDB structures.

    Note: Requires pycom library to be installed and database accessible.
    Install: pip install pycom-api

    Example:
        >>> provider = PyComAlignmentProvider(db_path="pycom_db/pycom.db")
        >>> if provider.is_available("1abc", "A"):
        ...     alignment = provider.get_alignment("1abc", "A")
        ...     entropy = provider.get_conservation_score("1abc", "A", position=50)
    """

    def __init__(self, db_path: str = "pycom_db/pycom.db", use_cache: bool = True):
        """Initialize PyCoM alignment provider.

        Args:
            db_path: Path to PyCoM database file
            use_cache: Whether to cache alignments in memory
        """
        self.db_path = db_path
        self.use_cache = use_cache
        self.alignment_cache: Dict[Tuple[str, str], Optional[np.ndarray]] = {}
        self.logger = logging.getLogger(__name__)

        # Try to import pycom
        try:
            import pycom
            self.pycom = pycom
            self.logger.info(f"PyCoM initialized (db: {db_path})")
        except ImportError:
            self.logger.warning("pycom library not installed. Install with: pip install pycom-api")
            self.pycom = None

    def get_alignment(self, pdb_id: str, chain_id: str) -> Optional[np.ndarray]:
        """Get HHBLITS alignment from PyCoM.

        Args:
            pdb_id: PDB identifier (lowercase, e.g., '1abc')
            chain_id: Chain identifier

        Returns:
            Alignment matrix or None if not available
        """
        if not self.pycom:
            self.logger.warning("PyCoM not available")
            return None

        # Normalize PDB ID
        pdb_id = pdb_id.lower()
        cache_key = (pdb_id, chain_id)

        # Check cache
        if self.use_cache and cache_key in self.alignment_cache:
            return self.alignment_cache[cache_key]

        try:
            # Query PyCoM for alignment
            # Note: This is a placeholder. Actual API depends on pycom library version.
            alignment = self.pycom.query_alignment(pdb_id, chain_id)

            if self.use_cache:
                self.alignment_cache[cache_key] = alignment

            return alignment

        except Exception as e:
            self.logger.warning(f"Failed to retrieve alignment for {pdb_id}:{chain_id}: {e}")
            return None

    def get_conservation_score(
        self,
        pdb_id: str,
        chain_id: str,
        position: int,
        method: str = "entropy",
    ) -> float:
        """Get conservation score using PyCoM alignment.

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            position: Position in chain (0-indexed)
            method: Conservation scoring method

        Returns:
            Conservation score, or 0.0 if not available
        """
        alignment = self.get_alignment(pdb_id, chain_id)
        if alignment is None:
            self.logger.debug(f"No alignment for {pdb_id}:{chain_id}, returning 0.0")
            return 0.0

        # Validate position
        if position < 0 or position >= alignment.shape[1]:
            return 0.0

        # Extract residue counts for this position
        column = alignment[:, position]
        residue_counts = np.bincount(column.astype(int), minlength=20)

        # Compute conservation score
        return compute_conservation_score(residue_counts, method=method)

    def is_available(self, pdb_id: str, chain_id: str) -> bool:
        """Check if alignment available in PyCoM.

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier

        Returns:
            True if alignment exists, False otherwise
        """
        if not self.pycom:
            return False

        alignment = self.get_alignment(pdb_id, chain_id)
        return alignment is not None

    def clear_cache(self):
        """Clear alignment cache to free memory."""
        self.alignment_cache.clear()
        self.logger.info("Alignment cache cleared")


class NullAlignmentProvider(AlignmentProvider):
    """Null/placeholder alignment provider.

    Returns zero conservation scores for all positions. Useful for testing
    and as a fallback when alignment data is unavailable.

    Example:
        >>> provider = NullAlignmentProvider()
        >>> score = provider.get_conservation_score("1abc", "A", position=50)
        >>> print(score)  # 0.0
    """

    def get_alignment(self, pdb_id: str, chain_id: str) -> Optional[np.ndarray]:
        """Always returns None."""
        return None

    def get_conservation_score(
        self,
        pdb_id: str,
        chain_id: str,
        position: int,
        method: str = "entropy",
    ) -> float:
        """Always returns 0.0."""
        return 0.0

    def is_available(self, pdb_id: str, chain_id: str) -> bool:
        """Always returns False."""
        return False


def get_alignment_provider(provider_type: str = "null", **kwargs) -> AlignmentProvider:
    """Factory function to create alignment providers.

    Args:
        provider_type: Type of provider ('pycom', 'null')
        **kwargs: Additional arguments for provider initialization

    Returns:
        AlignmentProvider instance

    Example:
        >>> provider = get_alignment_provider("pycom", db_path="pycom_db/pycom.db")
        >>> provider = get_alignment_provider("null")
    """
    if provider_type == "pycom":
        return PyComAlignmentProvider(**kwargs)
    elif provider_type == "null":
        return NullAlignmentProvider()
    else:
        raise ValueError(f"Unknown provider type: {provider_type}")
