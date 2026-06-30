"""PDU feature extraction module.

This module provides configurable feature extraction from PDU database with support
for multiple encoding strategies and optional conservation metrics.

Main Classes:
    - FeatureExtractor: Main class for extracting features
    - Encoding: Enum of available encodings (AA_20, FUNCTIONAL_5, AA_20_WITH_CONSERVATION)

Encoding Strategies:
    - AA20Encoding: Standard 20-amino-acid encoding (900D for 15Å)
    - Functional5Encoding: Functional class grouping (225D for 15Å)
    - AA20ConservationEncoding: 20-AA + conservation metrics (915D for 15Å)

Alignment Providers:
    - PyComAlignmentProvider: HHBLITS alignments from PyCoM database
    - NullAlignmentProvider: Placeholder (returns zeros)

Example:
    >>> from pdusearch.features import FeatureExtractor, Encoding
    >>> from pdusearch.config import Config
    >>>
    >>> config = Config(db="per_aa_sqlite/pdus_L.sqlite")
    >>> extractor = FeatureExtractor(config, encoding=Encoding.FUNCTIONAL_5)
    >>> X, pdu_ids = extractor.extract_for_aa("L")
    >>> print(f"Extracted {X.shape[0]:,} PDUs × {X.shape[1]} features")
"""

from .encodings import (
    Encoding,
    EncodingStrategy,
    AA20Encoding,
    Functional5Encoding,
    AA20ConservationEncoding,
    get_encoding_strategy,
    get_encoding_dimension_info,
)

from .alignment import (
    AlignmentProvider,
    PyComAlignmentProvider,
    NullAlignmentProvider,
    compute_shannon_entropy,
    compute_conservation_score,
    get_alignment_provider,
)

from .extractor import FeatureExtractor

__all__ = [
    # Main class
    "FeatureExtractor",
    # Encodings
    "Encoding",
    "EncodingStrategy",
    "AA20Encoding",
    "Functional5Encoding",
    "AA20ConservationEncoding",
    "get_encoding_strategy",
    "get_encoding_dimension_info",
    # Alignments
    "AlignmentProvider",
    "PyComAlignmentProvider",
    "NullAlignmentProvider",
    "compute_shannon_entropy",
    "compute_conservation_score",
    "get_alignment_provider",
]
