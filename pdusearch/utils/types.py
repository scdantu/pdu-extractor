"""Type definitions and constants for PDU extraction and analysis."""

from enum import Enum
from typing import List, Tuple

# Amino acid constants
AA_ORDER: List[str] = [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
]

AA_TO_INDEX: dict = {aa: i for i, aa in enumerate(AA_ORDER)}

# Secondary structure constants
SS_ORDER: Tuple[str, ...] = ("H", "E", "C")
SS_TO_INDEX: dict = {ss: i for i, ss in enumerate(SS_ORDER)}

# Functional residue classes (physicochemical grouping)
class FunctionalClass(Enum):
    """Functional residue classes for feature encoding."""
    HYDROPHOBIC = 0  # A, V, L, I, M, F, W, P
    POLAR = 1        # S, T, N, Q, C
    POSITIVE = 2     # K, R, H
    NEGATIVE = 3     # D, E
    SPECIAL = 4      # G, Y


# Mapping from amino acid to functional class
AA_TO_FUNCTIONAL: dict = {
    "A": FunctionalClass.HYDROPHOBIC,
    "V": FunctionalClass.HYDROPHOBIC,
    "L": FunctionalClass.HYDROPHOBIC,
    "I": FunctionalClass.HYDROPHOBIC,
    "M": FunctionalClass.HYDROPHOBIC,
    "F": FunctionalClass.HYDROPHOBIC,
    "W": FunctionalClass.HYDROPHOBIC,
    "P": FunctionalClass.HYDROPHOBIC,
    "S": FunctionalClass.POLAR,
    "T": FunctionalClass.POLAR,
    "N": FunctionalClass.POLAR,
    "Q": FunctionalClass.POLAR,
    "C": FunctionalClass.POLAR,
    "K": FunctionalClass.POSITIVE,
    "R": FunctionalClass.POSITIVE,
    "H": FunctionalClass.POSITIVE,
    "D": FunctionalClass.NEGATIVE,
    "E": FunctionalClass.NEGATIVE,
    "G": FunctionalClass.SPECIAL,
    "Y": FunctionalClass.SPECIAL,
}

FUNCTIONAL_CLASS_ORDER: List[str] = ["HYDROPHOBIC", "POLAR", "POSITIVE", "NEGATIVE", "SPECIAL"]
FUNCTIONAL_CLASS_ORDER_SHORT: List[str] = ["HY", "PO", "PO+", "PO-", "SP"]

# Constants for analysis
DISTANCE_BIN_WIDTH: float = 1.0
DEFAULT_RADIUS: float = 15.0
DEFAULT_MIN_CLUSTER_SIZE: int = 200
DEFAULT_MIN_SAMPLES: int = 5

# Residue class order (from kmers module)
RESIDUE_CLASS_ORDER: List[str] = [
    "Hydrophobic",
    "Polar_uncharged",
    "Polar_charged_positive",
    "Polar_charged_negative",
    "Special"
]
