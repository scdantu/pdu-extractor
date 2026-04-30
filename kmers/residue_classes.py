AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"

RESIDUE_CLASS_ORDER = (
    "hydrophobic",
    "aromatic",
    "polar",
    "positive",
    "negative",
    "special",
)

RESIDUE_CLASS = {
    "A": "hydrophobic",
    "V": "hydrophobic",
    "L": "hydrophobic",
    "I": "hydrophobic",
    "M": "hydrophobic",
    "F": "aromatic",
    "Y": "aromatic",
    "W": "aromatic",
    "S": "polar",
    "T": "polar",
    "N": "polar",
    "Q": "polar",
    "C": "polar",
    "K": "positive",
    "R": "positive",
    "H": "positive",
    "D": "negative",
    "E": "negative",
    "G": "special",
    "P": "special",
}


def residue_class(residue):
    return RESIDUE_CLASS.get(residue, "unknown")
