#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <unordered_map>
#include <unordered_set>
#include <string>

#define PDB_PARSING_CODES \
X(SUCCESS, "SUCCESS") \
X(RESOLUTION_TOO_LOW, "RESOLUTION_TOO_LOW") \
X(RESOLUTION_NOT_SPECIFIED, "RESOLUTION_NOT_SPECIFIED") \
X(MISSING_NON_TERMINAL_RESIDUES, "MISSING_NON_TERMINAL_RESIDUES") \
X(NO_ALPHA_CARBON_ATOMS_FOUND, "NO_ALPHA_CARBON_ATOMS_FOUND") \
X(IS_NOT_PROTEIN, "IS_NOT_PROTEIN") \
X(EXCLUDE_UNKNOWN_OR_RARE_AMINO_ACIDS, "EXCLUDE_UNKNOWN_OR_RARE_AMINO_ACIDS") \
X(HAS_UNKNOWN_RESIDUE, "HAS_UNKNOWN_RESIDUE") \
X(INVALID_SEQUENCE, "INVALID_SEQUENCE") \
X(NO_UNIPROT_ID, "NO_UNIPROT_ID") \

#define X(code, name) code,
enum PDBParsingCode : size_t {
    PDB_PARSING_CODES
    MAX_PDB_PARSING_CODES
};
#undef X

extern const char *code_name[MAX_PDB_PARSING_CODES];

enum ResidueConfirmation {
    RESIDUE_VALID,
    RESIDUE_DUPLICATE, // same residue multiple times (e.g. multiple confirmation)
    RESIDUE_OUT_OF_SEQUENCE // missing non-terminal residue
};

enum PDBType {PROTEIN, DNA, RNA, MISC};

extern const float MAX_RESOLUTION;
extern const std::unordered_map<std::string, char> aminoAcidLookup;
extern const std::unordered_set<char> invalidAminoAcids;

#endif // CONSTANTS_H

