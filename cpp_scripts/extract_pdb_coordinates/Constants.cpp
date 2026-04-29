#include "Constants.h"

#include <unordered_set>

#define X(code, name) name,
const char *code_name[] = {
    PDB_PARSING_CODES
};
#undef X

const float MAX_RESOLUTION = 2.5f;
const std::unordered_map<std::string, char> aminoAcidLookup = {
    {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'},
    {"CYS", 'C'}, {"GLN", 'Q'}, {"GLU", 'E'}, {"GLY", 'G'},
    {"HIS", 'H'}, {"HIP", 'H'}, {"HIE", 'H'}, {"ILE", 'I'},
    {"LEU", 'L'}, {"LYS", 'K'}, {"MET", 'M'}, {"PHE", 'F'},
    {"PRO", 'P'}, {"SER", 'S'}, {"THR", 'T'}, {"TYR", 'Y'},
    {"TRP", 'W'}, {"VAL", 'V'}, {"SEC", 'U'}, {"PYL", 'O'},
    {"XPL", 'O'}, // for pdb 1L2Q
    {"GLX", 'Z'}, // for pdb 1KP0 
    {"ASX", 'B'}, // for pdb 1KP0
    {"UNK", '.'} // unknown AA
};

// rare amino acids = SELENOCYSTEINE, PYRROLYSINE, ..., & unknown AA
const std::unordered_set<char> invalidAminoAcids{'U', 'O', 'Z', 'B', '.'};
