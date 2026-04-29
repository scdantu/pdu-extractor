#ifndef PDBCONTEXT_H
#define PDBCONTEXT_H

#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

struct PDBContext {
    // main data
    std::string pdbId;
    float resolution = -1.0f;
    std::unordered_set<std::string> uniprotIds;
    
    // input and output data
    std::vector<std::string> output;

    // sequence data
    std::stringstream parsedSequence;
    std::unordered_map<char, std::stringstream> sequenceStreams;

    // residue tracking
    int prevCAResiduePosition = -1;
    int firstCAResidue = 0;

    // condition flags
    bool hasResiduesOutOfOrder = true;
    bool anyCAAtomsPresent = false;
    bool isNotProtein = false;
    bool hasExcludedAminoAcid = false;

    // error tracking
    std::vector<std::string> errorOutput;

    void resetPDBOutput() {
        if (!anyCAAtomsPresent && output.size()) {
            anyCAAtomsPresent = true;
        }

        output.clear();
        hasResiduesOutOfOrder = true;
        hasExcludedAminoAcid = false;
        prevCAResiduePosition = -1;
        firstCAResidue = 0;
        parsedSequence.str("");
    }
};

#endif // PDBCONTEXT_H