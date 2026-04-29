#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include <unordered_set>

#include "PDBContext.h"
#include "AtomDataParser.h"
#include "Constants.h"
#include "Utils.h"


ResidueConfirmation validateAtomSequence(int &prevCAResiduePosition, const int &resSeq, int &firstCAResidue, std::vector<std::string> &errorOutput) {
    if (prevCAResiduePosition + 1 != resSeq) {
        if (prevCAResiduePosition == -1) { // initial value
            prevCAResiduePosition = resSeq;
            firstCAResidue = resSeq;
            return RESIDUE_VALID;
        }

        if (prevCAResiduePosition == resSeq)
            return RESIDUE_DUPLICATE;

        std::stringstream errorString;
        errorString << "missing residues; prev=" << prevCAResiduePosition << ", next=" << resSeq;
        errorOutput.push_back(errorString.str());

        prevCAResiduePosition = resSeq;
        return RESIDUE_OUT_OF_SEQUENCE;
    }
    prevCAResiduePosition = resSeq;
    return RESIDUE_VALID;
}

void processAtom(const std::string &line, PDBContext &con) {
    AtomData data(line, 0);
    if (!data.isValidAtom)
        return;

    switch(validateAtomSequence(con.prevCAResiduePosition, data.resSeq, con.firstCAResidue, con.errorOutput)) {
    case RESIDUE_VALID:
        break; // continue
    case RESIDUE_DUPLICATE:
        return; // skip to next
    case RESIDUE_OUT_OF_SEQUENCE:
        con.hasResiduesOutOfOrder = false;
        break;
    }

    char aminoAcid;
    try {
        aminoAcid = aminoAcidLookup.at(data.resName);
    } catch (std::out_of_range) { // should never throw if pdb is valid
        throw std::runtime_error("Unexpected atom type: " + data.resName);
    }

    // Selenocysteine, Pyrrolysine, GLX, ASX, or unknown
    if (invalidAminoAcids.find(aminoAcid) != invalidAminoAcids.end())
        con.hasExcludedAminoAcid = true;

    // construct output string
    std::stringstream ss;
    ss << aminoAcid << ' ' << data.resName << ' ' << data.chainId << ' ' << data.resSeq << ' '
       << data.insertionCode << ' ' << data.x << ' '  << data.y << ' ' << data.z;
    // std::cout << aminoAcid << ' ' << data.x << ' '  << data.y << ' ' << data.z << std::endl;

    // construct sequence string
    con.parsedSequence << aminoAcid;

    con.output.push_back(ss.str());
}

// Checks whether the parsed input, so far, produced a valid, sequential
// list of residues with coordinates.
// Returns PDBParsingCode.SUCCESS if successful, and a specific error code
// otherwise.
PDBParsingCode isPDBInvalid(PDBContext &con) {
    if (con.isNotProtein)
        return IS_NOT_PROTEIN;
    if (con.hasExcludedAminoAcid)
        return EXCLUDE_UNKNOWN_OR_RARE_AMINO_ACIDS;

    bool isResolutionValid = con.resolution < MAX_RESOLUTION;
    if (!isResolutionValid) // resolution too low
        return RESOLUTION_TOO_LOW;

    if (con.resolution == -1) // no valid resolution remark returned
        return RESOLUTION_NOT_SPECIFIED;

    if (!con.hasResiduesOutOfOrder) // missing non-terminal residues
        return MISSING_NON_TERMINAL_RESIDUES;

    if (con.prevCAResiduePosition == -1) { // no single CA atom found
        if (con.anyCAAtomsPresent) // if any model had, but last one didn't
            return MISSING_NON_TERMINAL_RESIDUES;
        return NO_ALPHA_CARBON_ATOMS_FOUND;
    }

    if (con.uniprotIds.size() == 0)
        return NO_UNIPROT_ID;

    return SUCCESS;
}

std::vector<std::string> processSequences(std::unordered_map<char, std::stringstream> &sequenceStreams, std::stringstream &parsedSequence) {
    std::unordered_set<std::string> uniqueSequences;
    std::vector<std::string> outputSeq;
    std::string matchedSequence = "N/A";

    if (sequenceStreams.size() == 0)
        return outputSeq;

    for (const auto & [_chainId, stream] : sequenceStreams) {
        auto sequence = stream.str();
        if (sequence.size() != 0) {
            auto sequencePosition = sequence.find(parsedSequence.str());
            if (parsedSequence.str().size() > 0 && sequencePosition != std::string::npos) {
                matchedSequence = sequence;
            } else {
                uniqueSequences.insert(sequence);
            }
        }
    }

    // line 5: matched sequence (parsed contained within matched)
    if (matchedSequence != "") {
        outputSeq.push_back("matched: " + matchedSequence);
    }

    // line 6: sequence parsed from ATOM records
    outputSeq.push_back("parsed:  " + parsedSequence.str());

    // sequences.push_back(matchedSequence);

    // line 7+: all other parsed sequences
    for (std::string seq: uniqueSequences) {
        outputSeq.push_back("other:   " + seq);
    }

    return outputSeq;
}

void printOutput(PDBContext &con, bool valid) {
    // line 1 -- validity (0=invalid, 1=valid)
    std::cout << "success: " << valid << std:: endl;

    // line 2 -- pdb id
    // "pdb_id:  201L" (printed in Utils.cpp)
    std::cout << "pdb_id:  " << con.pdbId << std::endl;

    // line 3 -- resolution
    std::cout << "resolut: " << con.resolution << std::endl;

    // line 4 -- uniprot IDs
    std::string allUniprotIds = concatenateString(con.uniprotIds);
    std::cout << "uniprot: " << allUniprotIds << std::endl;

    // line 5 -- matched sequence (atom record substring of reqres)
    // line 6 -- parsed sequence (atom records)
    // line 7-n -- other sequences (reqres sequence)
    for (auto lineSeq: processSequences(con.sequenceStreams, con.parsedSequence))
        std::cout << lineSeq << std::endl;

    if (!valid) // stop printing if invalid
        return;

    // line n+1: sequence number of initial residue (starts with 1)
    std::cout << "initres: " << con.firstCAResidue << std::endl;

    // line n+2: empty line
    std::cout << std::endl;

    // lines n+3 to end: coordinates in format
    // <one_letter_residue> <three_letter_residue> <chain_id> <residue_number> <insertion_code> <x> <y> <z>
    for (std::string pos : con.output) {
        std::cout << pos << std::endl;
    }
}

// Takes in a stream of a PDB file as input.
PDBParsingCode processPDBStream(std::istream& in) {
    PDBContext con;
    
    std::string line;
    while (getline(in, line)) {
        std::string param = line.substr(0, 6);

        if (param == "HEADER") {
            auto headerType = processHeader(line, con);
            if (headerType != PROTEIN) {
                con.isNotProtein = true;
                break;
            }
        } else if (param == "REMARK") {
            processRemark(line, con);
            if (con.resolution > -1 && con.resolution > MAX_RESOLUTION) // resolution bad
                break;
        } else if (param == "DBREF ") {
            processDBRef(line, con);
        } else if (param == "DBREF1") {
            processDBRef1(line, con);
        } else if (param == "SEQRES") {
            processSequence(line, con);
        } else if (param == "ATOM  ") { // HETATM residues are skipped
            processAtom(line, con);
        } else if (param == "TER   ") { // end of one chain
            auto pdbValidity = isPDBInvalid(con);
            // std::cout << code_name[pdbValidity] << std::endl;
            if (pdbValidity == SUCCESS)
                break; // terminate parser, output PDB

            // if at first you don't succeed, try, try again (parse next model)
            con.resetPDBOutput();
        } // else ignore line, until end is reached
    }
    
    auto pdbValidity = isPDBInvalid(con);
    if (pdbValidity != SUCCESS) {
        printOutput(con, false);
        std::cerr << code_name[pdbValidity] << std::endl;

        for (auto error : con.errorOutput)
            std::cout << error << std::endl;
        return pdbValidity;
    }

    printOutput(con, true);

    return SUCCESS;
}

// Part of a previouslly planned functionality
// std::stringstream decompressGZFile(const std::string& filename) {
//     gzFile gzfile = gzopen(filename.c_str(), "rb");
//     if (!gzfile) {
//         std::cout << "success: 0" << std::endl;
//         std::cerr << "Failed to open file " << filename << std::endl;
//         return std::stringstream();
//     }
    
//     char buffer[4096];
//     std::stringstream ss;
//     int num_read = 0;
//     while ((num_read = gzread(gzfile, buffer, sizeof(buffer))) > 0) {
//         ss.write(buffer, num_read);
//     }

//     gzclose(gzfile);
//     return ss;
// }

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    return processPDBStream(std::cin);
    
    // if(argc != 2 || std::string(argv[1]) == "-h") {
    //     std::cerr << "Usage: " << argv[0] << " (--process-file-stream | --process-files)" << std::endl;
    //     return 1;
    // }
    
    // std::string arg(argv[1]);
    // if(arg == "--process-file-stream") {
    //     return processPDBStream(std::cin);
    // } else if(arg == "--process-files") {
    //     std::string line;
    //     while(std::getline(std::cin, line) && line != "/") {
    //         std::cout << line << std::endl;
    //         std::stringstream ss = decompressGZFile(line);
    //         if(ss.rdbuf()->in_avail() > 0) {
    //             processPDBStream(ss);
    //         }
    //     }
    //     return 0;
    // } else {
    //     std::cerr << "Invalid argument. Usage: " << argv[0] << " (--process-file-stream | --process-files)" << std::endl;
    //     return 1;
    // }
}
