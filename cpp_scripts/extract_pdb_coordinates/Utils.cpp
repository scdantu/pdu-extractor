#include "Utils.h"
#include <sstream>
#include <iterator>
#include <unordered_set>
#include <iostream>

#include "Constants.h"
#include "PDBContext.h"

std::string concatenateString(const std::vector<std::string>& strings) {
    const char delim = ',';
    std::ostringstream oss;

    if (!strings.empty()) {
        // Convert all but the last element to avoid a trailing delimiter
        std::copy(strings.begin(), strings.end()-1,
            std::ostream_iterator<std::string>(oss, &delim));

        // Now add the last element with no delimiter
        oss << strings.back();
    }

    return oss.str();
}

std::string concatenateString(const std::unordered_set<std::string>& strings) {
    const char delim = ',';
    std::ostringstream oss;

    for (auto itr = strings.begin(); itr != strings.end(); ++itr) {
        if (itr != strings.begin()) {
            oss << delim;
        }
        oss << *itr;
    }

    return oss.str();
}

// Extracts the resolution from remark 2
float extractResolution(const std::string &line) {
    std::string res_section = line.substr(23, 7);
    // skip empty remark line
    if (res_section == "       ") {
        return -1;
    }
    try {
        return std::stof(line.substr(23, 7));
    } catch(std::invalid_argument) {
        return -2; // RESOLUTION. NOT APPLICABLE. (e.g. pdb 134D)
    }

    return std::stof(line.substr(23, 7));
}

/////////////////////
/// PROCESS ENTRY ///
/////////////////////

PDBType processHeader(const std::string &line, PDBContext &con) {
    std::string cls = line.substr(10, 40); // 11-50
    std::string pdbId = line.substr(62, 4); // 63-66

    con.pdbId = pdbId;
 
    if (cls.find("DNA") != std::string::npos) {
        if (cls.find("DNA BINDING PROTEIN") == std::string::npos)
            return DNA;
    }
    if (cls.find("RNA") != std::string::npos) {
        if (cls.find("RNA BINDING PROTEIN") == std::string::npos)
            return RNA;
    }

    return PROTEIN;
}

// Remark row
void processRemark(const std::string &line, PDBContext &con) {
    int remark_no = std::stoi(line.substr(7, 3));
    switch (remark_no) {
        case 2:
            int extractedRes = extractResolution(line);
            if (extractedRes != -1) {
                con.resolution = extractResolution(line);
            }
            break;
    }
}

// DBRef row
void processDBRef(const std::string &line, PDBContext &con) {
    std::string db = line.substr(26, 6); // 27 - 32
    if (db != "UNP   ") // only match uniprot
        return;

    std::string uniprotId = line.substr(33, 8); // 34 - 41
    std::stringstream parser(uniprotId);
    parser >> uniprotId;
    con.uniprotIds.insert(uniprotId);
}

void processDBRef1(const std::string &line, PDBContext &con) {
    // process 1 for uniprot
    std::string db = line.substr(26, 6); // 27 - 32
    if (db != "UNP   ") // only match uniprot
        return;

    // process 2 for id
    std::string nextLine;
    getline(std::cin, nextLine);

    std::string uniprotId = nextLine.substr(18, 22); // 19 - 40
    std::stringstream parser(uniprotId);
    parser >> uniprotId;
    con.uniprotIds.insert(uniprotId);
}

// SEQRES row
void processSequence(const std::string &line, PDBContext &con) {
    char chainId = line[11];
    std::string aa;
    std::string aaLine = line.substr(19, 51);
    std::stringstream aaStream(aaLine);
    while (aaStream >> aa) {
        try {
            con.sequenceStreams[chainId] << aminoAcidLookup.at(aa);
        } catch (std::out_of_range) {
            // replace non-standard AA with dot (.)
            con.sequenceStreams[chainId] << '.';
            return;
        }
    }
}