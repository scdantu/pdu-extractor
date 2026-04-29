#ifndef ATOMDATAPARSER_H
#define ATOMDATAPARSER_H

#include <string>
#include <sstream>

struct AtomData
{
    bool isValidAtom; // whether the atom is valid (CA)
    std::string resName; // residue name (AA)
    std::string chainId;
    int resSeq; // residue sequence number
    std::string insertionCode;
    float x, y, z;

    AtomData(const std::string& str, size_t offset = 0)
    {
        std::string atom_name = str.substr(-offset + 12, 4);
    
        isValidAtom = atom_name == " CA "; // whether the atom is ca
        if (!isValidAtom)
            return;

        resName = str.substr(-offset + 17, 3);
        chainId = str.substr(-offset + 21, 1);
        if (chainId == " ")
            chainId = ".";
        resSeq = std::stoi(str.substr(-offset + 22, 4));
        insertionCode = str.substr(-offset + 26, 1);
        if (insertionCode == " ")
            insertionCode = ".";
        x = std::stof(str.substr(-offset + 30, 8));
        y = std::stof(str.substr(-offset + 38, 8));
        z = std::stof(str.substr(-offset + 46, 8));
    }
};

#endif // ATOMDATAPARSER_H
