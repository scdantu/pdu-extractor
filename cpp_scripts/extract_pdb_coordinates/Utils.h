#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <unordered_set>

#include "Constants.h"
#include "PDBContext.h"

std::string concatenateString(const std::vector<std::string>& strings);
std::string concatenateString(const std::unordered_set<std::string>& strings);

PDBType processHeader(const std::string &line, PDBContext &con);
void processRemark(const std::string &line, PDBContext &con);
void processDBRef(const std::string &line, PDBContext &con);
void processDBRef1(const std::string &line, PDBContext &con);
void processSequence(const std::string &line, PDBContext &con);

#endif // UTILS_H
