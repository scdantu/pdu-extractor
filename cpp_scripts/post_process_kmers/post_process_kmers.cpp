#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <map>
#include <iomanip>

struct PdbInfo {
    std::string pdb_id;
    double resolution;
    int sequence_length;
};

namespace fs = std::filesystem;

std::unordered_map<std::string, int> global_kmers;
bool process_all_pdbs = false;
int kmer_size = 12;

std::vector<std::string> readUniprotFiles(const fs::path& uniprot_path) {
    std::vector<std::string> file_list;
    for(const auto& entry : fs::directory_iterator(uniprot_path)) {
        file_list.push_back(entry.path().string());
    }
    return file_list;
}

PdbInfo parseLine(const std::string& line) {
    std::istringstream iss(line);
    PdbInfo info;
    iss >> info.pdb_id >> info.resolution >> info.sequence_length;
    return info;
}

PdbInfo selectPdb(const std::vector<PdbInfo>& pdb_infos) {
    std::vector<PdbInfo> sorted_infos = pdb_infos;
    std::sort(sorted_infos.begin(), sorted_infos.end(), [](const PdbInfo& a, const PdbInfo& b) {
        return a.resolution < b.resolution || (a.resolution == b.resolution && a.sequence_length > b.sequence_length);
    });
    
    for (size_t i = 0; i < sorted_infos.size() - 1; ++i) {
        if (static_cast<double>(sorted_infos[i].sequence_length) >= 0.8 * sorted_infos[i + 1].sequence_length) {
            return sorted_infos[i];
        }
    }
    return sorted_infos.back();
}

PdbInfo parseInfoFile(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line;
    std::vector<PdbInfo> pdb_infos;
    while(std::getline(file, line)) {
        pdb_infos.push_back(parseLine(line));
    }

    return selectPdb(pdb_infos);
}

void parseKmersFile(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line;
    while(std::getline(file, line)) {
        if(line.length() >= kmer_size) {
            global_kmers[line.substr(0, kmer_size)]++;
        }
    }
}

int main(int argc, char** argv) {
    for(int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if(arg == "-a") {
            process_all_pdbs = true;
        } else if(arg == "-k" && i + 1 < argc) {
            kmer_size = std::stoi(argv[++i]);
        } else if(arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [OPTIONS]\n"
                      << "Options:\n"
                      << "  -a            Process all PDBs\n"
                      << "  -k <value>    Specify the size of the k-mers\n"
                      << "  -h, --help    Display this help message and exit\n";
            return 0;
        }
    }

    fs::path uniprot_path = "./pdb_output/uniprot";
    fs::path pdbs_path = "./pdb_output/pdbs";
    std::vector<std::string> file_list;

    if(process_all_pdbs) {
        for(const auto& entry : fs::directory_iterator(pdbs_path)) {
            file_list.push_back(entry.path().c_str());
        }
    } else {
        file_list = readUniprotFiles(uniprot_path);
    }

    int total_files = file_list.size();

    for(int i = 0; i < total_files; ++i) {
        if(process_all_pdbs) {
            parseKmersFile(file_list[i]);
        } else {
            PdbInfo selected_pdb = parseInfoFile(file_list[i]);

            std::string kmers_file_path = (pdbs_path / (selected_pdb.pdb_id + ".kmers")).c_str();
            if(fs::exists(kmers_file_path)) {
                parseKmersFile(kmers_file_path);
            }
        }

        if((i + 1) % 100 == 0 || i + 1 == total_files) {
            std::cerr << "\rProcessed " << (i + 1) << " / " << total_files << "; " 
                      << std::fixed << std::setprecision(2) << static_cast<double>(i + 1) / total_files * 100 << "%";
            std::cerr.flush();
        }
    }

    std::cerr << std::endl << "Prepairing results..." << std::endl;

    std::vector<std::pair<std::string, int>> sorted_kmers(global_kmers.begin(), global_kmers.end());

    std::sort(sorted_kmers.begin(), sorted_kmers.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    for(const auto& [kmer, freq] : sorted_kmers) {
        std::cout << kmer << " " << freq << std::endl;
    }

    return 0;
}
