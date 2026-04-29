#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>

// using namespace std;

void printHelp() {
    std::cout << "Usage: ./program [-k <kmer size>] [-f] [-s] [--help]\n"
         << "Options:\n"
         << "  -k <kmer size>     Set the size of the kmers (default=12)\n"
         << "  -f                 Output frequencies only (default=false)\n"
         << "  -s                 Output statistics and missing kmers (default=false)\n"
         << "  --help, -h         Show this help message\n";
}

// Generate synthetic set using recursion
// void generateSyntheticSet(const string& kmer, int k, const unordered_map<string, int>& naturalSet, const function<void(const string&)>& callback) {
//     static string aminoAcids = "ACDEFGHIKLMNPQRSTVWY";

//     if (k == 0) {
//         if (naturalSet.find(kmer) == naturalSet.end()) {
//             callback(kmer);
//         }
//         return;
//     }

//     for (char c : aminoAcids) {
//         generateSyntheticSet(kmer + c, k - 1, naturalSet, callback);
//     }
// }

// Main function
int main(int argc, char* argv[]) {
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(NULL);

	// Default values
	int kmerSize = 12;
	bool outputFrequenciesOnly = false;
	bool outputStatsAndMissingKmers = false;


	// Check command-line arguments
	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "-k") {
			if (i + 1 < argc) { // Make sure we aren't at the end of argv!
				kmerSize = std::stoi(argv[++i]); // Increment 'i' so we don't get the arguments confused.
			}
			else { // Uh-oh, there was no argument to the kmer option.
				std::cerr << "-k option requires one argument." << std::endl;
				return 1;
			}  
		}
		else if (arg == "-f") {
			outputFrequenciesOnly = true;
		}
		else if (arg == "-s") {
			outputStatsAndMissingKmers = true;
		}
		else if (arg == "--help" || arg == "-h") {
		    printHelp();
		    return 0;
		}
	}

	std::unordered_map<std::string, int> kmerCounts;

	std::string inputLine;
	while(getline(std::cin, inputLine) && !inputLine.empty()) {
		if(inputLine.size() >= kmerSize) {
			std::string kmer = inputLine.substr(0,kmerSize);
			kmerCounts[kmer]++;
		}
	}

	if (kmerCounts.empty()) {
	    printHelp();
	    return 1;
	}

	// Transfer the unordered_map to vector of pairs for sorting.
	std::vector<std::pair<std::string, int>> sortedKmers(kmerCounts.begin(), kmerCounts.end());

	std::sort(sortedKmers.begin(), sortedKmers.end(), 
		[](const std::pair<std::string, int> &a, const std::pair<std::string, int> &b) {
			return a.second > b.second;
		});

	for(const auto &kmer : sortedKmers) {
		if(outputFrequenciesOnly) {
			std::cout << kmer.second << std::endl;
		} else {
			std::cout << kmer.second << " " << kmer.first << std::endl;
		}
	}

	// Output statistics and missing kmers
	// if (outputStatsAndMissingKmers) {
	// 	cout << "Coverage: " << static_cast<double>(kmerCounts.size()) / pow(20, kmerSize) * 100 << "% (" << kmerCounts.size() << "/" << static_cast<int>(pow(20, kmerSize)) << ")\n";
	// 	generateSyntheticSet("", kmerSize, kmerCounts, [](const string& missingKmer) {
	// 		cout << missingKmer << "\n";
	// 	});
	// }

	return 0;
}
