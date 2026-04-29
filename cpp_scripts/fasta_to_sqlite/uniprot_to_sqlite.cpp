#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <sqlite3.h>

#define CHUNK_SIZE 10000

struct Record {
    std::string id;
    std::string sequence;
};

void chunked_insert(sqlite3* db, std::vector<Record>& records) {
    char* zErrMsg = 0;
    sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &zErrMsg);
    
    for(auto& record : records){
        std::string sql = "INSERT INTO sequences (id, sequence) VALUES ('" + record.id + "', '" + record.sequence + "');";
        sqlite3_exec(db, sql.c_str(), NULL, NULL, &zErrMsg);
    }
    
    sqlite3_exec(db, "COMMIT;", NULL, NULL, &zErrMsg);
    records.clear();
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    sqlite3* db;
    sqlite3_open("uniprotkb/uniprot_sequences.db", &db);

    char* zErrMsg = 0;
    sqlite3_exec(db, "CREATE TABLE IF NOT EXISTS sequences (id TEXT PRIMARY KEY, sequence TEXT NOT NULL);", NULL, NULL, &zErrMsg);
    // sqlite3_exec(db, "CREATE INDEX IF NOT EXISTS idx_id ON sequences (id);", NULL, NULL, &zErrMsg);
    // sqlite3_exec(db, "CREATE INDEX IF NOT EXISTS idx_sequence ON sequences (sequence);", NULL, NULL, &zErrMsg);
    
    std::vector<Record> records;
    std::string line, uniprot_id, sequence;
    while(std::getline(std::cin, line)) {
        if(line[0] == '>') {
            if(!uniprot_id.empty() && !sequence.empty()) {
                records.push_back({uniprot_id, sequence});
                sequence.clear();
                
                if(records.size() >= CHUNK_SIZE){
                    chunked_insert(db, records);
                }
            }

            // Extracting Uniprot ID
            int first_pipe = line.find("|");
            int second_pipe = line.find("|", first_pipe + 1);
            uniprot_id = line.substr(first_pipe + 1, second_pipe - first_pipe - 1);
        } else {
            sequence += line;
        }
    }
    
    if(!uniprot_id.empty() && !sequence.empty()){
        records.push_back({uniprot_id, sequence});
    }
    
    if(!records.empty()){
        chunked_insert(db, records);
    }
    
    sqlite3_close(db);
    return 0;
}
