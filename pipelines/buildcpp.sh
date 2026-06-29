#!/bin/bash

# script for compiling binaries

# Check if g++ is installed
command -v g++ >/dev/null 2>&1 || { echo >&2 "g++ required but it's not installed. Aborting."; exit 1; }

# Check if sqlite3 is installed
command -v sqlite3 >/dev/null 2>&1 || { echo >&2 "sqlite3 required but it's not installed (apt/brew install sqlite3 [on debian/macos]). Aborting."; exit 1; }

# Compile C++ program
echo "Compiling C++ binaries."

mkdir -p bin

arch -x86_64 g++ -std=c++17 -o "bin/fasta_to_sqlite" cpp_scripts/fasta_to_sqlite/*.cpp -lsqlite3
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

g++ -std=c++17 -o "bin/post_process_kmers" cpp_scripts/post_process_kmers/*.cpp
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

g++ -std=c++17 -o "bin/extract_pdb_coordinates" cpp_scripts/extract_pdb_coordinates/*.cpp
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

echo "Done compiling binaries."
