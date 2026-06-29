#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

PDB_DIR="${PDB_DIR:-pdb}"
PDU_DB="${PDU_DB:-pdu_output/pdus.sqlite}"
RADIUS_ANGSTROM="${RADIUS_ANGSTROM:-15.0}"
ANNOTATION_CSV="${ANNOTATION_CSV:-}"
HANDLE_ALL_PDBS="${HANDLE_ALL_PDBS:-true}"
RESET_OUTPUT="${RESET_OUTPUT:-false}"
LOG_DIR="${LOG_DIR:-logs}"
LOG_FILE="${LOG_FILE:-${LOG_DIR}/extract_pdus.log}"
LOG_LEVEL="${LOG_LEVEL:-INFO}"

mkdir -p "${LOG_DIR}"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "$(date '+%Y-%m-%d %H:%M:%S') INFO extract_pdus: starting"
echo "$(date '+%Y-%m-%d %H:%M:%S') INFO extract_pdus: PDB_DIR=${PDB_DIR} PDU_DB=${PDU_DB} RADIUS_ANGSTROM=${RADIUS_ANGSTROM}"

if ! find "${PDB_DIR}" -name '*.ent.gz' -print -quit | grep -q .; then
    echo "No .ent.gz PDB files found under ${PDB_DIR}."
    echo "Download a subset with: python scripts/download_pdb_subset.py --count 1000 --output-dir pdb"
    exit 1
fi

if [[ ! -x "bin/extract_pdb_coordinates" ]]; then
    scripts/buildcpp.sh
fi

if [[ "${RESET_OUTPUT}" == "true" ]]; then
    rm -rf pdb_output "${PDU_DB}"
fi

mkdir -p "$(dirname "${PDU_DB}")"
rm -rf pdb_output

cmd=(
    python kmers/pipeline.py
    --handle_all_pdbs "${HANDLE_ALL_PDBS}"
    --pdu_db "${PDU_DB}"
    --pdu_radius_angstrom "${RADIUS_ANGSTROM}"
    --log-file "${LOG_DIR}/pipeline.log"
    --log-level "${LOG_LEVEL}"
)

if [[ -n "${ANNOTATION_CSV}" ]]; then
    cmd+=(--annotation_csv "${ANNOTATION_CSV}")
fi

PYTHONPATH="${REPO_DIR}:${PYTHONPATH:-}" "${cmd[@]}"

echo "PDU database written to ${PDU_DB}"
