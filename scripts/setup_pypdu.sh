#!/usr/bin/env bash
set -euo pipefail

VENV_DIR="${1:-pypdu}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
LOG_DIR="${LOG_DIR:-logs}"
LOG_FILE="${LOG_FILE:-${LOG_DIR}/setup_pypdu.log}"

mkdir -p "${LOG_DIR}"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "$(date '+%Y-%m-%d %H:%M:%S') INFO setup_pypdu: creating ${VENV_DIR} with ${PYTHON_BIN}"

"${PYTHON_BIN}" -m venv "${VENV_DIR}"
"${VENV_DIR}/bin/python" -m pip install --upgrade pip setuptools wheel
"${VENV_DIR}/bin/python" -m pip install -r requirements.txt

mkdir -p .cache/matplotlib
"${VENV_DIR}/bin/python" scripts/check_environment.py --log-file "${LOG_DIR}/check_environment.log"
