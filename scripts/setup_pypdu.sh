#!/usr/bin/env bash
set -euo pipefail

VENV_DIR="${1:-pypdu}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

"${PYTHON_BIN}" -m venv "${VENV_DIR}"
"${VENV_DIR}/bin/python" -m pip install --upgrade pip setuptools wheel
"${VENV_DIR}/bin/python" -m pip install -r requirements.txt

mkdir -p .cache/matplotlib
"${VENV_DIR}/bin/python" scripts/check_environment.py
