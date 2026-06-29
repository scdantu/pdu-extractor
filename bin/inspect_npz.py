#!/usr/bin/env python3
"""Inspect NPZ file structure."""

import sys
import numpy as np
from pathlib import Path

if len(sys.argv) < 2:
    print("Usage: python3 inspect_npz.py <npz_file>")
    sys.exit(1)

npz_file = sys.argv[1]

if not Path(npz_file).exists():
    print(f"File not found: {npz_file}")
    sys.exit(1)

print(f"Inspecting: {npz_file}\n")

data = np.load(npz_file)

print("Keys in NPZ:")
for key in data.files:
    arr = data[key]
    print(f"  {key}: shape={arr.shape}, dtype={arr.dtype}, size={arr.nbytes / 1e6:.1f} MB")

print("\nFirst few values of each array:")
for key in data.files:
    arr = data[key]
    if len(arr) > 0:
        print(f"\n{key}[0:5]:")
        print(f"  {arr[0:5]}")
