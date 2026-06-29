#!/usr/bin/env python3
import importlib
import os
import platform
import sys
from pathlib import Path

from kmers.logging_utils import add_logging_args, configure_logging


REQUIRED_MODULES = [
    "numpy",
    "pandas",
    "sklearn",
    "scipy",
    "umap",
    "hdbscan",
    "matplotlib",
    "seaborn",
    "tqdm",
    "torch",
]


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Check the pypdu Python environment.")
    add_logging_args(parser)
    args = parser.parse_args()
    configure_logging(args.log_file, args.log_level)
    configure_local_caches()
    print(f"Python: {sys.version.split()[0]}")
    print(f"Platform: {platform.platform()}")
    missing = []
    for module in REQUIRED_MODULES:
        try:
            imported = importlib.import_module(module)
            version = getattr(imported, "__version__", "unknown")
            print(f"{module}: {version}")
        except Exception as exc:
            missing.append((module, exc))
            print(f"{module}: missing ({exc})")

    if missing:
        raise SystemExit(1)

    import torch

    if torch.cuda.is_available():
        device = "cuda"
        detail = torch.cuda.get_device_name(0)
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        device = "mps"
        detail = "Apple Metal Performance Shaders"
    else:
        device = "cpu"
        detail = "CPU"

    print(f"Selected PyTorch device: {device} ({detail})")


def configure_local_caches():
    cache_dir = Path(".cache")
    cache_dir.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))


if __name__ == "__main__":
    main()
