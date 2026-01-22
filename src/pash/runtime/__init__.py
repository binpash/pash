"""
PaSh Runtime Module

This module provides utilities for finding compiled runtime binaries.
The binaries are compiled from C source during pip install.

Available binaries:
- eager: Eager evaluation node
- split: Data splitting
- r_split: Record-aware splitting
- r_merge: Record-aware merging
- r_wrap: Record wrapping
- r_unwrap: Record unwrapping
- set-diff: Set difference
"""

from pathlib import Path


def get_runtime_dir() -> Path:
    """Get the runtime directory containing binaries."""
    return Path(__file__).parent


def get_binary_path(name: str) -> Path:
    """Get the path to a runtime binary.

    Args:
        name: Name of the binary (e.g., 'eager', 'split')

    Returns:
        Path to the binary

    Raises:
        FileNotFoundError: If the binary doesn't exist
    """
    binary = get_runtime_dir() / name
    if not binary.exists():
        raise FileNotFoundError(
            f"Runtime binary '{name}' not found at {binary}. "
            f"You may need to reinstall pash with 'pip install --force-reinstall pash'"
        )
    return binary


# Convenience functions for each binary
def get_eager_binary() -> Path:
    """Get path to the eager binary."""
    return get_binary_path("eager")


def get_split_binary() -> Path:
    """Get path to the split binary."""
    return get_binary_path("split")


def get_r_split_binary() -> Path:
    """Get path to the r_split binary."""
    return get_binary_path("r_split")


def get_r_merge_binary() -> Path:
    """Get path to the r_merge binary."""
    return get_binary_path("r_merge")


def get_r_wrap_binary() -> Path:
    """Get path to the r_wrap binary."""
    return get_binary_path("r_wrap")


def get_r_unwrap_binary() -> Path:
    """Get path to the r_unwrap binary."""
    return get_binary_path("r_unwrap")


def get_set_diff_binary() -> Path:
    """Get path to the set-diff binary."""
    return get_binary_path("set-diff")
