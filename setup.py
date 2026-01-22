"""
Setup script for PaSh with C binary compilation hooks.

This script handles the compilation of C runtime binaries during pip install.
The actual package configuration is in pyproject.toml.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop


def get_runtime_src_dir():
    """Get the runtime source directory."""
    return Path(__file__).parent / "src" / "pash" / "runtime"


def compile_runtime_binaries(build_lib=None):
    """Compile C runtime binaries using make.

    Args:
        build_lib: If provided, copy compiled binaries to this directory.
                   If None, compile in-place (for development mode).
    """
    runtime_src = get_runtime_src_dir()

    if not runtime_src.exists():
        print(f"Warning: Runtime source directory not found: {runtime_src}", file=sys.stderr)
        return

    # Check for required tools
    for tool in ["gcc", "make", "git"]:
        if shutil.which(tool) is None:
            raise RuntimeError(
                f"Required tool '{tool}' not found. "
                f"Please install build tools:\n"
                f"  Ubuntu/Debian: sudo apt install build-essential git\n"
                f"  Fedora/RHEL: sudo dnf install gcc make git\n"
                f"  macOS: xcode-select --install && brew install git"
            )

    # Run make in the runtime directory
    print(f"Compiling PaSh runtime binaries in {runtime_src}...")

    # Binary targets (including dgsh-tee which requires git to clone during build)
    targets = ["eager", "split", "r-merge", "r-wrap", "r-split", "r-unwrap", "set-diff", "dgsh-tee"]

    try:
        # Clean first to ensure fresh build
        subprocess.run(["make", "clean"], cwd=runtime_src, check=False, capture_output=True)

        # Compile
        result = subprocess.run(
            ["make"] + targets,
            cwd=runtime_src,
            check=True,
            capture_output=True,
            text=True,
        )
        if result.stdout:
            print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Make failed with output:\n{e.stdout}\n{e.stderr}", file=sys.stderr)
        raise RuntimeError(f"Failed to compile runtime binaries: {e}")

    # Copy binaries to target if different from source
    if build_lib:
        target_dir = Path(build_lib) / "pash" / "runtime"
        target_dir.mkdir(parents=True, exist_ok=True)

        # Binary names (note: make targets use hyphens, binaries use underscores)
        binaries = ["eager", "split", "r_split", "r_merge", "r_wrap", "r_unwrap", "set-diff", "dgsh-tee"]

        for binary in binaries:
            src_binary = runtime_src / binary
            if src_binary.exists():
                shutil.copy2(src_binary, target_dir / binary)
                print(f"  Copied {binary} to {target_dir}")
            else:
                print(f"  Warning: Binary {binary} not found at {src_binary}", file=sys.stderr)


class BuildPyWithRuntime(build_py):
    """Custom build_py that compiles C runtime binaries."""

    def run(self):
        # Run standard build_py first
        build_py.run(self)
        # Then compile runtime binaries
        compile_runtime_binaries(self.build_lib)


class DevelopWithRuntime(develop):
    """Custom develop that compiles C runtime binaries in place."""

    def run(self):
        # Compile runtime binaries in-place for development
        compile_runtime_binaries(None)
        # Run standard develop
        develop.run(self)


setup(
    cmdclass={
        "build_py": BuildPyWithRuntime,
        "develop": DevelopWithRuntime,
        "install": InstallWithRuntime,
    },
    # Don't mark as having ext_modules - we compile from source at install time
    # This creates a pure Python wheel that works on any platform
)
