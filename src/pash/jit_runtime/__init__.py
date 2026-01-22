"""
PaSh JIT Runtime Module

This module provides utilities for finding JIT runtime shell scripts.

The JIT runtime handles:
- Just-in-time compilation of shell regions
- Shell state management during parallel execution
- Communication with the compilation server
"""

from pathlib import Path


def get_jit_runtime_dir() -> Path:
    """Get the JIT runtime directory containing shell scripts."""
    return Path(__file__).parent


def get_script_path(name: str) -> Path:
    """Get the path to a JIT runtime script.

    Args:
        name: Name of the script (e.g., 'jit.sh')

    Returns:
        Path to the script

    Raises:
        FileNotFoundError: If the script doesn't exist
    """
    script = get_jit_runtime_dir() / name
    if not script.exists():
        raise FileNotFoundError(
            f"JIT runtime script '{name}' not found at {script}."
        )
    return script


# Convenience functions for each script
def get_jit_script() -> Path:
    """Get path to the main jit.sh script."""
    return get_script_path("jit.sh")


def get_prepare_call_compiler_script() -> Path:
    """Get path to pash_prepare_call_compiler.sh."""
    return get_script_path("pash_prepare_call_compiler.sh")


def get_restore_state_script() -> Path:
    """Get path to pash_restore_state_and_execute.sh."""
    return get_script_path("pash_restore_state_and_execute.sh")


def get_set_from_to_script() -> Path:
    """Get path to pash_set_from_to.sh."""
    return get_script_path("pash_set_from_to.sh")


def get_source_declare_vars_script() -> Path:
    """Get path to pash_source_declare_vars.sh."""
    return get_script_path("pash_source_declare_vars.sh")
