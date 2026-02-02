#!/usr/bin/env python3
"""
PaSh CLI Entry Point

This module provides the `pash` command-line interface as a wrapper
around pa.sh. When installed via pip, users can run `pash` instead
of `./pa.sh`.

The actual execution is delegated to pa.sh which handles:
1. Environment setup (shell functions, FIFOs, etc.)
2. Starting the compilation server
3. Running the preprocessor
4. Executing the preprocessed script
5. Cleanup
"""

import os
import subprocess
import sys
from pathlib import Path

from . import __version__


def get_pash_top() -> Path:
    """Get the PASH_TOP directory (package installation path)."""
    return Path(__file__).parent


def get_pa_sh_path() -> Path:
    """Get the path to pa.sh.

    pa.sh is located in the same directory as this file (src/pash/).
    """
    pash_top = get_pash_top()
    pa_sh = pash_top / "pa.sh"

    if pa_sh.exists():
        return pa_sh

    # Fallback: check PASH_TOP environment variable
    pash_top_env = os.environ.get("PASH_TOP")
    if pash_top_env:
        env_pa_sh = Path(pash_top_env) / "pa.sh"
        if env_pa_sh.exists():
            return env_pa_sh

    raise FileNotFoundError(
        "Could not find pa.sh. Make sure PaSh is properly installed."
    )


def main():
    """Main entry point for PaSh CLI.

    This delegates to pa.sh for actual execution, passing through all arguments.
    """
    # Handle --version specially
    if len(sys.argv) == 2 and sys.argv[1] in ("--version", "-V"):
        print(f"pash {__version__}")
        sys.exit(0)

    try:
        pa_sh_path = get_pa_sh_path()
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Set PASH_TOP to point to the src/pash directory
    env = os.environ.copy()
    env["PASH_TOP"] = str(get_pash_top())

    # Pass all arguments to pa.sh
    args = ["bash", str(pa_sh_path)] + sys.argv[1:]

    try:
        result = subprocess.run(args, env=env)
        sys.exit(result.returncode)
    except KeyboardInterrupt:
        sys.exit(130)
    except Exception as e:
        print(f"Error running PaSh: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
