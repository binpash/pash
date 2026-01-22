"""
Configuration module for PaSh compiler.

Constants (set at import time from environment or package location):
- PASH_TOP: Root directory of PaSh installation/package
- PASH_TMP_PREFIX: Temporary directory for this session
- BASH_VERSION: Tuple of bash version numbers
- HDFS_PREFIX: Prefix for HDFS paths

Runtime state (set during initialization):
- config: Dict loaded from config.json
- pash_args: Parsed command-line arguments
- LOGGING_PREFIX: Prefix for log messages
"""

import json
import logging
import os
import sys
from pathlib import Path


def get_package_root() -> Path:
    """Get the root directory of the pash package.

    This is the directory containing compiler/, preprocessor/, etc.
    When installed via pip, this is the site-packages/pash directory.
    """
    return Path(__file__).parent.parent


def get_pash_top() -> str:
    """Get PASH_TOP, preferring environment variable if set.

    This allows both pip-installed usage (uses package location)
    and development usage (uses PASH_TOP from pa.sh).
    """
    env_pash_top = os.getenv("PASH_TOP")
    if env_pash_top:
        return env_pash_top
    return str(get_package_root())


# === Constants from environment or package location ===

PASH_TOP = get_pash_top()

# PASH_TMP_PREFIX must be set by the entry point (pa.sh or cli.py)
_pash_tmp_prefix = os.getenv("PASH_TMP_PREFIX")
if _pash_tmp_prefix is None:
    # Fallback for when module is imported outside of normal execution
    import tempfile
    _pash_tmp_prefix = tempfile.mkdtemp(prefix="pash_") + "/"
PASH_TMP_PREFIX = _pash_tmp_prefix

# BASH_VERSION with fallback
_bash_version_str = os.getenv("PASH_BASH_VERSION", "5 0 0")
BASH_VERSION = tuple(int(i) for i in _bash_version_str.split(" "))

HDFS_PREFIX = "$HDFS_DATANODE_DIR/"

# === Runtime state ===

LOGGING_PREFIX = ""
config = {}
pash_args = None

# Increase recursion limit for parser/unparser
sys.setrecursionlimit(10000)


def set_config_globals_from_pash_args(given_pash_args):
    """Initialize runtime state from parsed arguments."""
    global pash_args
    pash_args = given_pash_args

    # Configure logging
    if given_pash_args.log_file == "":
        logging.basicConfig(format="%(message)s")
    else:
        logging.basicConfig(
            format="%(message)s",
            filename=os.path.abspath(given_pash_args.log_file),
            filemode="w",
        )

    # Set debug level
    if given_pash_args.debug == 0:
        logging.getLogger().setLevel(logging.ERROR)
    elif given_pash_args.debug == 1:
        logging.getLogger().setLevel(logging.WARNING)
    elif given_pash_args.debug == 2:
        logging.getLogger().setLevel(logging.INFO)
    elif given_pash_args.debug >= 3:
        logging.getLogger().setLevel(logging.DEBUG)


def load_config(config_file_path=""):
    """Load configuration from JSON file."""
    global config

    if config_file_path == "":
        # Look for config.json relative to this module (in compiler directory)
        config_file_path = str(Path(__file__).parent / "config.json")

    with open(config_file_path) as config_file:
        config = json.load(config_file)

    if not config:
        raise Exception(f"No valid configuration could be loaded from {config_file_path}")
    if "distr_planner" not in config:
        raise Exception(f"Missing `distr_planner` config in {config_file_path}")


def set_vars_file(var_file_path: str, var_dict: dict):
    """Set shell variables in config."""
    global config
    config["shell_variables"] = var_dict
    config["shell_variables_file_path"] = var_file_path
