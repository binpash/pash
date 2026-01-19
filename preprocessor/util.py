"""
Utility functions for the PaSh preprocessor.

This is a simplified version of compiler/util.py that doesn't depend on config.py.
Configuration values are read from environment variables.
"""

from datetime import timedelta
import functools
import logging
import os
import tempfile

# Configuration from environment variables (set by pa.sh or pash_runtime.sh)
PASH_TMP_PREFIX = os.environ.get("PASH_TMP_PREFIX", "/tmp/pash_tmp/")
OUTPUT_TIME = os.environ.get("pash_output_time_flag", "1") == "1"
LOGGING_PREFIX = "PaSh: "


def unzip(lst):
    """Unzip a list of pairs into two separate lists."""
    res = [[i for i, j in lst], [j for i, j in lst]]
    return res


def print_time_delta(prefix, start_time, end_time):
    """Output timing information to the log."""
    time_difference = (end_time - start_time) / timedelta(milliseconds=1)
    if OUTPUT_TIME:
        log("{} time:".format(prefix), time_difference, " ms", level=0)
    else:
        log("{} time:".format(prefix), time_difference, " ms")


def logging_prefix(logging_prefix_str):
    """Decorator to add logging prefix to a function."""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            global LOGGING_PREFIX
            old_prefix = LOGGING_PREFIX
            LOGGING_PREFIX = logging_prefix_str
            result = func(*args, **kwargs)
            LOGGING_PREFIX = old_prefix
            return result
        return wrapper
    return decorator


def log(*args, end="\n", level=1):
    """Wrapper for logging."""
    if level >= 1:
        concatted_args = " ".join([str(a) for a in list(args)])
        logging.info(f"{LOGGING_PREFIX} {concatted_args}")


def ptempfile():
    """Create a temporary file in the PaSh temp directory."""
    fd, name = tempfile.mkstemp(dir=PASH_TMP_PREFIX)
    os.close(fd)
    return name


def make_kv(key, val):
    """Make a key-value pair in AST JSON format."""
    return [key, val]
