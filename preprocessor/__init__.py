"""
PaSh Preprocessor Package

This package handles preprocessing of shell scripts for PaSh parallelization.
It identifies candidate dataflow regions and replaces them with calls to the
PaSh runtime for compilation and optimization.
"""

from preprocessor import preprocess

__all__ = ["preprocess"]
