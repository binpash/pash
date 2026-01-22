"""
PaSh Compiler Package

Entry point:
- compilation_server.py: Compilation daemon/server (started by pa.sh)

Core API:
- compilation.compile_ir(): Main compilation function
- compilation.CompilerConfig: Compilation configuration

Argument Parsing:
- arg_parser.BaseParser: Base argument parser
- arg_parser.CompilerParser: Full compiler argument parser
"""

from .compilation import (
    compile_ir,
    CompilerConfig,
    compile_optimize_df_region,
    compile_optimize_output_script,
    load_df_region,
    optimize_irs,
)
from .ir import IR, FileIdGen

__all__ = [
    'compile_ir',
    'CompilerConfig',
    'compile_optimize_df_region',
    'compile_optimize_output_script',
    'load_df_region',
    'optimize_irs',
    'IR',
    'FileIdGen',
]
