"""
DEPRECATED: This module exists for backward compatibility.
Use 'from compilation import ...' instead.

This module re-exports all symbols from compilation.py for backward compatibility.
It also maintains the main() function for CLI usage.
"""

import sys
import traceback

from sh_expand import env_vars_util

import config
from pash_graphviz import maybe_generate_graphviz
from util import log

from arg_parser import CompilerParser

# Re-export everything from compilation for backward compatibility
from compilation import (
    CompilerConfig,
    compile_ir,
    compile_optimize_output_script,
    load_df_region,
    compile_optimize_df_region,
    maybe_log_optimized_script,
    compile_candidate_df_region,
    optimize_irs,
    print_graph_statistics,
    choose_and_apply_parallelizing_transformations,
    choose_parallelizing_transformations,
    choose_parallelizing_transformation,
    apply_parallelizing_transformations,
    split_hdfs_cat_input,
    add_eager,
    add_eager_nodes,
)

# Import runtime_config from compilation (it's a module-level variable there)
import compilation
runtime_config = compilation.runtime_config


## We want to catch all exceptions here so that they are logged correctly
## and not just printed to the stderr.
def main():
    try:
        main_body()
    except Exception:
        log("Compiler failed, no need to worry, executing original script...")
        log(traceback.format_exc())
        sys.exit(1)


def main_body():
    ## Parse arguments
    args = parse_args()
    config.set_config_globals_from_pash_args(args)

    ## Load the configuration
    if not config.config:
        config.load_config(args.config_path)

    # Set runtime_config in compilation module
    compilation.runtime_config = config.config["distr_planner"]

    ## Read any shell variables files if present
    vars_dict = env_vars_util.read_vars_file(args.var_file, config.BASH_VERSION)
    config.set_vars_file(args.var_file, vars_dict)

    log("Input:", args.input_ir, "Compiled file:", args.compiled_script_file)

    ## Call the main procedure
    compiler_config = CompilerConfig(args.width)
    ast_or_ir = compile_optimize_output_script(
        args.input_ir, args.compiled_script_file, args, compiler_config
    )
    maybe_generate_graphviz(ast_or_ir, args)


def parse_args():
    parser = CompilerParser()
    args, _ = parser.parse_known_args()
    return args


if __name__ == "__main__":
    main()
