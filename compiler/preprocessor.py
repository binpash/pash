import argparse
from datetime import datetime
import os

import config
import ast_to_ast
from ir import FileIdGen
from parse import parse_shell_to_asts, from_ast_objects_to_shell
from util import *

LOGGING_PREFIX = "PaSh Preprocessor: "

@logging_prefix(LOGGING_PREFIX)
def preprocess(input_script_path, args, mode):
    ## 1. Execute the POSIX shell parser that returns the AST in JSON
    preprocessing_parsing_start_time = datetime.now()
    ast_objects = parse_shell_to_asts(input_script_path)
    preprocessing_parsing_end_time = datetime.now()
    print_time_delta("Preprocessing -- Parsing", preprocessing_parsing_start_time, preprocessing_parsing_end_time)

    ## 2. Preprocess ASTs by replacing possible candidates for compilation
    ##    with calls to the PaSh runtime.
    preprocessing_pash_start_time = datetime.now()
    preprocessed_asts = preprocess_asts(ast_objects, mode)
    preprocessing_pash_end_time = datetime.now()
    print_time_delta("Preprocessing -- PaSh", preprocessing_pash_start_time, preprocessing_pash_end_time)

    ## 3. Translate the new AST back to shell syntax
    preprocessing_unparsing_start_time = datetime.now()
    preprocessed_shell_script = from_ast_objects_to_shell(preprocessed_asts)
    
    preprocessing_unparsing_end_time = datetime.now()
    print_time_delta("Preprocessing -- Unparsing", preprocessing_unparsing_start_time, preprocessing_unparsing_end_time)
    return preprocessed_shell_script


def preprocess_asts(ast_objects, mode):

    ## TODO: Add a potential selection here to decide what kind of transformation to apply

    ## Preprocess ASTs by replacing AST regions with calls to PaSh's runtime.
    ## Then the runtime will do the compilation and optimization with additional
    ## information.
    preprocessed_asts = ast_to_ast.replace_ast_regions(ast_objects, mode)

    return preprocessed_asts

##
## This is the command line interface for the preprocessor
##
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the script to be preprocessed")
    config.add_general_config_arguments(parser)

    args = parser.parse_args()
    config.set_config_globals_from_pash_args(args)

    ## Initialize the log file
    ## TODO: Can we move this somewhere where there is no need for copy paste?
    config.init_log_file()
    mode = ast_to_ast.TransformationType('spec')
    preprocessed_shell_script = preprocess(args.input, args, mode)
    print(preprocessed_shell_script)

if __name__ == '__main__':
    main()
