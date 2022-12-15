import argparse
from datetime import datetime
import os

import config
from ast_to_ir import replace_ast_regions
from ir import FileIdGen
from parse import parse_shell_to_asts, from_ast_objects_to_shell
from util import *

LOGGING_PREFIX = "PaSh Preprocessor: "

@logging_prefix(LOGGING_PREFIX)
def preprocess(input_script_path, args):
    print(input_script_path)
    ## 1. Execute the POSIX shell parser that returns the AST in JSON
    preprocessing_parsing_start_time = datetime.now()
    ast_objects = parse_shell_to_asts(input_script_path)
    preprocessing_parsing_end_time = datetime.now()
    print_time_delta("Preprocessing -- Parsing", preprocessing_parsing_start_time, preprocessing_parsing_end_time, args)

    ## 2. Preprocess ASTs by replacing possible candidates for compilation
    ##    with calls to the PaSh runtime.
    preprocessing_pash_start_time = datetime.now()
    preprocessed_asts = preprocess_asts(ast_objects)
    preprocessing_pash_end_time = datetime.now()
    print_time_delta("Preprocessing -- PaSh", preprocessing_pash_start_time, preprocessing_pash_end_time, args)

    ## 3. Translate the new AST back to shell syntax
    preprocessing_unparsing_start_time = datetime.now()
    preprocessed_shell_script = from_ast_objects_to_shell(preprocessed_asts)
    if(args.output_preprocessed):
        log("Preprocessed script:")
        log(preprocessed_shell_script)
    
    preprocessing_unparsing_end_time = datetime.now()
    print_time_delta("Preprocessing -- Unparsing", preprocessing_unparsing_start_time, preprocessing_unparsing_end_time, args)
    return preprocessed_shell_script


def preprocess_asts(ast_objects):

    ## TODO: Add a potential selection here to decide what kind of transformation to apply

    ## Preprocess ASTs by replacing AST regions with calls to PaSh's runtime.
    ## Then the runtime will do the compilation and optimization with additional
    ## information.
    preprocessed_asts = replace_ast_regions(ast_objects)

    return preprocessed_asts

##
## This is the command line interface for the preprocessor
##
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the script to be preprocessed")
    parser.add_argument("--output_preprocessed",
                        help=" output the preprocessed script",
                        action="store_true")
    config.add_general_config_arguments(parser)

    args = parser.parse_args()
    config.set_config_globals_from_pash_args(args)

    ## Initialize the log file
    ## TODO: Can we move this somewhere where there is no need for copy paste?
    config.init_log_file()
    preprocess(args.input, args)

if __name__ == '__main__':
    main()
