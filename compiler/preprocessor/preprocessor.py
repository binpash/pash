import argparse
from datetime import datetime
import os

import config
from shell_ast import ast_to_ast
from ir import FileIdGen
from parse import parse_shell_to_asts, from_ast_objects_to_shell
from util import *
import server_util
from speculative import util_spec

LOGGING_PREFIX = "PaSh Preprocessor: "

@logging_prefix(LOGGING_PREFIX)
def preprocess(input_script_path, args):
    ## 1. Execute the POSIX shell parser that returns the AST in JSON
    preprocessing_parsing_start_time = datetime.now()
    ast_objects = parse_shell_to_asts(input_script_path)
    preprocessing_parsing_end_time = datetime.now()
    print_time_delta("Preprocessing -- Parsing", preprocessing_parsing_start_time, preprocessing_parsing_end_time)

    ## 2. Preprocess ASTs by replacing possible candidates for compilation
    ##    with calls to the PaSh runtime.
    preprocessing_pash_start_time = datetime.now()
    preprocessed_asts = preprocess_asts(ast_objects, args)
    preprocessing_pash_end_time = datetime.now()
    print_time_delta("Preprocessing -- PaSh", preprocessing_pash_start_time, preprocessing_pash_end_time)

    ## 3. Translate the new AST back to shell syntax
    preprocessing_unparsing_start_time = datetime.now()
    preprocessed_shell_script = from_ast_objects_to_shell(preprocessed_asts)
    
    preprocessing_unparsing_end_time = datetime.now()
    print_time_delta("Preprocessing -- Unparsing", preprocessing_unparsing_start_time, preprocessing_unparsing_end_time)
    return preprocessed_shell_script


def preprocess_asts(ast_objects, args):
    trans_mode = ast_to_ast.TransformationType(args.preprocess_mode)
    if trans_mode is ast_to_ast.TransformationType.SPECULATIVE:
        trans_options = ast_to_ast.SpeculativeTransformationState(mode=trans_mode,
                                                                  po_file=args.partial_order_file)
        util_spec.initialize(trans_options)
    else:
        trans_options = ast_to_ast.TransformationState(mode=trans_mode)

    ## Preprocess ASTs by replacing AST regions with calls to PaSh's runtime.
    ## Then the runtime will do the compilation and optimization with additional
    ## information.
    preprocessed_asts = ast_to_ast.replace_ast_regions(ast_objects, trans_options)

    ## Let the scheduler know that we are done with the partial_order file
    ## TODO: We could stream the partial_order_file to the scheduler
    if trans_mode is ast_to_ast.TransformationType.SPECULATIVE:
        ## First complete the partial_order file
        util_spec.serialize_partial_order(trans_options)

        ## Then inform the scheduler that it can read it
        unix_socket_file = os.getenv("PASH_SPEC_SCHEDULER_SOCKET")
        msg = util_spec.scheduler_server_init_po_msg(trans_options.get_partial_order_file())
        server_util.unix_socket_send_and_forget(unix_socket_file, msg)

    return preprocessed_asts

##
## This is the command line interface for the preprocessor
##
def main():
    parser = argparse.ArgumentParser()
    config.add_general_config_arguments(parser)

    subparsers = parser.add_subparsers(help='sub-command help')

    # create the parser for the "a" command
    parser_pash = subparsers.add_parser('pash', help='Preprocess the script so that it can be run with PaSh')
    config.add_common_arguments(parser_pash)
    parser_pash.add_argument("input", help="the script to be preprocessed")
    parser_pash.set_defaults(preprocess_mode='pash')

    # create the parser for the "b" command
    parser_spec = subparsers.add_parser('spec', help='Preprocess the script so that it can be run with speculation')
    parser_spec.add_argument("input", help="the script to be preprocessed")

    ## TODO: When we better integrate, this should be automatically set.
    parser_spec.add_argument("partial_order_file", help="the file to store the partial order (currently just a sequence)")
    parser_spec.set_defaults(preprocess_mode='spec')
    
    args = parser.parse_args()
    config.set_config_globals_from_pash_args(args)

    ## Initialize the log file
    ## TODO: Can we move this somewhere where there is no need for copy paste?
    config.init_log_file()

    ## Choose the transformation node based on the argument given
    log("Preprocesing mode:", args.preprocess_mode)
    preprocessed_shell_script = preprocess(args.input, args)
    print(preprocessed_shell_script)

if __name__ == '__main__':
    main()
