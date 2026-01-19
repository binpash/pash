from datetime import datetime
import os
import socket

from shell_ast import transformation_options, ast_to_ast
from parse import parse_shell_to_asts, from_ast_objects_to_shell
from util import *
from speculative import util_spec

LOGGING_PREFIX = "PaSh Preprocessor: "


@logging_prefix(LOGGING_PREFIX)
def preprocess(input_script_path, args):
    ## 1. Execute the POSIX shell parser that returns the AST in JSON
    preprocessing_parsing_start_time = datetime.now()
    ast_objects = parse_shell_to_asts(input_script_path, bash_mode=args.bash)
    preprocessing_parsing_end_time = datetime.now()
    print_time_delta(
        "Preprocessing -- Parsing",
        preprocessing_parsing_start_time,
        preprocessing_parsing_end_time,
    )

    ## 2. Preprocess ASTs by replacing possible candidates for compilation
    ##    with calls to the PaSh runtime.
    preprocessing_pash_start_time = datetime.now()
    preprocessed_asts = preprocess_asts(ast_objects, args)
    preprocessing_pash_end_time = datetime.now()
    print_time_delta(
        "Preprocessing -- PaSh",
        preprocessing_pash_start_time,
        preprocessing_pash_end_time,
    )

    ## 3. Translate the new AST back to shell syntax
    preprocessing_unparsing_start_time = datetime.now()
    preprocessed_shell_script = from_ast_objects_to_shell(preprocessed_asts)

    preprocessing_unparsing_end_time = datetime.now()
    print_time_delta(
        "Preprocessing -- Unparsing",
        preprocessing_unparsing_start_time,
        preprocessing_unparsing_end_time,
    )
    return preprocessed_shell_script


def preprocess_asts(ast_objects, args):
    trans_mode = transformation_options.TransformationType(args.preprocess_mode)
    if trans_mode is transformation_options.TransformationType.SPECULATIVE:
        trans_options = transformation_options.SpeculativeTransformationState(
            po_file=args.partial_order_file
        )
        util_spec.initialize(trans_options)
    elif trans_mode is transformation_options.TransformationType.AIRFLOW:
        trans_options = transformation_options.AirflowTransformationState()
    else:
        trans_options = transformation_options.TransformationState()

    ## Preprocess ASTs by replacing AST regions with calls to PaSh's runtime.
    ## Then the runtime will do the compilation and optimization with additional
    ## information.
    preprocessed_asts = ast_to_ast.replace_ast_regions(ast_objects, trans_options)

    ## Let the scheduler know that we are done with the partial_order file
    ## TODO: We could stream the partial_order_file to the scheduler
    if trans_mode is transformation_options.TransformationType.SPECULATIVE:
        ## First complete the partial_order file
        util_spec.serialize_partial_order(trans_options)

        ## Then inform the scheduler that it can read it
        unix_socket_file = os.getenv("PASH_SPEC_SCHEDULER_SOCKET")
        msg = util_spec.scheduler_server_init_po_msg(
            trans_options.get_partial_order_file()
        )
        _unix_socket_send_and_forget(unix_socket_file, msg)

    return preprocessed_asts


def _unix_socket_send_and_forget(socket_file, msg):
    """Send a message to a Unix socket without waiting for a response."""
    if socket_file is None:
        return
    try:
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        sock.connect(socket_file)
        sock.sendall(msg.encode())
        sock.close()
    except Exception as e:
        log(f"Warning: Failed to send message to scheduler socket: {e}")
