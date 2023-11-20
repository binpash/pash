import os
from datetime import datetime

from util import *
from shell_ast.ast_util import *
from parse import from_ast_objects_to_shell
import config

RM_PASH_FIFOS_NAME = "rm_pash_fifos"
MKFIFO_PASH_FIFOS_NAME = "mkfifo_pash_fifos"


def to_shell(ir, args):
    backend_start_time = datetime.now()

    ## First call an IR to AST compilation pass
    output_asts = ir2ast(ir, args)

    ## Then just call the parser.
    output_script = from_ast_objects_to_shell(output_asts)

    backend_end_time = datetime.now()
    print_time_delta("Backend", backend_start_time, backend_end_time)

    return output_script


def ir2ast(ir, args):
    clean_up_graph = False
    drain_streams = False
    if args.termination == "clean_up_graph":
        clean_up_graph = True
    elif args.termination == "drain_stream":
        drain_streams = True

    ## NOTE: We first need to make the main body because it might create additional ephemeral fids.

    ## TODO: If we have just a single node maybe we should just instantiate that without anything else.
    ## This is implemented below, but there is no need to turn it on now.
    ##
    ## If we only have a single node, then we don't need to make prologues, epilogues, etc
    # if len(ir.nodes) == 1:
    #     nodes = list(ir.nodes.values())
    #     assert(len(nodes) == 1)
    #     node = nodes[0]
    #     body = [node.to_ast(ir.edges, drain_streams)]
    #     return body

    ## Make the main body
    body = ir.to_ast(drain_streams)

    all_fids = ir.all_fids()

    # log("All fids:", all_fids)
    ## Find all the ephemeral fids and turn them to ASTs
    ephemeral_fids = [fid for fid in all_fids if fid.is_ephemeral()]

    # log("Ephemeral fids:", ephemeral_fids)

    ## Call the prologue that creates fifos for all ephemeral fids
    prologue = make_ir_prologue(ephemeral_fids)

    ## Call the epilogue that removes all ephemeral fids
    epilogue = make_ir_epilogue(ephemeral_fids, clean_up_graph, args.log_file)

    final_asts = prologue + body + epilogue

    return final_asts


def make_rms_f_prologue_epilogue(ephemeral_fids):
    asts = []
    ## Create an `rm -f` for each ephemeral fid
    for eph_fid in ephemeral_fids:
        args = [eph_fid.to_ast()]
        command = make_rm_f_ast(args)
        asts.append(command)
    return asts


def make_ir_prologue(ephemeral_fids) -> "list[AstNode]":
    asts = []
    ## Create an `rm -f` for each ephemeral fid
    rm_asts = make_rms_f_prologue_epilogue(ephemeral_fids)
    defun_rm_pash_fifos = make_defun(RM_PASH_FIFOS_NAME, make_semi_sequence(rm_asts))
    asts.append(defun_rm_pash_fifos)

    ## Create a `mkfifo` for each ephemeral fid
    mkfifo_asts = []
    for eph_fid in ephemeral_fids:
        args = [eph_fid.to_ast()]
        command = make_mkfifo_ast(args)
        mkfifo_asts.append(command)

    defun_mkfifos = make_defun(MKFIFO_PASH_FIFOS_NAME, make_semi_sequence(mkfifo_asts))
    asts.append(defun_mkfifos)

    call_rm_pash_fifos = make_command([string_to_argument(RM_PASH_FIFOS_NAME)])
    asts.append(call_rm_pash_fifos)

    call_mkfifos = make_command([string_to_argument(MKFIFO_PASH_FIFOS_NAME)])
    asts.append(call_mkfifos)

    class_asts = [to_ast_node(ast) for ast in asts]
    return class_asts


def make_ir_epilogue(ephemeral_fids, clean_up_graph, log_file) -> "list[AstNode]":
    asts = []
    if clean_up_graph:
        ## TODO: Wait for all output nodes not just one
        pids = [[standard_var_ast("!")]]
        clean_up_path_script = os.path.join(
            config.PASH_TOP, config.config["runtime"]["clean_up_graph_binary"]
        )
        com_args = [
            string_to_argument("source"),
            string_to_argument(clean_up_path_script),
        ] + pids
        if log_file == "":
            com = make_command(com_args)
        else:
            redirection = redir_append_stderr_to_string_file(log_file)
            com = make_command(com_args, redirections=[redirection])
        asts.append(com)
    else:
        ## Otherwise we just wait for all processes to die.
        wait_com = make_command([string_to_argument("wait")])
        exit_status = make_command([string_to_argument("internal_exec_status=$?")])
        asts.extend([wait_com, exit_status])

    ## Create an `rm -f` for each ephemeral fid
    call_rm_pash_funs = make_command([string_to_argument(RM_PASH_FIFOS_NAME)])
    asts.append(call_rm_pash_funs)

    ## Make the following command:
    #    (exit $internal_exec_status)
    exit_ec_ast = make_exit_ec_ast()
    asts.append(exit_ec_ast)

    class_asts = [to_ast_node(ast) for ast in asts]
    return class_asts


def make_exit_ec_ast():
    command = make_command(
        [string_to_argument("exit"), [make_quoted_variable("internal_exec_status")]]
    )
    ast = make_subshell(command)
    return ast


def make_rm_f_ast(arguments):
    all_args = [string_to_argument("rm"), string_to_argument("-f")] + arguments
    return make_command(all_args)


def make_mkfifo_ast(arguments):
    all_args = [string_to_argument("mkfifo")] + arguments
    return make_command(all_args)
