import os
import time
import subprocess
from datetime import datetime
from util import *
from ir_utils import *

from json_ast import save_asts_json
from parse import from_ir_to_shell
import config

def to_shell(ir, output_dir, args):
    backend_start_time = datetime.now()

    ## First call an IR to AST compilation pass
    output_asts = ir2ast(ir, args)

    ## Then just call the parser.
    temp_filename = os.path.join("/tmp", get_random_string())
    save_asts_json(output_asts, temp_filename)
    output_script = from_ir_to_shell(temp_filename)

    backend_end_time = datetime.now()
    print_time_delta("Backend", backend_start_time, backend_end_time, args)

    return output_script


def ir2ast(ir, args):
    clean_up_graph = False
    drain_streams = False
    if(args.termination == "clean_up_graph"):
        clean_up_graph = True
    elif(args.termination == "drain_stream"):
        drain_streams = True

    ## NOTE: We first need to make the main body because it might create additional ephemeral fids.

    ## Make the main body
    body = ir.to_ast(drain_streams)

    all_fids = ir.all_fids()

    log("All fids:", all_fids)
    ## Find all the ephemeral fids and turn them to ASTs
    ephemeral_fids = [fid for fid in all_fids
                      if fid.is_ephemeral()]

    log("Ephemeral fids:", ephemeral_fids)

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

def make_ir_prologue(ephemeral_fids):
    asts = []
    ## Create an `rm -f` for each ephemeral fid
    asts += make_rms_f_prologue_epilogue(ephemeral_fids)

    ## Create a `mkfifo` for each ephemeral fid
    for eph_fid in ephemeral_fids:
        args = [eph_fid.to_ast()]
        command = make_mkfifo_ast(args)
        asts.append(command)
    
    return asts

def make_ir_epilogue(ephemeral_fids, clean_up_graph, log_file):
    asts = []
    if (clean_up_graph):
        ## TODO: Wait for all output nodes not just one
        pids = [[standard_var_ast('!')]]
        clean_up_path_script = os.path.join(config.PASH_TOP, config.config['runtime']['clean_up_graph_binary'])
        com_args = [string_to_argument('source'), string_to_argument(clean_up_path_script)] + pids
        if (log_file == ""):
            com = make_command(com_args)
        else:
            redirection = redir_append_stderr_to_string_file(log_file)
            com = make_command(com_args, redirections=[redirection])
        asts.append(com)
    else:
        ## Otherwise we just wait for all processes to die.
        wait_com = make_command(string_to_argument('wait'))
        asts.append(wait_com)

    ## Create an `rm -f` for each ephemeral fid
    asts += make_rms_f_prologue_epilogue(ephemeral_fids)
    return asts

def make_rm_f_ast(arguments):
    all_args = [string_to_argument("rm"), string_to_argument("-f")] + arguments
    return make_command(all_args)

def make_mkfifo_ast(arguments):
    all_args = [string_to_argument("mkfifo")] + arguments
    return make_command(all_args)

## Keeping the old code below just in case we want to reimplement it

# def new_split(input_file, outputs, batch_size, auto_split):
#     if(auto_split):
#         auto_split_bin = '{}/{}'.format(config.PASH_TOP, config.config['runtime']['auto_split_binary'])
#         command_no_outputs = '{} "{}"'.format(auto_split_bin, input_file, batch_size)
#     else:
#         split_bin = '{}/{}'.format(config.PASH_TOP, config.config['runtime']['split_binary'])
#         command_no_outputs = '{} "{}" {}'.format(split_bin, input_file, batch_size)
#     return ' '.join([command_no_outputs] + outputs)

# def drain_stream():
#     script = '{}/{}'.format(config.PASH_TOP, config.config['distr_planner']['drain_stream_executable_path'])
#     return [script]
