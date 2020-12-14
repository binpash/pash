import os
import time
import subprocess
from datetime import datetime
from util import *
from ir_utils import *

from json_ast import save_asts_json
from parse import from_ir_to_shell
import config

##
## TODO: Make this take the real graph as input (and not the serializable one).
##       Taking the serialized graph was once done to interface with a JVM runtime.
##       However, we can make a much cleaner IR -> AST transformation as a backend here.
##

def to_shell(ir, output_dir, args):
    backend_start_time = datetime.now()

    ## First call an IR to AST compilation pass
    output_asts = ir2ast(ir, args)

    ## Then just call the parser.
    temp_filename = os.path.join("/tmp", get_random_string())
    save_asts_json(output_asts, temp_filename)
    output_script = from_ir_to_shell(temp_filename)

    ## TODO: Delete, obsolete
    # graph_json = ir.serialize_as_JSON()
    # log("JSON IR:", graph_json)
    # output_script = shell_backend(graph_json, output_dir, args)

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

    all_fids = ir.all_fids()

    ## Find all the ephemeral fids and turn them to ASTs
    ephemeral_fids = [fid for fid in all_fids
                      if fid.resource is None]

    ## Call the prologue that creates fifos for all ephemeral fids    
    prologue = make_ir_prologue(ephemeral_fids)
    
    ## Make the main body
    body = ir.to_ast(drain_streams)

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


## TODO: Delete all the obsolete code below here

def shell_backend(graph_json, output_dir, args):
    ## TODO: Remove output_dir since it is not used

    clean_up_graph = False
    drain_streams = False
    if(args.termination == "clean_up_graph"):
        clean_up_graph = True
    elif(args.termination == "drain_stream"):
        drain_streams = True
    if(args.split_fan_out > 1):
        auto_split = True
    else:
        auto_split = False
    # print("Translate:")
    # print(graph_json)
    fids = graph_json["fids"]
    in_fids = graph_json["in"]
    out_fids = graph_json["out"]
    nodes = graph_json["nodes"]

    output_script_commands = []

    ## TODO: These are osbolete since we now redirect the output outisde
    ## Make the output directory if it doesn't exist
    # output_script_commands.append('rm -rf {}'.format(output_dir))
    # output_script_commands.append('mkdir -p {}'.format(output_dir))

    ## Setup pipes
    rm_com = remove_fifos(fids)
    mkfifo_com = make_fifos(fids)
    output_script_commands.append(rm_com)
    output_script_commands.append(mkfifo_com)

    ## Redirect stdin
    ## TODO: Assume that only stdin can be an in_fid.
    assert(len(in_fids) <= 1)
    for in_fid in in_fids:
        ## TODO: Is this a hack?
        in_com = '{{ cat > "{}" <&3 3<&- & }} 3<&0'.format(in_fid)
        output_script_commands.append(in_com)

    ## Execute nodes
    processes = [execute_node(node, drain_streams, auto_split)
                 for node_id, node in nodes.items()]
    output_script_commands += processes

    ## Collect outputs
    ## TODO: Make this work for more than one output. 
    ##       For now it is fine to only have stdout as output
    ##
    ## TODO: Make this not cat if the output is a real file.
    assert(len(out_fids) == 1)
    output_com = 'cat "{}" &'.format(out_fids[0])
    output_script_commands.append(output_com)

    ## If the option to clean up the graph is enabled, we should only
    ## wait on the final pid and kill the rest using SIGPIPE.
    if (clean_up_graph):
        suffix = ""
        if (not config.pash_args.log_file == ""):
            suffix = " 2>> " + config.pash_args.log_file
        output_script_commands.append('wait $!')
        output_script_commands.append('pids_to_kill="$(ps --ppid $$ |' 
                                      "awk '{print $1}' | "
                                      "grep -E '[0-9]'" + suffix + ')"') 
        output_script_commands.append('for pid in $pids_to_kill')
        output_script_commands.append('do')
        output_script_commands.append('  (kill -SIGPIPE $pid || true)' + suffix)
        output_script_commands.append('done')
    else:
        ## Otherwise we just wait for all processes to die.
        output_script_commands.append('wait')

    ## Kill pipes
    output_script_commands.append(rm_com)

    return ("\n".join(output_script_commands) + "\n")

def remove_fifos(fids):
    ## We remove one fifo at a time because the big benchmarks crash
    ## with "too many arguments" error
    rms = []
    for fid in fids:
        assert(fid[0] == "#")
        rms.append('rm -f "{}"'.format(fid))
    return "\n".join(rms)
    # return 'rm -f {}'.format(" ".join(['"{}"'.format(fid) for fid in fids]))

def make_fifos(fids):
    ## We make one fifo at a time because the big benchmarks crash
    ## with "too many arguments" error
    mkfifos = []
    for fid in fids:
        assert(fid[0] == "#")
        mkfifos.append('mkfifo "{}"'.format(fid))
    return "\n".join(mkfifos)
    # return 'mkfifo {}'.format(" ".join(['"{}"'.format(fid) for fid in fids]))

def execute_node(node, drain_streams, auto_split):
    script = node_to_script(node, drain_streams, auto_split)
    return "{} &".format(script)

def node_to_script(node, drain_streams, auto_split):
    # print(node)

    inputs = node["in"]
    outputs = node["out"]
    command = node["command"]

    script = []
    ## Split file
    if(command.split(" ")[0] == "split_file"):
        assert(len(inputs) == 1)
        batch_size = command.split(" ")[1]
        split_outputs = command.split(" ")[2:]
        split_string = new_split(inputs[0], split_outputs, batch_size, auto_split)
        script.append(split_string)
    else:
        ## All the other nodes
        if (len(inputs) > 0):
            script.append("cat")
            for fid in inputs:
                ## Hack to not output double double quotes
                string_fid = str(fid)
                if(string_fid.startswith('"')):
                    script.append('{}'.format(string_fid))
                else:
                    script.append('"{}"'.format(string_fid))
            script.append("|")

        ## Add the command together with a drain stream
        if (drain_streams):
            script.append("(")
            script += command.split(" ")
            script.append(";")
            script += drain_stream()
            script.append(")")
        else:
            script += command.split(" ")

        if(len(outputs) == 1):
            script.append(">")
            script.append('"{}"'.format(outputs[0]))
        elif(len(outputs) == 0):
            pass
        else:
            assert(False)

        # print("Script:", script)
    return " ".join(script)

def new_split(input_file, outputs, batch_size, auto_split):
    if(auto_split):
        auto_split_bin = '{}/{}'.format(config.PASH_TOP, config.config['runtime']['auto_split_binary'])
        command_no_outputs = '{} "{}"'.format(auto_split_bin, input_file, batch_size)
    else:
        split_bin = '{}/{}'.format(config.PASH_TOP, config.config['runtime']['split_binary'])
        command_no_outputs = '{} "{}" {}'.format(split_bin, input_file, batch_size)
    return ' '.join([command_no_outputs] + outputs)

def drain_stream():
    script = '{}/{}'.format(config.PASH_TOP, config.config['distr_planner']['drain_stream_executable_path'])
    return [script]

def old_split(inputs, outputs, batch_size):
    script = []
    script.append('{')
    script.append('head -n {}'.format(batch_size))
    script.append('"{}"'.format(inputs[0]))
    script.append(">")
    script.append('"{}"'.format(outputs[0]))
    script.append(";")
    script.append("cat")
    script.append('"{}"'.format(inputs[0]))
    script.append(">")
    script.append('"{}"'.format(outputs[1]))
    script.append('; }')
    return script
