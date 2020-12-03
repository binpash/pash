import os
import time
import subprocess
from datetime import datetime
from util import *
import config

def to_shell(graph_json, output_dir, args):
    backend_start_time = datetime.now()

    output_script = shell_backend(graph_json, output_dir, args)

    backend_end_time = datetime.now()
    print_time_delta("Backend", backend_start_time, backend_end_time, args)

    return output_script



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

    start_time = time.time()

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

    ## Execute nodes
    processes = [execute_node(node, drain_streams, auto_split)
                 for node_id, node in nodes.items()]
    output_script_commands += processes

    ## Collect outputs
    ## TODO: Make this work for more than one output. 
    ##       For now it is fine to only have stdout as output
    assert(len(out_fids) == 1)
    output_com = 'cat "{}" &'.format(out_fids[0])
    output_script_commands.append(output_com)
    # ## Old moving of outputs to a temporary directory
    # for i, out_fid in enumerate(out_fids):
    #     output_com = 'cat "{}" > {}/{} &'.format(out_fid, output_dir, i)
    #     output_script_commands.append(output_com)

    ## If the option to clean up the graph is enabled, we should only
    ## wait on the final pid and kill the rest using SIGPIPE.
    if (clean_up_graph):
        output_script_commands.append('wait $!')
        output_script_commands.append("ps --ppid $$ | awk '{print $1}'"
                                      " | grep -E '[0-9]' | xargs -n 1 kill -SIGPIPE")
    else:
        ## Otherwise we just wait for all processes to die.
        output_script_commands.append('wait')

    ## Kill pipes
    final_rm_com = remove_fifos(fids)
    output_script_commands.append(rm_com)

    end_time = time.time()

    ## TODO: Cat all outputs together if a specific flag is given

    # print("Distributed translation execution time:", end_time - start_time)
    return "\n".join(output_script_commands)

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
