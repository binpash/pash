import os
import time
import subprocess
from datetime import datetime
from util import *
import config

def execute(graph_json, output_dir, output_script_name, output_optimized, args):
    backend_start_time = datetime.now()

    output_script = shell_backend(graph_json, output_dir)

    ## Output the optimized shell script for inspection
    if(output_optimized):
        with open(output_script_name, "w") as output_script_file:
            print("Optimized script:")
            print(output_script)
            output_script_file.write(output_script)

    backend_end_time = datetime.now()
    print_time_delta("Backend", backend_start_time, backend_end_time, args)

    if(not args.compile_optimize_only):
        execution_start_time = datetime.now()

        ## TODO: Handle stdout, stderr, errors
        exec_obj = subprocess.run(output_script, shell=True, executable="/bin/bash")
        exec_obj.check_returncode()

        execution_end_time = datetime.now()
        print_time_delta("Execution", execution_start_time, execution_end_time, args)



def shell_backend(graph_json, output_dir):
    # print("Translate:")
    # print(graph_json)
    fids = graph_json["fids"]
    in_fids = graph_json["in"]
    out_fids = graph_json["out"]
    nodes = graph_json["nodes"]

    start_time = time.time()

    output_script_commands = []

    shared_memory_dir = '/dev/shm/dish'
    ## Make the output directory if it doesn't exist
    output_script_commands.append('rm -rf {}'.format(output_dir))
    output_script_commands.append('mkdir -p {}'.format(output_dir))
    output_script_commands.append('mkdir -p {}'.format(shared_memory_dir))

    ## Setup pipes
    rm_com = remove_fifos(fids)
    mkfifo_com = make_fifos(fids)
    output_script_commands.append(rm_com)
    output_script_commands.append(mkfifo_com)

    ## Execute nodes
    processes = [execute_node(node, shared_memory_dir) for node_id, node in nodes.items()]
    output_script_commands += processes

    ## Collect outputs
    # collect_output_args = ["cat"]
    # for out_fid in out_fids:
    #     collect_output_args.append('"{}"'.format(out_fid))
    # collect_output_args.append(">")
    # collect_output_args.append('"{}"'.format(output_file))
    # output_script = " ".join(collect_output_args)
    # print("Collect output:")
    # print(output_script)
    # out_p = subprocess.Popen(output_script, shell=True,
    #                          executable="/bin/bash")
    for i, out_fid in enumerate(out_fids):
        output_com = 'cat "{}" > {}/{} &'.format(out_fid, output_dir, i)
        output_script_commands.append(output_com)

    ## Wait for all processes to die
    # for proc in processes:
    #     ret = proc.wait()
    #     if(not ret == 0):
    #         print("-- Error!", proc, ret)
    # output_script_commands.append('for job in `jobs -p` \ndo \n echo $job\n wait $job \ndone')
    output_script_commands.append('wait')

    ## Kill pipes
    final_rm_com = remove_fifos(fids)
    output_script_commands.append(rm_com)

    output_script_commands.append('rm -rf "{}"'.format(shared_memory_dir))
    end_time = time.time()

    ## TODO: Cat all outputs together if a specific flag is given

    # print("Distributed translation execution time:", end_time - start_time)
    return "\n".join(output_script_commands)

def remove_fifos(fids):
    for fid in fids:
        assert(fid[0] == "#")
    return 'rm -f {}'.format(" ".join(['"{}"'.format(fid) for fid in fids]))

def make_fifos(fids):
    for fid in fids:
        assert(fid[0] == "#")
    return 'mkfifo {}'.format(" ".join(['"{}"'.format(fid) for fid in fids]))

def execute_node(node, shared_memory_dir):
    script = node_to_script(node, shared_memory_dir)
    return "{} &".format(script)

def node_to_script(node, shared_memory_dir):
    # print(node)

    inputs = node["in"]
    outputs = node["out"]
    command = node["command"]

    script = []
    ## Split file
    if(len(outputs) == 2 and
       command.split(" ")[0] == "split_file"):
        assert(len(inputs) == 1)
        batch_size = command.split(" ")[1]
        if (len(inputs) > 0):
            script.append("cat")
            for fid in inputs:
                script.append('"{}"'.format(fid))
            script.append("|")
        script += new_split(outputs, batch_size, shared_memory_dir)
        # script += old_split(inputs, outputs, batch_size)
        # print(script)
        # print(node)
    else:
        ## All the other nodes
        if (len(inputs) > 0):
            script.append("cat")
            for fid in inputs:
                script.append('"{}"'.format(fid))
            script.append("|")

        ## Add the command together with a drain stream
        script.append("(")
        script += command.split(" ")
        script.append(";")
        script += drain_stream()
        script.append(")")

        if(len(outputs) == 1):
            script.append(">")
            script.append('"{}"'.format(outputs[0]))
        elif(len(outputs) == 0):
            pass
        else:
            assert(False)

        # print("Script:", script)
    return " ".join(script)

def new_split(outputs, batch_size, shared_memory_dir):
    script = []
    shared_mem_file = '"{}/{}"'.format(shared_memory_dir, outputs[0])
    script.append('tee >(')
    script.append('head -n {}'.format(batch_size))
    script.append('> {};'.format(shared_mem_file))
    # script.append('dd of="{}" > /dev/null 2>&1'.format(outputs[1]))
    script.append('dd of=/dev/null > /dev/null 2>&1 & cat {} > "{}")'.format(shared_mem_file, outputs[0]))
    script.append('| (tail -n +{}'.format(int(batch_size)+1))
    script.append('> "{}";'.format(outputs[1]))
    script.append('dd of=/dev/null > /dev/null 2>&1)')
    # script.append('"{}"'.format(outputs[1]))
    # script.append('; }')
    return script

def drain_stream():
    script = '{}/{}'.format(config.DISH_TOP, config.config["drain_stream_executable_path"])
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
