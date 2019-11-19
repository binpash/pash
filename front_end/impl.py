import os
import time
import subprocess

def execute(graph_json, output_file):
    ## TODO: Remove environment hardcoding
    env={"IN": "../scripts/input/i1M.txt"}

    print("Execute:")
    # print(graph_json)
    fids = graph_json["fids"]
    in_fids = graph_json["in"]
    out_fids = graph_json["out"]
    nodes = graph_json["nodes"]


    start_time = time.time()

    ## Setup pipes
    for fid in fids:
        # print(fid)
        remove_fifo_if_exists(fid)
        make_fifo(fid)

    ## Execute nodes
    processes = [execute_node(node, env) for node_id, node in nodes.items()]

    ## Collect outputs
    collect_output_args = ["cat"]
    for out_fid in out_fids:
        collect_output_args.append('"{}"'.format(out_fid))
    collect_output_args.append(">")
    collect_output_args.append('"{}"'.format(output_file))
    output_script = " ".join(collect_output_args)
    print("Collect output:")
    print(output_script)
    out_p = subprocess.Popen(output_script, shell=True,
                             executable="/bin/bash")

    ## Wait for all processes to die
    for proc in processes:
        ret = proc.wait()
        if(not ret == 0):
            print("Error!", proc, ret)

    ## Kill pipes
    for fid in fids:
        # print(fid)
        remove_fifo_if_exists(fid)

    end_time = time.time()

    print("Distributed graph execution time:", end_time - start_time)

def remove_fifo_if_exists(pipe):
    assert(pipe[0] == "#")
    if os.path.exists(pipe):
        os.remove(pipe)

def make_fifo(pipe):
    assert(pipe[0] == "#")
    os.mkfifo(pipe)

def execute_node(node, env):
    # print(node)
    script = node_to_script(node)
    print(script)
    p = subprocess.Popen(script, shell=True,
                         executable="/bin/bash", env=env)
    return p

def node_to_script(node):
    # print(node)

    inputs = node["in"]
    outputs = node["out"]
    command = node["command"]

    script = []
    if (len(inputs) > 0):
        script.append("cat")
        for fid in inputs:
            script.append('"{}"'.format(fid))
        script.append("|")

    script += command.split(" ")

    if(len(outputs) == 1):
        script.append(">")
        script.append('"{}"'.format(outputs[0]))
    else:
    ## TODO: Implement split
        assert(False)

    # print("Script:", script)
    return " ".join(script)
