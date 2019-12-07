import os
import time
import subprocess

def execute(graph_json, output_dir):
    ## TODO: Remove environment hardcoding
    env={"IN": "../scripts/input/i1M.txt"}

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
    for fid in fids:
        # print(fid)
        rm_com = remove_fifo_if_exists(fid)
        mkfifo_com = make_fifo(fid)
        output_script_commands.append(rm_com)
        output_script_commands.append(mkfifo_com)

    ## Execute nodes
    processes = [execute_node(node, env, shared_memory_dir) for node_id, node in nodes.items()]
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
    for fid in fids:
        # print(fid)
        rm_com = remove_fifo_if_exists(fid)
        output_script_commands.append(rm_com)

    output_script_commands.append('rm -rf "{}"'.format(shared_memory_dir))
    end_time = time.time()

    ## TODO: Cat all outputs together if a specific flag is given

    # print("Distributed translation execution time:", end_time - start_time)
    return "\n".join(output_script_commands)

def remove_fifo_if_exists(pipe):
    assert(pipe[0] == "#")
    # print('rm -f "{}"'.format(pipe))
    # if os.path.exists(pipe):
    #     os.remove(pipe)
    return 'rm -f "{}"'.format(pipe)

def make_fifo(pipe):
    assert(pipe[0] == "#")
    # print ('mkfifo "{}"'.format(pipe))
    # os.mkfifo(pipe)
    return 'mkfifo "{}"'.format(pipe)

def execute_node(node, env, shared_memory_dir):
    # print(node)
    script = node_to_script(node, shared_memory_dir)
    # print("{} &".format(script))
    # p = subprocess.Popen(script, shell=True,
    #                      executable="/bin/bash", env=env)
    # return p
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

        script += command.split(" ")

        if(len(outputs) == 1):
            script.append(">")
            script.append('"{}"'.format(outputs[0]))
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
