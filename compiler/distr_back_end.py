import os
import time
import subprocess
from ir import *

def distr_execute(graph, output_dir, output_script_name, output_optimized, compile_optimize_only, distributed_nodes):

    output_scripts = shell_backend(graph, output_dir, distributed_nodes)

    ## TODO: Do something with the output_scripts

    # ## Output the optimized shell script for inspection
    # if(output_optimized):
    #     with open(output_script_name, "w") as output_script_file:
    #         output_script_file.write(output_script)

    # if(not compile_optimize_only):
    #     ## TODO: Handle stdout, stderr, errors
    #     exec_obj = subprocess.run(output_script, shell=True, executable="/bin/bash")
    #     exec_obj.check_returncode()


def map_graph_to_physical_nodes(graph, number_nodes):

    ## This keeps a list of file ids that will be connected to remote
    ## channels in each node. It keeps triples of file id, previous,
    ## and next physical node.
    remote_channels = []

    ## This will be a list of graphs, one for each physical node. Graphs are in
    ## the exact same format as the initial graph.
    distributed_graphs = []
    curr_graph_nodes = []
    curr_graph_inputs = []
    curr_graph_outputs = []

    ## The DFS returns distributed graphs.
    ## 1. Given the input fids we can find the source nodes of the graph
    ## 2. We keep a working stack of the node identifiers and every new node
    ##    that we visit, we add it to the current physical node.
    ## 3. Whenever we encounter a node of the graph that is already assigned to 
    ##    a physical node we stop, and complete the graph for the current node.
    curr_physical_node = 0
    work_stack = graph.source_nodes()

    ## WARNING: At the moment the work stack is not exactly correct,
    ## as it contains reduce nodes too. This won't be a correctness
    ## problem though, as all of the nodes will be assigned to some
    ## physical one.
    # print(work_stack)
    while (len(work_stack) > 0):
        curr = work_stack.pop(0)
        next_nodes = graph.get_next_nodes(curr)
        # print("Next:", curr, next_nodes)
        work_stack = next_nodes + work_stack

        ## If the node is a sink of the graph, then we have to stop,
        ## so we fix its physical node to trigger the creation of the
        ## distributed graph.
        if(len(next_nodes) == 0 and not hasattr(curr, 'physical_node')):
            curr.physical_node = curr_physical_node
            curr_graph_nodes.append(curr)
            ## Since this is the last node, its output is also the
            ## output of the current graph
            curr_graph_outputs += curr.get_flat_output_file_ids()
            if(len(curr_graph_inputs) == 0):
                curr_graph_inputs += curr.get_flat_input_file_ids()

        if(hasattr(curr, 'physical_node')):
            ## Stop and close the distributed graph for
            ## curr_physical_node. The final nodes are all being put
            ## into the central (last) node.
            if(curr_physical_node + 1 < number_nodes and len(curr_graph_nodes) > 0):
                curr_physical_node += 1
                curr_graph = IR(curr_graph_nodes)
                distributed_graphs.append((curr_graph, curr_graph_inputs, curr_graph_outputs))
                curr_graph_nodes = []
                curr_graph_inputs = []
                curr_graph_outputs = []


            input_nodes_and_edges = graph.get_previous_nodes_and_edges(curr)
            new_remote_channels = [(edge, node.physical_node, curr.physical_node)
                                   for node, edge in input_nodes_and_edges
                                   if hasattr(node, 'physical_node') and
                                   not node.physical_node == curr.physical_node]
            remote_channels += new_remote_channels
        else:
            curr.physical_node = curr_physical_node
            curr_graph_nodes.append(curr)
            if(len(curr_graph_inputs) == 0):
                curr_graph_inputs += curr.get_flat_input_file_ids()

        # print("Physical:", curr, curr.physical_node)

    ## Remove duplicates
    remote_channels = list(set(remote_channels))

    return distributed_graphs, remote_channels

def json_graph_to_sh(graph_json, graph_inputs, graph_outputs, output_dir, remote_channels, distributed_nodes):
    fids = graph_json["fids"]

    print("Graph inputs:", graph_inputs)
    print("Graph outputs:", graph_outputs)

    ## TODO: Replcae those with graph_inputs, graph_outputs
    in_fids = graph_json["in"]
    out_fids = graph_json["out"]
    nodes = graph_json["nodes"]

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
    processes = [execute_node(node, shared_memory_dir) for node_id, node in nodes.items()]
    output_script_commands += processes

    ## TODO: Instead of collecting outputs with cat, we need to send
    ## them back to the central node
    for i, out_fid in enumerate(out_fids):
        output_com = 'cat "{}" > {}/{} &'.format(out_fid, output_dir, i)
        output_script_commands.append(output_com)

    ## Wait for all processes to die
    output_script_commands.append('wait')
    ## Kill pipes
    for fid in fids:
        # print(fid)
        rm_com = remove_fifo_if_exists(fid)
        output_script_commands.append(rm_com)

    output_script_commands.append('rm -rf "{}"'.format(shared_memory_dir))

    return "\n".join(output_script_commands)

def shell_backend(graph, output_dir, distributed_nodes):
    # print("Translate:")
    # print(graph)

    start_time = time.time()

    number_nodes = len(distributed_nodes)
    distributed_graphs, remote_channels = map_graph_to_physical_nodes(graph, number_nodes)
    distributed_graphs_json = [(dgraph.serialize_as_JSON(), g_ins, g_outs) for dgraph, g_ins, g_outs in distributed_graphs]
    print(distributed_graphs_json)
    print(remote_channels)
    sh_scripts = [json_graph_to_sh(graph_json, g_ins, g_outs, output_dir, remote_channels, distributed_nodes)
                  for graph_json, g_ins, g_outs in distributed_graphs_json]
    for i, sh_script in enumerate(sh_scripts):
        print("Node:", i, ", ip:", distributed_nodes[i])
        print(sh_script)

    end_time = time.time()


    # print("Distributed translation execution time:", end_time - start_time)
    return sh_scripts

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

def execute_node(node, shared_memory_dir):
    # print(node)
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
