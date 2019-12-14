import sys
import pickle
import subprocess
import jsonpickle
import yaml

from ir import *
from json_ast import *
from impl import execute

## This file receives the name of a file that holds an IR, reads the
## IR, read some configuration file with node information, and then
## should make a distribution plan for it.

config = {}

def load_config():
    global config
    dish_config = {}
    CONFIG_KEY = 'distr_planner'

    ## TODO allow this to be passed as an argument
    config_file_path = 'config.yaml'
    with open(config_file_path) as config_file:
        dish_config = yaml.load(config_file, Loader=yaml.FullLoader)

    if not dish_config:
        raise Exception('No valid configuration could be loaded from {}'.format(config_file_path))

    if CONFIG_KEY not in dish_config:
        raise Exception('Missing `{}` config in {}'.format(CONFIG_KEY, config_file_path))

    config = dish_config[CONFIG_KEY]

def main():
    load_config()
    optimize_script()

def optimize_script():
    global config
    if not config:
        load_config()

    with open(config['ir_filename'], "rb") as ir_file:
        ir_node = pickle.load(ir_file)

    print("Retrieving IR: {} ...".format(config['ir_filename']))
    shell_string = ast_to_shell(ir_node.ast)
    print(shell_string)

    distributed_graph = simpl_file_distribution_planner(ir_node, config['fan_out'], config['batch_size'])

    ## Notes:
    ##
    ## Each node in the graph has one input stream and one output
    ## stream.
    ##
    ## However, the input (and the output) might be in several
    ## different "locations". So the input (and output) can be
    ## represented as a sequence of resources. These can be files,
    ## urls, blocks of files in HDFS, etc.
    ##
    ## The distribution planner takes a graph, where the edges
    ## arriving to a node are uniquely numbered depending on the
    ## position that they have in the sequence of inputs to the
    ## incoming node.
    ##
    ## The distribution planner also takes information about the
    ## available computational resources, such as nodes in the system,
    ## and their cores, etc.

    ## Call the backend that executes the optimized dataflow graph
    execute(distributed_graph.serialize_as_JSON(), config['output_dir'], config['output_script_name'], config['output_optimized'])

## This is a simplistic planner, that pushes the available
## parallelization from the inputs in file stateless commands. The
## planner starts from the sources of the graph, and pushes
## file parallelization as far as possible.
##
## It returns a maximally expanded (regarding files) graph, that can
## be scheduled depending on the available computational resources.
def simpl_file_distribution_planner(graph, fan_out, batch_size):
    # print("Intermediate representation:")
    # print(graph)

    ## We assume that the file identifiers that have been added in the
    ## intermediate representation show the edges between different IR
    ## nodes.
    ##
    ## We assume that all nodes have an in_stream and an out_stream
    ## list, and that these are the ones which will be used to create
    ## the graph.

    ## Starting from the sources of the graph, if they are stateless,
    ## duplicate the command as many times as the number of
    ## identifiers in its in_stream. Then connect their outputs in
    ## order to next command.
    ##
    ## Note: The above has to be done in a BFS fashion (starting from
    ## all sources simultaneously) so that we don't have to iterate to
    ## reach a fixpoint.
    naive_parallelize_stateless_nodes_bfs(graph, fan_out, batch_size)
    # print("Parallelized graph:")
    # print(graph)

    # Output the graph as json
    # frozen = jsonpickle.encode(graph)
    # f = open("minimal2_ir.json", "w")
    # f.write(frozen)
    # f.close()

    # print(graph.serialize())
    # f = open("serialized_ir", "w")
    # f.write(graph.serialize_as_JSON_string())
    # f.close()

    ## The result of the above steps should be an expanded
    ## intermediate representation graph, that can be then mapped to
    ## real nodes.
    return graph


def naive_parallelize_stateless_nodes_bfs(graph, fan_out, batch_size):
    source_nodes = graph.source_nodes()
    # print("Source nodes:")
    # print(source_nodes)

    ## Generate a fileIdGen from a graph, that doesn't clash with the
    ## current graph fileIds.
    fileIdGen = graph.get_file_id_gen()

    ## If the source nodes only have one file input, then split it in
    ## partial files.


    # commands_to_split_input = ["cat"]
    # for source_node in source_nodes:
    #     input_file_ids = source_node.get_flat_input_file_ids()
    #     ## TODO: Also split when we have more than one input file
    #     if(len(input_file_ids) == 1 and
    #        str(source_node.command) in commands_to_split_input):
    #         input_file_id = input_file_ids[0]
    #         input_file_id.split_resource(2, batch_size, fileIdGen)

    nodes = source_nodes
    while (len(nodes) > 0):
        curr = nodes.pop(0)

        next_nodes = graph.get_next_nodes(curr)
        nodes += next_nodes

        ## Question: What does it mean for a command to have more
        ## than one next_node? Does it mean that it duplicates its
        ## output to all of them? Or does it mean that it writes
        ## some to the first and some to the second? Both are not
        ## very symmetric, but I think I would prefer the first.
        assert(len(next_nodes) <= 1)

        ## TODO: Remove hardcoded
        graph = split_command_input(curr, graph, fileIdGen, fan_out, batch_size)
        graph = parallelize_command(curr, graph, fileIdGen)

## Optimizes several commands by splitting its input
def split_command_input(curr, graph, fileIdGen, fan_out, batch_size):
    input_file_ids = curr.get_flat_input_file_ids()
    if (curr.category in ["stateless", "pure"] and
        len(input_file_ids) == 1 and
        curr.in_stream[0] == "stdin" and
        fan_out > 1):
        ## We can split command input in several files, as long as the
        ## input is not a file, and it makes sense to do so, if the
        ## command is stateless or pure.

        ## This will still be the output file id of the previous
        ## node
        input_file_id = input_file_ids[0]

        ## Generate the split file ids
        split_file_ids = [fileIdGen.next_file_id() for i in range(fan_out)]

        ## Generate a new file id that has all the split file ids
        ## as children, in place of the current one
        output_file_id = fileIdGen.next_file_id()
        output_file_id.set_children(split_file_ids)

        ## Set this new file id to be the input ot xargs.
        curr.stdin = output_file_id

        ## Add a new node that executes split_file and takes
        ## the input_file_id as input, split_file_ids as output
        ## and the batch_size as an argument
        split_file_commands = make_split_files(input_file_id, output_file_id, batch_size, fileIdGen)
        [graph.add_node(split_file_command) for split_file_command in split_file_commands]

    return graph

def parallelize_command(curr, graph, fileIdGen):
    ## If the command is stateless or pure parallelizable and has more
    ## than one inputs, it can be parallelized.
    input_file_ids = curr.get_flat_input_file_ids()
    if ((curr.category == "stateless" or curr.is_pure_parallelizable())
        and len(input_file_ids) > 1):

        ## We assume that every command has one output_file_id for
        ## now. I am not sure if this must be lifted. This seems
        ## to be connected to assertion regarding the next_nodes.
        out_edge_file_ids = curr.get_output_file_ids()
        assert(len(out_edge_file_ids) == 1)
        out_edge_file_id = out_edge_file_ids[0]

        ## Add children to the output file (thus also changing the
        ## input file of the next command to have children)
        new_output_file_ids = [fileIdGen.next_file_id() for in_fid in input_file_ids]

        if(curr.category == "stateless"):
            out_edge_file_id.set_children(new_output_file_ids)

            ## For each new input and output file id, make a new command
            new_commands = curr.stateless_duplicate()
        elif(curr.is_pure_parallelizable()):
            intermediate_output_file_id = fileIdGen.next_file_id()
            intermediate_output_file_id.set_children(new_output_file_ids)
            curr.stdout = intermediate_output_file_id

            ## TODO: Generalize this part to not only work for sort
            ## For each new input and output file id, make a new command
            new_commands = sort_duplicate(curr, fileIdGen)
            merge_command = create_sort_merge_command(curr, new_output_file_ids, out_edge_file_id)
            graph.add_node(merge_command)
        else:
            print("Unreachable code reached :(")
            assert(False)
            ## This should be unreachable

        ## WARNING: In order for the above to not mess up
        ## anything, there must be no other node that writes to
        ## the same output as the curr node. Otherwise, the above
        ## procedure will mess this up.
        ##
        ## TODO: Either make an assertion to catch any case that
        ## doesn't satisfy the above assumption here, or extend
        ## the intermediate representation and the above procedure
        ## so that this assumption is lifted (either by not
        ## parallelizing, or by properly handling this case)

        graph.remove_node(curr)
        # print("New commands:")
        # print(new_commands)
        for new_com in new_commands:
            graph.add_node(new_com)

        return graph


## TODO: Find a better place to have this and the function below.
def sort_duplicate(curr, fileIdGen):
    assert(str(curr.command) == "sort")
    input_file_ids = curr.get_flat_input_file_ids()
    output_file_ids = curr.get_flat_output_file_ids()

    in_out_file_ids = zip(input_file_ids, output_file_ids)

    new_commands = [curr.make_duplicate_command(in_fid, out_fid) for in_fid, out_fid in in_out_file_ids]

    return new_commands

## TODO: This must be generated using some file information
def create_sort_merge_command(curr, new_output_file_ids, out_edge_file_id):
    ## Create a merge sort command
    ##
    ## WARNING: Since the merge has to take the files as arguments, we
    ## pass the pipe names as its arguments and nothing in stdin
    old_options = curr.get_non_file_options()
    input_file_ids = curr.get_flat_input_file_ids()
    ## TODO: Implement a proper version of parallel sort -m, instead
    ## of using the parallel flag.
    options = [string_to_argument("-m"), string_to_argument("--parallel={}".format(len(input_file_ids)))] + old_options
    # options = []
    options += [string_to_argument(fid.pipe_name()) for fid in new_output_file_ids]
    opt_indices = [("option", i) for i in range(len(options))]
    # in_stream = [("option", i + len(opt_indices)) for i in range(len(new_output_file_ids))]
    in_stream = []
    merge_command = Command(None, # TODO: Make a proper AST
                            curr.command,
                            options,
                            in_stream,
                            ["stdout"],
                            opt_indices,
                            "pure",
                            # intermediate_output_file_id,
                            [],
                            out_edge_file_id)
    return merge_command

def parallelize_pure_wc(curr, graph, fileIdGen):
    input_file_ids = curr.get_flat_input_file_ids()
    assert(len(input_file_ids) > 1)

    ## We assume that every command has one output_file_id for
    ## now. I am not sure if this must be lifted. This seems
    ## to be connected to assertion regarding the next_nodes.
    out_edge_file_ids = curr.get_output_file_ids()
    assert(len(out_edge_file_ids) == 1)
    out_edge_file_id = out_edge_file_ids[0]

    ## Add children to the output file (thus also changing the
    ## input file of the next command to have children)
    new_output_file_ids = [fileIdGen.next_file_id() for in_fid in input_file_ids]
    intermediate_output_file_id = fileIdGen.next_file_id()
    intermediate_output_file_id.set_children(new_output_file_ids)

    curr.stdout = intermediate_output_file_id
    ## For each new input and output file id, make a new command
    new_commands = sort_duplicate(curr, fileIdGen)
    # print("New commands:")
    # print(new_commands)

    ## Create a merge sort command
    ##
    ## WARNING: Since the merge has to take the files as arguments, we
    ## pass the pipe names as its arguments and nothing in stdin
    old_options = curr.get_non_file_options()
    ## TODO: Implement a proper version of parallel sort -m, instead
    ## of using the parallel flag.
    options = [string_to_argument("-m"), string_to_argument("--parallel={}".format(len(input_file_ids)))]
    # options = []
    options += [string_to_argument(fid.pipe_name()) for fid in new_output_file_ids]
    opt_indices = [("option", i) for i in range(len(options))]
    # in_stream = [("option", i + len(opt_indices)) for i in range(len(new_output_file_ids))]
    in_stream = []
    merge_command = Command(None, # TODO: Make a proper AST
                            curr.command,
                            options,
                            in_stream,
                            ["stdout"],
                            opt_indices,
                            "pure",
                            # intermediate_output_file_id,
                            [],
                            out_edge_file_id)

    graph.remove_node(curr)
    for new_com in new_commands:
        graph.add_node(new_com)

    graph.add_node(merge_command)

    return graph



    ## TODO: In order to be able to execute it, we either have to
    ## execute it in the starting shell (so that we have its state),
    ## or we should somehow pass the parent shell's state to the the
    ## distribution planner, and then the implementation environment.
    ## In general, we probably have to find a way to pass around a
    ## shell's state, as this will be essential for the distributed
    ## setting too.
    ##
    ## Note: A way to do this is by using set > temp_file. Source:
    ## https://arstechnica.com/civis/viewtopic.php?f=16&t=805521

    ## TODO: We have to handle xargs in a special way. First of all,
    ## in order to parallelize the command that xarg runs, we have to
    ## do xargs -L 1 (or some other number) so that for every -L
    ## lines, it calls a different instance of the command. Then it
    ## will be parallelizable. In addition, we have to somehow
    ## statically decide how much we will parallelize xarg, and how
    ## many lines are going to be sent to each operator.
    ##
    ## This can probably be solved if we allow partial files in files
    ## without resources.

    ## TODO: There is slight problem with *, and other expansions in
    ## the shell. The normal shell semantics is to expand the strings
    ## in a command that is in a pipeline after the different
    ## subshells have been spawned. However, we would like to have all
    ## strings expanded as much as possible, so that we can statically
    ## make choices about how much to distribute each command.
    ##
    ## Maybe we should run expansions on our own, before calling the
    ## distribution planner? Or in the distribution planner itself? It
    ## seems that the distribution planner should be able to do some
    ## expansion itself though

    ## TODO: There is a problem when given an unexpanded string. It
    ## might be many files, so spliting the file up in different
    ## pieces might be wrong.

    ## BIG TODO: Extend the file class so that it supports tee etc.

if __name__ == "__main__":
    main()
