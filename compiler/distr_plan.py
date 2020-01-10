import os
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

DISH_TOP_VAR = "DISH_TOP"

config = {}

def load_config():
    global config
    dish_config = {}
    CONFIG_KEY = 'distr_planner'

    ## TODO: allow this to be passed as an argument
    config_file_path = '{}/compiler/config.yaml'.format(os.environ[DISH_TOP_VAR])
    with open(config_file_path) as config_file:
        dish_config = yaml.load(config_file, Loader=yaml.FullLoader)

    if not dish_config:
        raise Exception('No valid configuration could be loaded from {}'.format(config_file_path))

    if CONFIG_KEY not in dish_config:
        raise Exception('Missing `{}` config in {}'.format(CONFIG_KEY, config_file_path))

    config = dish_config[CONFIG_KEY]

def optimize_script(output_script_path, compile_optimize_only):
    global config
    if not config:
        load_config()

    with open(config['ir_filename'], "rb") as ir_file:
        ir_node = pickle.load(ir_file)

    print("Retrieving IR: {} ...".format(config['ir_filename']))
    shell_string = ast_to_shell(ir_node.ast)
    print(shell_string)

    print(ir_node)
    distributed_graph = naive_parallelize_stateless_nodes_bfs(ir_node, config['fan_out'], config['batch_size'])
    print(distributed_graph)
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

    ## Call the backend that executes the optimized dataflow graph
    execute(distributed_graph.serialize_as_JSON(), config['output_dir'], output_script_path, config['output_optimized'], compile_optimize_only)

## This is a simplistic planner, that pushes the available
## parallelization from the inputs in file stateless commands. The
## planner starts from the sources of the graph, and pushes
## file parallelization as far as possible.
##
## It returns a maximally expanded (regarding files) graph, that can
## be scheduled depending on the available computational resources.

## We assume that the file identifiers that have been added in the
## intermediate representation show the edges between different IR
## nodes.
##
## We assume that all nodes have an in_stream and an out_stream
## list, and that these are the ones which will be used to create
## the graph.
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

    ## Starting from the sources of the graph, if they are stateless,
    ## duplicate the command as many times as the number of
    ## identifiers in its in_stream. Then connect their outputs in
    ## order to next command.
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
    return graph

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
        elif(curr.is_pure_parallelizable()):
            ## If the command is pure, there is need for a merging
            ## command to be added at the end.
            intermediate_output_file_id = fileIdGen.next_file_id()
            intermediate_output_file_id.set_children(new_output_file_ids)
            curr.stdout = intermediate_output_file_id

            ## Make a merge command that joins the results of all the
            ## duplicated commands
            merge_commands = create_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen)
            for merge_command in merge_commands:
                graph.add_node(merge_command)
        else:
            print("Unreachable code reached :(")
            assert(False)
            ## This should be unreachable

        ## For each new input and output file id, make a new command
        new_commands = curr.stateless_duplicate()
        graph.remove_node(curr)
        # print("New commands:")
        # print(new_commands)
        for new_com in new_commands:
            graph.add_node(new_com)

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

    return graph

## TODO: Extend the merge creation process to work with a quasi-binary
## tree, instead of just one step.

## TODO: Make the sort and bigram_aux map/reduce commands as subtypes
## of command

## Creates a merge command for all pure commands that can be
## parallelized using a map and a reduce/merge step
def create_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen):
    if(str(curr.command) == "sort"):
        return create_sort_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen)
    else:
        assert(False)

## TODO: This must be generated using some file information
##
## TODO: Find a better place to put this function
def create_sort_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen):

    tree = create_reduce_tree(SortGReduce, [[fid] for fid in new_output_file_ids], out_edge_file_id, fileIdGen)
    print(tree)
    # ## Create a merge sort command
    # ##
    # ## WARNING: Since the merge has to take the files as arguments, we
    # ## pass the pipe names as its arguments and nothing in stdin
    # old_options = curr.get_non_file_options()
    # input_file_ids = curr.get_flat_input_file_ids()
    # ## TODO: Implement a proper version of parallel sort -m, instead
    # ## of using the parallel flag.
    # options = [string_to_argument("-m"), string_to_argument("--parallel={}".format(len(input_file_ids)))] + old_options
    # # options = []
    # options += [string_to_argument(fid.pipe_name()) for fid in new_output_file_ids]
    # opt_indices = [("option", i) for i in range(len(options))]
    # # in_stream = [("option", i + len(opt_indices)) for i in range(len(new_output_file_ids))]
    # in_stream = []
    # merge_command = Command(None, # TODO: Make a proper AST
    #                         curr.command,
    #                         options,
    #                         in_stream,
    #                         ["stdout"],
    #                         opt_indices,
    #                         "pure",
    #                         # intermediate_output_file_id,
    #                         [],
    #                         out_edge_file_id)
    # return [merge_command]
    return tree

## This function creates the reduce tree. Both input and output file
## ids must be lists of lists, as the input file ids and the output
## file ids might contain auxiliary files.
def create_reduce_tree(init_func, input_file_ids, output_file_id, fileIdGen):
    tree = []
    curr_file_ids = input_file_ids
    while(len(curr_file_ids) > 1):
        new_level, curr_file_ids = create_reduce_tree_level(init_func, curr_file_ids, fileIdGen)
        tree += new_level
    ## Get the main output file identifier from the reduce and union
    ## it with the wanted output file identifier
    curr_file_ids[0][0].union(output_file_id)
    return tree


## This function creates a level of the reduce tree. Both input and
## output file ids must be lists of lists, as the input file ids and
## the output file ids might contain auxiliary files.
def create_reduce_tree_level(init_func, input_file_ids, fileIdGen):
    if(len(input_file_ids) % 2 == 0):
        output_file_ids = []
        even_input_file_ids = input_file_ids
    else:
        output_file_ids = [input_file_ids[0]]
        even_input_file_ids = input_file_ids[1:]

    level = []
    for i in range(0, len(even_input_file_ids), 2):
        new_out_file_ids = [fileIdGen.next_file_id() for _ in input_file_ids[i]]
        output_file_ids.append(new_out_file_ids)
        new_node = create_reduce_node(init_func, even_input_file_ids[i:i+2], new_out_file_ids)
        level.append(new_node)
    return (level, output_file_ids)

## This function creates one node of the reduce tree
def create_reduce_node(init_func, input_file_ids, output_file_ids):
    return init_func(flatten_list(input_file_ids) + output_file_ids)

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
