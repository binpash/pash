import os
import argparse
import sys
import pickle
import subprocess
import jsonpickle
import yaml

from ir import *
from json_ast import *
from impl import execute
from distr_back_end import distr_execute

## This file receives the name of a file that holds an IR, reads the
## IR, read some configuration file with node information, and then
## should make a distribution plan for it.

GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
if 'DISH_TOP' in os.environ:
    DISH_TOP = os.environ['DISH_TOP']
else:
    DISH_TOP = subprocess.run(GIT_TOP_CMD, capture_output=True,
            text=True).stdout.rstrip()
    
PARSER_BINARY = os.path.join(DISH_TOP, "parser/parse_to_json.native")


config = {}

def load_config():
    global config
    dish_config = {}
    CONFIG_KEY = 'distr_planner'

    ## TODO: allow this to be passed as an argument
    config_file_path = '{}/compiler/config.yaml'.format(DISH_TOP)
    with open(config_file_path) as config_file:
        dish_config = yaml.load(config_file, Loader=yaml.FullLoader)

    if not dish_config:
        raise Exception('No valid configuration could be loaded from {}'.format(config_file_path))

    if CONFIG_KEY not in dish_config:
        raise Exception('Missing `{}` config in {}'.format(CONFIG_KEY, config_file_path))

    config = dish_config[CONFIG_KEY]

## There are two ways to enter the distributed planner, either by
## calling dish (which straight away calls the distributed planner),
## or by calling the distributed planner with the name of an ir file
## to execute.
def main():
    ## Parse arguments
    args = parse_args()

    ## Call the main procedure
    optimize_script(args.input_ir, args.compile_optimize_only)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_ir", help="the file containing the dataflow graph to be optimized and executed")
    parser.add_argument("--compile_optimize_only",
                        help="only compile and optimize the input script and not execute it",
                        action="store_true")
    args = parser.parse_args()
    return args

def optimize_script(ir_filename, compile_optimize_only):
    global config
    if not config:
        load_config()

    with open(ir_filename, "rb") as ir_file:
        ir_node = pickle.load(ir_file)

    print("Retrieving IR: {} ...".format(ir_filename))
    shell_string = ast_to_shell(ir_node.ast)
    print(shell_string)

    print(ir_node)
    distributed_graph = naive_parallelize_stateless_nodes_bfs(ir_node, config['fan_out'], config['batch_size'])
    print(distributed_graph)
    # print("Parallelized graph:")
    # print(graph)

    ## Call the backend that executes the optimized dataflow graph
    output_script_path = config['optimized_script_filename']
    if(config['distr_backend']):
        distr_execute(distributed_graph, config['output_dir'], output_script_path,
                      config['output_optimized'], compile_optimize_only, config['nodes'])
    else:
        execute(distributed_graph.serialize_as_JSON(), config['output_dir'],
                output_script_path, config['output_optimized'], compile_optimize_only)

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
    #     input_file_ids = source_node.get_input_file_ids()
    #     ## TODO: Also split when we have more than one input file
    #     if(len(input_file_ids) == 1 and
    #        str(source_node.command) in commands_to_split_input):
    #         input_file_id = input_file_ids[0]
    #         input_file_id.split_resource(2, batch_size, fileIdGen)

    ## Starting from the sources of the graph, if they are stateless,
    ## duplicate the command as many times as the number of
    ## identifiers in its in_stream. Then connect their outputs in
    ## order to next command.
    workset= source_nodes
    while (len(workset) > 0):
        curr = workset.pop(0)

        next_nodes = graph.get_next_nodes(curr)
        workset += next_nodes

        ## Question: What does it mean for a command to have more
        ## than one next_node? Does it mean that it duplicates its
        ## output to all of them? Or does it mean that it writes
        ## some to the first and some to the second? Both are not
        ## very symmetric, but I think I would prefer the first.
        # print(curr, next_nodes)
        # assert(len(next_nodes) <= 1)

        ## If the current command is a split file, then we will
        ## recursively add a huge amount of splits, as several nodes
        ## are revisited due to the suboptimal adding of all new nodes
        ## to the workset.
        ##
        ## TODO: Remove hardcoded
        if(not str(curr.command) == "split_file"):
            for next_node in next_nodes:
                graph = split_command_input(next_node, graph, fileIdGen, fan_out, batch_size)
        graph, new_nodes = parallelize_cat(curr, graph, fileIdGen)

        ## Add new nodes to the workset depending on the optimization.
        ##
        ## WARNING: There is an assumption here that if there are new
        ## nodes there was an optimization that happened and these new
        ## nodes should ALL be added to the workset. Even if that is
        ## correct, that is certainly non-optimal.
        ##
        ## TODO: Fix that
        if(len(new_nodes) > 0):
            workset += new_nodes

    return graph


## TODO: Instead of setting children, we have to use the new way of
## having a cat command before the node so that we enable
## optimizations.

## Optimizes several commands by splitting its input
def split_command_input(curr, graph, fileIdGen, fan_out, batch_size):
    ## At the moment this only works for nodes that have one
    ## input. TODO: Extend it to work for nodes that have more than
    ## one input.
    ##
    ## TODO: Change the test to check if curr is instance of Cat and
    ## not check its name.
    previous_nodes = graph.get_previous_nodes(curr)
    if (curr.category in ["stateless", "pure"] and
        not str(curr.command) == "cat" and
        len(previous_nodes) == 1 and
        fan_out > 1):
        ## If the previous command is either a cat with one input, or
        ## if it something else
        previous_node = previous_nodes[0]
        if(not isinstance(previous_node, Cat) or
           (isinstance(previous_node, Cat) and
            len(previous_node.get_input_file_ids()) == 1)):

            if(not isinstance(previous_node, Cat)):
                input_file_ids = curr.get_input_file_ids()
                assert(len(input_file_ids) == 1)
                input_file_id = input_file_ids[0]
            else:
                ## If the previous node is a cat, we need its input
                input_file_ids = previous_node.get_input_file_ids()
                assert(len(input_file_ids) == 1)
                input_file_id = input_file_ids[0]

            ## First we have to make the split file commands.
            split_file_commands, output_fids = make_split_files(input_file_id, fan_out,
                                                                batch_size, fileIdGen)
            [graph.add_node(split_file_command) for split_file_command in split_file_commands]

            ## With the split file commands in place and their output
            ## fids (and this commands new input ids, we have to
            ## create a new Cat node (or modify the existing one) to
            ## have these as inputs, and connect its output to our
            ## input.

            ## Generate a new file id for the input of the current
            ## command.
            new_input_file_id = fileIdGen.next_file_id()
            new_cat = make_cat_node(output_fids, new_input_file_id)
            graph.add_node(new_cat)

            if(not isinstance(previous_node, Cat)):
                ## Change the current node's input with the new_input
                index = curr.find_file_id_in_in_stream(input_file_id)
                chunk = curr.in_stream[index]
                curr.set_file_id(chunk, new_input_file_id)
            else:
                ## Change the current node's input with the new_input
                curr_input_file_ids = curr.get_input_file_ids()
                assert(len(curr_input_file_ids) == 1)
                curr_input_file_id = curr_input_file_ids[0]
                chunk = curr.find_file_id_in_in_stream(curr_input_file_id)
                curr.set_file_id(chunk, new_input_file_id)
                graph.remove_node(previous_node)
    return graph

## If the current command is a cat, and is followed by a node that
## is either stateless or pure parallelizable, commute the cat
## after the node.
def parallelize_cat(curr, graph, fileIdGen):
    new_nodes_for_workset = []
    if(isinstance(curr, Cat)):
        next_nodes_and_edges = graph.get_next_nodes_and_edges(curr)

        ## Cat can only have one output
        assert(len(next_nodes_and_edges) <= 1)
        if(len(next_nodes_and_edges) == 1):
            next_node = next_nodes_and_edges[0][0]

            ## Get cat inputs and output. Note that there is only one
            ## output.
            cat_input_file_ids = curr.get_input_file_ids()
            cat_output_file_id = next_nodes_and_edges[0][1]
            graph, new_nodes = parallelize_command(next_node, cat_input_file_ids,
                                                   cat_output_file_id, graph, fileIdGen)
            new_nodes_for_workset += new_nodes

            ## If there are no new nodes, it means that no
            ## parallelization happened. If there were, we should
            ## delete the original cat.
            if(len(new_nodes) > 0):
                graph.remove_node(curr)

    return graph, new_nodes_for_workset

def parallelize_command(curr, new_input_file_ids, old_input_file_id, graph, fileIdGen):
    ## If the next command is stateless or pure parallelizable and has more
    ## than one inputs, it can be parallelized.
    new_nodes = []
    if ((curr.category == "stateless" or curr.is_pure_parallelizable())
        and len(new_input_file_ids) > 1):

        ## We assume that every stateless and pure command has
        ## one output_file_id for now.
        ##
        ## TODO: Check if this can be lifted. This seems to be
        ## connected to assertion regarding the next_nodes.
        node_out_file_ids = curr.get_output_file_ids()
        assert(len(node_out_file_ids) == 1)
        node_out_file_id = node_out_file_ids[0]

        ## Get the file names of the outputs of the map commands. This
        ## differs if the command is stateless, pure that can be
        ## written as a map and a reduce, and a pure that can be
        ## written as a generalized map and reduce.
        new_output_file_ids = curr.get_map_output_files(new_input_file_ids, fileIdGen)

        ## Make a merge command that joins the results of all the
        ## duplicated commands
        if(curr.is_pure_parallelizable()):
            merge_commands = create_merge_commands(curr, new_output_file_ids,
                                                   node_out_file_id, fileIdGen)
            for merge_command in merge_commands:
                graph.add_node(merge_command)
            new_nodes += merge_commands

        ## For each new input and output file id, make a new command
        new_commands = curr.duplicate(old_input_file_id, new_input_file_ids,
                                      new_output_file_ids, fileIdGen)


        new_nodes += new_commands

        # print("New commands:")
        # print(new_commands)
        graph.remove_node(curr)
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
    return graph, new_nodes


## Creates a merge command for all pure commands that can be
## parallelized using a map and a reduce/merge step
def create_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen):
    if(str(curr.command) == "sort"):
        return create_sort_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen)
    elif(str(curr.command) == "bigrams_aux"):
        return create_bigram_aux_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen)
    elif(str(curr.command) == "alt_bigrams_aux"):
        return create_alt_bigram_aux_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen)
    else:
        assert(False)

## TODO: These must be generated using some file information
##
## TODO: Find a better place to put these functions
def create_sort_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen):
    tree = create_reduce_tree(SortGReduce, new_output_file_ids, out_edge_file_id, fileIdGen)
    return tree

def create_bigram_aux_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen):
    tree = create_reduce_tree(BigramGReduce, new_output_file_ids, out_edge_file_id, fileIdGen)
    return tree

def create_alt_bigram_aux_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen):
    tree = create_reduce_tree(AltBigramGReduce, new_output_file_ids, out_edge_file_id, fileIdGen)
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

if __name__ == "__main__":
    main()
