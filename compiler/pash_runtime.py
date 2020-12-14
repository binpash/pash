import cProfile
import os
import argparse
import sys
import pickle
import subprocess
import jsonpickle
from datetime import datetime

from ir import *
from ast_to_ir import compile_asts
from json_ast import *
from ir_to_ast import to_shell
from parse import from_ir_to_shell
from util import *
import config

from definitions.ir.nodes.alt_bigram_g_reduce import *
from definitions.ir.nodes.bigram_g_map import *
from definitions.ir.nodes.bigram_g_reduce import *
from definitions.ir.nodes.sort_g_reduce import *
from definitions.ir.nodes.eager import *
from definitions.ir.nodes.pash_split import *

runtime_config = {}
## There are two ways to enter the distributed planner, either by
## calling pash_runtime.py (which straight away calls the distributed planner),
## or by calling the distributed planner with the name of an ir file
## to execute.
def main():
    ## Parse arguments
    args = parse_args()
    config.pash_args = args

    ## TODO: Instead of just calling optimize we first need to call compile.

    ## Call the main procedure
    compile_optimize_script(args.input_ir, args.compiled_script_file, args)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("compiled_script_file", 
                        help="the file in which to output the compiled script")
    parser.add_argument("input_ir", 
                        help="the file containing the dataflow graph to be optimized and executed")
    parser.add_argument("--var_file",
                        help="determines the path of a file containing all shell variables.",
                        default=None)
    config.add_common_arguments(parser)
    args = parser.parse_args()
    return args

def compile_optimize_script(ir_filename, compiled_script_file, args):
    global runtime_config
    if not config.config:
        config.load_config(args.config_path)

    ## Load annotations
    config.annotations = load_annotation_files(config.config['distr_planner']['annotations_dir'])

    runtime_config = config.config['distr_planner']

    log("Retrieving candidate DF region: {} ... ".format(ir_filename), end='')
    with open(ir_filename, "rb") as ir_file:
        candidate_df_region = pickle.load(ir_file)
    log("Done!")

    ## Read any shell variables files if present
    config.read_vars_file(args.var_file)

    ## Compile the candidate DF regions
    compilation_start_time = datetime.now()
    asts_and_irs = compile_candidate_df_region(candidate_df_region, config.config)
    compilation_end_time = datetime.now()
    print_time_delta("Compilation", compilation_start_time, compilation_end_time, args)

    ## Optimize all the IRs that can be optimized
    if(args.compile_only):
        optimized_asts_and_irs = asts_and_irs
    else:
        optimized_asts_and_irs = optimize_irs(asts_and_irs, args)

    ## TODO: Normally this could return more than one compiled ASTs (containing IRs in them). 
    ##       To correctly handle that we would need to really replace the optimized IRs
    ##       with the final parallel corresponding scripts.
    ##
    ##       However, for now we just assume that there is one IR that we can execute as is.
    ##
    ## TODO: This might bite us with the quick-abort. 
    ##       It might complicate things having a script whose half is compiled to a graph and its other half not.
    assert(len(optimized_asts_and_irs) == 1)
    optimized_ast_or_ir = optimized_asts_and_irs[0]

    ## Call the backend that executes the optimized dataflow graph
    output_script_path = runtime_config['optimized_script_filename']
    ## TODO: Should never be the case for now. This is obsolete.
    assert(not runtime_config['distr_backend'])
    
    ## If the candidate DF region was indeed a DF region then we have an IR
    ## which should be translated to a parallel script.
    if(isinstance(optimized_ast_or_ir, IR)):
        script_to_execute = to_shell(optimized_ast_or_ir, 
                                     runtime_config['output_dir'], args)

        ## TODO: Merge this write with the one below. Maybe even move this logic in `pash_runtime.sh`
        ## Output the optimized shell script for inspection
        if(args.output_optimized):
            with open(output_script_path, "w") as output_script_file:
                log("Optimized script:")
                log(script_to_execute)
                output_script_file.write(script_to_execute)

        with open(compiled_script_file, "w") as f:
                f.write(script_to_execute)

    else:
        ## Instead of outputing the script here, we just want to exit with a specific exit code
        ## TODO: Figure out the code and save it somewhere
        exit(120)
    

def compile_candidate_df_region(candidate_df_region, config):
    ## This is for the files in the IR
    fileIdGen = FileIdGen()

    ## If the candidate DF region is not from the top level then 
    ## it won't be a list and thus we need to make it into a list to compile it.
    if(not isinstance(candidate_df_region, list)):
        candidate_df_region = [candidate_df_region]

    ## Compile the asts
    ## TODO: Since compilation happens at runtime, we can now expand everything accordingly.
    ##       We can do that using a shell for start: 
    ##         if a word is safe to expand, then call a shell to expand it.
    compiled_asts = compile_asts(candidate_df_region, fileIdGen, config)

    return compiled_asts

def optimize_irs(asts_and_irs, args):
    global runtime_config

    optimization_start_time = datetime.now()

    optimized_asts_and_irs = []
    for ast_or_ir in asts_and_irs:
        if(isinstance(ast_or_ir, IR)):
            # log(ir_node)
            # with cProfile.Profile() as pr:
            distributed_graph = naive_parallelize_stateless_nodes_bfs(ast_or_ir, args.split_fan_out,
                                                                      runtime_config['batch_size'])
            # pr.print_stats()
            # log(distributed_graph)

            if(not args.no_eager):
                eager_distributed_graph = add_eager_nodes(distributed_graph)
            else:
                eager_distributed_graph = distributed_graph

            ## Print statistics of output nodes
            print_graph_statistics(eager_distributed_graph)
            # log(eager_distributed_graph)

            optimized_asts_and_irs.append(eager_distributed_graph)
        else:
            optimized_asts_and_irs.append(ast_or_ir)

    optimization_end_time = datetime.now()
    print_time_delta("Optimization", optimization_start_time, optimization_end_time, args)

    return optimized_asts_and_irs


def print_graph_statistics(graph):
    total_nodes = graph.nodes
    cat_nodes = [node for node in total_nodes if isinstance(node, Cat)]
    eager_nodes = [node for node in total_nodes if isinstance(node, Eager)]
    log("Total nodes after optimization:", len(total_nodes))
    log(" -- out of which:")
    log("Cat nodes:", len(cat_nodes))
    log("Eager nodes:", len(eager_nodes))

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
    # log("Source nodes:")
    # log(source_nodes)

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
    visited = set()
    while (len(workset) > 0):

        curr = workset.pop(0)
        ## Node must not be in visited, but must also be in the graph
        ## (because it might have been deleted after some
        ## optimization).
        if(not curr in visited
           and curr in graph.nodes):
            visited.add(curr)
            next_nodes = graph.get_next_nodes(curr)
            workset += next_nodes

            ## Question: What does it mean for a command to have more
            ## than one next_node? Does it mean that it duplicates its
            ## output to all of them? Or does it mean that it writes
            ## some to the first and some to the second? Both are not
            ## very symmetric, but I think I would prefer the first.
            # log(curr, next_nodes)
            # assert(len(next_nodes) <= 1)

            graph, new_nodes = parallelize_cat(curr, graph, fileIdGen, fan_out, batch_size)

            ## Add new nodes to the workset depending on the optimization.
            ##
            ## WARNING: There is an assumption here that if there are new
            ## nodes there was an optimization that happened and these new
            ## nodes should ALL be added to the workset. Even if that is
            ## correct, that is certainly non-optimal.
            ##
            ## TODO: Fix that
            if(len(new_nodes) > 0):
                # log("New nodes:", new_nodes)
                workset += new_nodes

    return graph


## TODO: Instead of setting children, we have to use the new way of
## having a cat command before the node so that we enable
## optimizations.

## Optimizes several commands by splitting its input
def split_command_input(curr, previous_node, graph, fileIdGen, fan_out, batch_size):
    ## At the moment this only works for nodes that have one
    ## input. TODO: Extend it to work for nodes that have more than
    ## one input.
    ##
    ## TODO: Change the test to check if curr is instance of Cat and
    ## not check its name.
    number_of_previous_nodes = curr.get_number_of_inputs()
    assert(curr.category == "stateless" or curr.is_pure_parallelizable())
    # curr.category in ["stateless", "pure"] and

    new_cat = None
    if (not isinstance(curr, Cat)
        and number_of_previous_nodes == 1
        and fan_out > 1):
        ## If the previous command is either a cat with one input, or
        ## if it something else
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
    return (graph, new_cat)

## If the current command is a cat, and is followed by a node that
## is either stateless or pure parallelizable, commute the cat
## after the node.
def parallelize_cat(curr, graph, fileIdGen, fan_out, batch_size):
    new_nodes_for_workset = []

    ## Get next nodes in the graph
    next_nodes_and_edges = graph.get_next_nodes_and_edges(curr)

    ## If there is only one node afterwards (meaning that we reached a
    ## thin part of the graph), we need to try to parallelize
    if(len(next_nodes_and_edges) == 1):
        next_node = next_nodes_and_edges[0][0]

        ## If the next node can be parallelized, then we should try to
        ## parallelize
        if(next_node.category == "stateless" or next_node.is_pure_parallelizable()):

            ## If the current node is not a cat, it means that we need
            ## to generate a cat using a split
            if(not isinstance(curr, Cat)):
                graph, new_cat = split_command_input(next_node, curr, graph, fileIdGen, fan_out, batch_size)
                if(not new_cat is None):
                    curr = new_cat
                    next_nodes_and_edges = graph.get_next_nodes_and_edges(curr)
                    assert(len(next_nodes_and_edges) == 1)
                    next_node = next_nodes_and_edges[0][0]

            ## If curr is cat, it means that split suceeded, or it was
            ## already a cat. In any case, we can proceed with the
            ## parallelization
            if(isinstance(curr, Cat)):

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
    assert(curr.category == "stateless" or curr.is_pure_parallelizable())
    ## If the next command is stateless or pure parallelizable and has more
    ## than one inputs, it can be parallelized.
    new_nodes = []
    if (len(new_input_file_ids) > 1):

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

            ## Add the merge commands in the graph
            for merge_command in merge_commands:
                graph.add_node(merge_command)
            new_nodes += merge_commands

        ## For each new input and output file id, make a new command
        new_commands = duplicate(curr, old_input_file_id, new_input_file_ids,
                                      new_output_file_ids, fileIdGen)


        new_nodes += new_commands

        # log("New commands:")
        # log(new_commands)
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

def duplicate(command_node, old_input_file_id, new_input_file_ids, new_output_file_ids, fileIdGen):
    assert(command_node.category == "stateless" or command_node.is_pure_parallelizable())
    if (command_node.category == "stateless"):
        return stateless_duplicate(command_node, old_input_file_id, new_input_file_ids, new_output_file_ids)
    elif (command_node.is_pure_parallelizable()):
        return pure_duplicate(command_node, old_input_file_id, new_input_file_ids, new_output_file_ids, fileIdGen)

    ## This should be unreachable
    log("Unreachable code reached :(")
    assert(False)

def stateless_duplicate(command_node, old_input_file_id, input_file_ids, output_file_ids):
    assert(command_node.category == "stateless")

    out_edge_file_ids = command_node.get_output_file_ids()
    assert(len(out_edge_file_ids) == 1)
    out_edge_file_id = out_edge_file_ids[0]

    ## Make a new cat command to add after the current command.
    new_cat = make_cat_node(output_file_ids, out_edge_file_id)

    in_out_file_ids = zip(input_file_ids, output_file_ids)

    new_commands = [command_node.find_in_out_and_make_duplicate_command(old_input_file_id, in_fid,
        out_edge_file_id, out_fid) for in_fid, out_fid in in_out_file_ids]

    return new_commands + [new_cat]

def pure_duplicate(command_node, old_input_file_id, input_file_ids, output_file_ids, fileIdGen):
    assert(command_node.is_pure_parallelizable())

    in_out_file_ids = zip(input_file_ids, output_file_ids)

    simple_map_pure_commands = ["sort",
                                "alt_bigrams_aux",
                                "uniq"]

    ## This is the category of all commands that don't need a
    ## special generalized map
    if(str(command_node.command) in simple_map_pure_commands):

        ## make_duplicate_command duplicates a node based on its
        ## output file ids, so we first need to assign them.
        new_output_file_ids = [fids[0] for fids in output_file_ids]

        new_commands = [command_node.find_in_out_and_make_duplicate_command(old_input_file_id, in_fid,
            command_node.stdout, out_fids[0]) for in_fid, out_fids in in_out_file_ids]
    elif(str(command_node.command) == "bigrams_aux"):
        new_commands = [BigramGMap([in_fid] + out_fids)
                        for in_fid, out_fids in in_out_file_ids]
    else:
        ## This should be unreachable
        log("Unreachable code reached :(")
        assert(False)

    return new_commands

## Creates a merge command for all pure commands that can be
## parallelized using a map and a reduce/merge step
def create_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen):
    if(str(curr.command) == "sort"):
        return create_sort_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen)
    elif(str(curr.command) == "bigrams_aux"):
        return create_bigram_aux_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen)
    elif(str(curr.command) == "alt_bigrams_aux"):
        return create_alt_bigram_aux_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen)
    elif(str(curr.command) == "uniq"):
        return create_uniq_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen)
    else:
        raise NotImplementedError()

## TODO: These must be generated using some file information
##
## TODO: Find a better place to put these functions
def create_sort_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen):
    tree = create_reduce_tree(lambda file_ids: SortGReduce(curr, file_ids),
                              new_output_file_ids, out_edge_file_id, fileIdGen)
    return tree

def create_bigram_aux_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen):
    tree = create_reduce_tree(lambda file_ids: BigramGReduce(curr, file_ids),
                              new_output_file_ids, out_edge_file_id, fileIdGen)
    return tree

def create_alt_bigram_aux_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen):
    tree = create_reduce_tree(lambda file_ids: AltBigramGReduce(curr, file_ids),
                              new_output_file_ids, out_edge_file_id, fileIdGen)
    return tree

## Instead of creating a tree, we just create a single level reducer for uniq
def create_uniq_merge_commands(curr, new_output_file_ids, out_edge_file_id, fileIdGen):
    ## Add a cat node that takes all inputs
    intermediate_file_id = fileIdGen.next_file_id()
    new_cat = make_cat_node(flatten_list(new_output_file_ids), intermediate_file_id)

    ## Add a uniq node after it
    command = string_to_argument("uniq")
    ## TODO: Pass the options of `curr` correctly
    options = []
    in_stream = ["stdin"]
    out_stream = ["stdout"]
    stdin = intermediate_file_id
    stdout = out_edge_file_id
    opt_indices = []
    ## TODO: Maybe hack this and change this to pure to not reparallelize
    category = "pure_parallelizable"
    ## TODO: Fill the AST
    ast = None
    uniq_node = Command(ast, command, options, in_stream, out_stream, opt_indices, 
                        category, stdin=stdin, stdout=stdout)
    
    return [new_cat, uniq_node]

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
    ##
    ## TODO: If the union below is done in the reverse order, then
    ## redirections are not transfered. FIX this disgusting union-find
    ## structure.
    output_file_id.union(curr_file_ids[0][0])

    ## Drain the final auxiliary outputs
    final_auxiliary_outputs = curr_file_ids[0][1:]
    drain_file_id = fileIdGen.next_file_id()
    drain_file_id.set_resource('/dev/null')
    drain_cat_commands = [make_cat_node([final_auxiliary_output], drain_file_id)
                          for final_auxiliary_output in final_auxiliary_outputs]
    return (tree + drain_cat_commands)


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



## This function adds eager nodes wherever the width of graph is
## becoming smaller.
def add_eager_nodes(graph):
    source_nodes = graph.source_nodes()
    # log("Source nodes:")
    # log(source_nodes)

    eager_exec_path = '{}/{}'.format(config.PASH_TOP, runtime_config['eager_executable_path'])
    ## Generate a fileIdGen from a graph, that doesn't clash with the
    ## current graph fileIds.
    fileIdGen = graph.get_file_id_gen()
    intermediateFileIdGen = FileIdGen(0, runtime_config['eager_intermediate_prefix'])

    ## Get the next nodes
    workset = [node for source_node in source_nodes for node in graph.get_next_nodes(source_node)]
    visited = set()
    while (len(workset) > 0):
        curr = workset.pop(0)
        if (not curr in visited):
            visited.add(curr)
            next_nodes = graph.get_next_nodes(curr)
            workset += next_nodes

            ## Add eager nodes if the node has more than one input
            curr_input_file_ids = curr.get_input_file_ids()
            if (len(curr_input_file_ids) > 1):
                new_input_file_ids = [fileIdGen.next_file_id() for _ in curr_input_file_ids]
                intermediate_file_ids = [intermediateFileIdGen.next_file_id() for _ in curr_input_file_ids]

                file_ids = list(zip(curr_input_file_ids, new_input_file_ids, intermediate_file_ids))
                eager_nodes = [make_eager_node(old_file_id, new_file_id,
                                               intermediate_file_id, eager_exec_path)
                               for old_file_id, new_file_id, intermediate_file_id in file_ids]

                for eager_node in eager_nodes:
                    graph.add_node(eager_node)

                ## Update input file ids
                for curr_input_file_id, new_input_file_id, _ in file_ids:
                    chunk_i = curr.find_file_id_in_in_stream(curr_input_file_id)
                    chunk = curr.in_stream[chunk_i]
                    curr.set_file_id(chunk, new_input_file_id)

            ## TODO: Make sure that we don't add duplicate eager nodes

            ## TODO: Refactor this and the above as they are very symmetric
            if(isinstance(curr, Split)):
                curr_output_file_ids = curr.get_output_file_ids()
                output_file_ids_to_eager = curr_output_file_ids[:-1]

                new_output_file_ids = [fileIdGen.next_file_id()
                                       for _ in output_file_ids_to_eager]
                intermediate_file_ids = [intermediateFileIdGen.next_file_id()
                                         for _ in output_file_ids_to_eager]

                file_ids = list(zip(output_file_ids_to_eager, new_output_file_ids, intermediate_file_ids))
                eager_nodes = [make_eager_node(new_file_id, old_file_id,
                                               intermediate_file_id, eager_exec_path)
                               for old_file_id, new_file_id, intermediate_file_id in file_ids]

                for eager_node in eager_nodes:
                    graph.add_node(eager_node)

                ## Update input file ids
                for curr_output_file_id, new_output_file_id, _ in file_ids:
                    chunk_i = curr.find_file_id_in_out_stream(curr_output_file_id)
                    chunk = curr.out_stream[chunk_i]
                    curr.set_file_id(chunk, new_output_file_id)



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
