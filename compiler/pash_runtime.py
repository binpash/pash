import argparse
import sys
import pickle
import traceback
from datetime import datetime

from annotation_generation_new.datatypes.parallelizability.AggregatorKind import AggregatorKindEnum

import config
from ir import *
from ast_to_ir import compile_asts
from json_ast import *
from ir_to_ast import to_shell
from pash_graphviz import maybe_generate_graphviz
from util import *

from definitions.ir.aggregator_node import *

from definitions.ir.dfg_node import DFGNode
from definitions.ir.nodes.eager import *
from definitions.ir.nodes.pash_split import *

import definitions.ir.nodes.r_merge as r_merge
import definitions.ir.nodes.r_split as r_split
import definitions.ir.nodes.r_unwrap as r_unwrap
import definitions.ir.nodes.dgsh_tee as dgsh_tee
import definitions.ir.nodes.dfs_split_reader as dfs_split_reader
# Distirbuted Exec
import dspash.hdfs_utils as hdfs_utils 

runtime_config = {}
## We want to catch all exceptions here so that they are logged correctly
## and not just printed to the stderr.
def main():
    try:
        main_body()
    except Exception:
        log("Compiler failed, no need to worry, executing original script...")
        log(traceback.format_exc())
        sys.exit(1)

def main_body():
    global runtime_config

    ## Parse arguments
    args = parse_args()
    config.pash_args = args

    ## Load the configuration
    if not config.config:
        config.load_config(args.config_path)

    ## Load annotations
    config.annotations = load_annotation_files(config.config['distr_planner']['annotations_dir'])

    runtime_config = config.config['distr_planner']

    ## Read any shell variables files if present
    config.read_vars_file(args.var_file)

    log("Input:", args.input_ir, "Compiled file:", args.compiled_script_file)

    ## Call the main procedure
    compiler_config = CompilerConfig(args.width)
    ast_or_ir = compile_optimize_output_script(args.input_ir, args.compiled_script_file, args, compiler_config)
    maybe_generate_graphviz(ast_or_ir, args)


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
    args, unknown_args = parser.parse_known_args()
    return args

## TODO: Add more fields from args in this
class CompilerConfig:
    def __init__(self, width):
        self.width = width
    
    def __repr__(self):
        return f'CompilerConfig(Width:{self.width})'

def compile_ir(ir_filename, compiled_script_file, args, compiler_config):
    """
    Return IR object for compilation success. None otherwise.
    """
    ret = None
    try:
        ret = compile_optimize_output_script(ir_filename, compiled_script_file, args, compiler_config)
    except Exception as e:
        log("WARNING: Exception caught:", e)
        # traceback.print_exc()

    return ret

def compile_optimize_output_script(ir_filename, compiled_script_file, args, compiler_config):
    global runtime_config
    
    ret = None

    ## Load the df_region from a file
    candidate_df_region = load_df_region(ir_filename)
    
    ## Compile it
    optimized_ast_or_ir = compile_optimize_df_region(candidate_df_region, args, compiler_config)

    ## Call the backend that executes the optimized dataflow graph
    ## TODO: Should never be the case for now. This is obsolete.
    assert(not runtime_config['distr_backend'])

    ## If the candidate DF region was indeed a DF region then we have an IR
    ## which should be translated to a parallel script.
    if(isinstance(optimized_ast_or_ir, IR)):
        if args.distributed_exec:
            ir_filename = ptempfile()
            script_to_execute = f"$PASH_TOP/compiler/dspash/remote_exec_graph.sh {ir_filename}\n"
            ## This might not be needed anymore (since the output script is output anyway)
            ## TODO: This is probably useless, remove
            maybe_log_optimized_script(script_to_execute, args)

            with open(ir_filename, "wb") as f:
                obj = (optimized_ast_or_ir, config.config['shell_variables'])
                pickle.dump(obj, f)
        else:
            script_to_execute = to_shell(optimized_ast_or_ir, args)
            
        log("Optimized script saved in:", compiled_script_file)
        with open(compiled_script_file, "w") as f:
            f.write(script_to_execute)
    
        ret = optimized_ast_or_ir
    else:
        raise Exception("Script failed to compile!")
    
    return ret

def load_df_region(ir_filename):
    log("Retrieving candidate DF region: {} ... ".format(ir_filename), end='')
    with open(ir_filename, "rb") as ir_file:
        candidate_df_region = pickle.load(ir_file)
    log("Done!")
    return candidate_df_region

def compile_optimize_df_region(df_region, args, compiler_config):
    ## Compile the candidate DF regions
    compilation_start_time = datetime.now()
    asts_and_irs = compile_candidate_df_region(df_region, config.config)
    compilation_end_time = datetime.now()
    print_time_delta("Compilation", compilation_start_time, compilation_end_time, args)

    ## Optimize all the IRs that can be optimized
    if(args.no_optimize):
        optimized_asts_and_irs = asts_and_irs
    else:
        optimized_asts_and_irs = optimize_irs(asts_and_irs, args, compiler_config)

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
    
    return optimized_ast_or_ir

def maybe_log_optimized_script(script_to_execute, args):
    ## TODO: Merge this write with the one below. Maybe even move this logic in `pash_runtime.sh`
    ## Output the optimized shell script for inspection
    if(args.output_optimized):
        output_script_path = runtime_config['optimized_script_filename']
        with open(output_script_path, "w") as output_script_file:
            log("Optimized script:")
            log(script_to_execute)
            output_script_file.write(script_to_execute)

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

## TODO: Switch args to compiler_config
def optimize_irs(asts_and_irs, args, compiler_config):
    global runtime_config

    optimization_start_time = datetime.now()

    optimized_asts_and_irs = []
    for ast_or_ir in asts_and_irs:
        if(isinstance(ast_or_ir, IR)):
            ## Assert that the graph that was returned from compilation is valid
            assert(ast_or_ir.valid())

            # log(ir_node)
            # with cProfile.Profile() as pr:
            distributed_graph = choose_and_apply_parallelizing_transformations(ast_or_ir, compiler_config.width,
                                                                      runtime_config['batch_size'],
                                                                      args.no_cat_split_vanish,
                                                                      args.r_split, args.r_split_batch_size)
            # pr.print_stats()

            # Eagers are added in remote notes when using distributed exec
            if(not args.no_eager and not args.distributed_exec): 
                eager_distributed_graph = add_eager_nodes(distributed_graph, args.dgsh_tee)
            else:
                eager_distributed_graph = distributed_graph

            ## Assert that the graph stayed valid after all transformations
            assert(eager_distributed_graph.valid())

            ## Print statistics of output nodes
            print_graph_statistics(eager_distributed_graph)

            optimized_asts_and_irs.append(eager_distributed_graph)
        else:
            optimized_asts_and_irs.append(ast_or_ir)

    optimization_end_time = datetime.now()
    print_time_delta("Optimization", optimization_start_time, optimization_end_time, args)

    return optimized_asts_and_irs


def print_graph_statistics(graph):
    total_nodes = graph.nodes
    cat_nodes = [node for node in total_nodes.values() if isinstance(node, Cat)]
    eager_nodes = [node for node in total_nodes.values() if isinstance(node, Eager)]
    log("Total nodes after optimization:", len(total_nodes))
    log(" -- out of which:")
    log("Cat nodes:", len(cat_nodes))
    log("Eager nodes:", len(eager_nodes))


def choose_and_apply_parallelizing_transformations(graph, fan_out, batch_size, no_cat_split_vanish,
                                                   r_split_flag, r_split_batch_size):
    parallelizer_map = choose_parallelizing_transformations(graph, r_split_flag)
    apply_parallelizing_transformations(graph, parallelizer_map, fan_out, batch_size, no_cat_split_vanish,
                                                   r_split_flag, r_split_batch_size)
    return graph


def choose_parallelizing_transformations(graph, r_split_flag): # shall return map
    source_node_ids = graph.source_nodes()
    parallelizer_map = {}
    workset = source_node_ids
    visited = set()
    # We apply a modified BFS such that we ensure that we know which parallelizer was chosen for all previous nodes
    # and assume that the decision for any subsequent node will exploit any potential synergy effects
    while (len(workset) > 0):
        curr_id = workset.pop(0)
        assert(isinstance(curr_id, int))
        all_previous_nodes_visited = all(prev in visited for prev in graph.get_previous_nodes(curr_id))
        if not all_previous_nodes_visited:
            workset.append(curr_id)
        elif not curr_id in visited:
            next_node_ids = graph.get_next_nodes(curr_id)
            workset += next_node_ids
            parallelizer_map[curr_id] = choose_parallelizing_transformation(curr_id, graph, r_split_flag)
            visited.add(curr_id)
    return parallelizer_map


def choose_parallelizing_transformation(curr_id, graph, r_split_flag): # shall return map entry
    # here we can implement more sophisticated techniques to decide how to parallelize
    curr = graph.get_node(curr_id)
    # we ignore `r_split_flag` here as we want to exploit r_merge followed by commutative command
    # which only works if the a parallelizer for the latter is chosen (sort does not have RR-parallelizer)
    # we prioritize round robin over round robin with unwrap over consecutive chunks:
    list_all_parallelizers_in_priority = [curr.get_option_implemented_round_robin_parallelizer(),
                                          curr.get_option_implemented_round_robin_with_unwrap_parallelizer(),
                                          curr.get_option_implemented_consecutive_chunks_parallelizer()]
    return next((item for item in list_all_parallelizers_in_priority if item is not None), None)
    # When `r_split_flag` should be used:
    # if r_split_flag:
    #     option_parallelizer = curr.get_option_implemented_round_robin_parallelizer()
    # else:
    #     option_parallelizer = curr.get_option_implemented_consecutive_chunks_parallelizer()
    # return option_parallelizer


def apply_parallelizing_transformations(graph, parallelizer_map, fan_out, batch_size, no_cat_split_vanish,
                                        r_split_flag, r_split_batch_size):
    fileIdGen = graph.get_file_id_gen()
    node_id_non_none_parallelizer_list = [(node_id, parallelizer) for (node_id, parallelizer) in parallelizer_map.items()
                                                                  if parallelizer is not None]
    for (node_id, parallelizer) in node_id_non_none_parallelizer_list:
        graph.apply_parallelization_to_node(node_id, parallelizer, fileIdGen, fan_out,
                                            batch_size, no_cat_split_vanish, r_split_batch_size)
## This is a simplistic planner, that pushes the available
## parallelization from the inputs in file stateless commands. The
## planner starts from the sources of the graph, and pushes
## file parallelization as far as possible.
##
## It returns a maximally expanded (regarding files) graph, that can
## be scheduled depending on the available computational resources.
def naive_parallelize_stateless_nodes_bfs(graph, fan_out, batch_size, no_cat_split_vanish,
                                          r_split_flag, r_split_batch_size):
    assert(False)
    source_node_ids = graph.source_nodes()

    ## Generate a fileIdGen from a graph, that doesn't clash with the
    ## current graph fileIds.
    fileIdGen = graph.get_file_id_gen()

    ## Starting from the sources of the graph traverse the whole graph using a
    ## node_id workset. Every iteration we add the next nodes to the workset as
    ## well as any newly added nodes due to optimizations.
    workset = source_node_ids
    visited = set()
    while (len(workset) > 0):
        curr_id = workset.pop(0)
        assert(isinstance(curr_id, int))
        ## Node must not be in visited, but must also be in the graph
        ## (because it might have been deleted after some
        ## optimization).
        if(not curr_id in visited
           and curr_id in graph.nodes):
            # log("Curr id:", curr_id)
            visited.add(curr_id)
            next_node_ids = graph.get_next_nodes(curr_id)
            workset += next_node_ids

            # function application has side effects on graphs
            new_nodes = parallelize_node(curr_id, graph, fileIdGen,
                                         fan_out, batch_size, no_cat_split_vanish,
                                         r_split_flag, r_split_batch_size)

            ## Assert that the graph stayed valid after the transformation
            ## TODO: Do not run this everytime in the loop if we are not in debug mode.
            # log("Graph nodes:", graph.nodes)
            # log("Graph edges:", graph.edges)
            # assert(graph.valid())

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
                workset += [node.get_id() for node in new_nodes]

    return graph


## Optimizes several commands by splitting its input
def split_command_input(curr, graph, fileIdGen, fan_out, _batch_size, r_split_flag, r_split_batch_size):
    assert(curr.is_parallelizable())
    assert(not isinstance(curr, Cat))
    assert(fan_out > 1)

    ## At the moment this only works for nodes that have one standard input.
    standard_input_ids = curr.get_standard_inputs()
    new_merger = None
    if (len(standard_input_ids) == 1):
        ## If the previous command is either a cat with one input, or
        ## if it something else
        
        input_id = standard_input_ids[0]

        ## First we have to make the split file commands.
        split_file_commands, output_fids = make_split_files(input_id, fan_out, fileIdGen, r_split_flag, r_split_batch_size)
        for output_fid in output_fids:
            output_fid.make_ephemeral()
            graph.add_edge(output_fid)

        for split_file_command in split_file_commands:
            graph.add_node(split_file_command)

        ## With the split file commands in place and their output
        ## fids (and this commands new input ids, we have to
        ## create a new Cat node (or modify the existing one) to
        ## have these as inputs, and connect its output to our
        ## input.

        ## Generate a new file id for the input of the current command.
        new_input_fid = fileIdGen.next_file_id()
        new_input_fid.make_ephemeral()
        graph.add_edge(new_input_fid)
        new_input_id = new_input_fid.get_ident()

        output_ids = [fid.get_ident() for fid in output_fids]

        ## If we add r_split, then the merger is actually an r_merge
        if(r_split_flag):
            new_merger = r_merge.make_r_merge_node(output_ids, new_input_id)
        else:
            new_merger = make_cat_node(output_ids, new_input_id)
        graph.add_node(new_merger)

        ## Replace the previous input edge with the new input edge that is after the cat.
        curr.replace_edge(input_id, new_input_id)
        graph.set_edge_to(new_input_id, curr.get_id())

        # log("graph nodes:", graph.nodes)
        # log("graph edges:", graph.edges)

    return new_merger

def split_hdfs_cat_input(hdfs_cat, next_node, graph, fileIdGen):
    """
    Replaces hdfs cat with a cat per block, each cat uses has an HDFSResource input fid
    Returns: A normal Cat that merges the blocks (will be removed when parallizing next_node)
    """
    assert(isinstance(hdfs_cat, HDFSCat))

    ## At the moment this only works for nodes that have one standard input.
    if len(next_node.get_standard_inputs()) != 1:
        return

    hdfscat_input_id = hdfs_cat.get_standard_inputs()[0]
    hdfs_fid = graph.get_edge_fid(hdfscat_input_id)
    hdfs_filepath = str(hdfs_fid.get_resource())
    output_ids = []

    # Create a cat command per file block
    file_config = hdfs_utils.get_file_config(hdfs_filepath)
    dummy_config_path = ptempfile() # Dummy config file, should be updated by workers
    for split_num, block in enumerate(file_config.blocks):
        resource = DFSSplitResource(file_config.dumps(), dummy_config_path, split_num, block.hosts)
        block_fid = fileIdGen.next_file_id()
        block_fid.set_resource(resource)
        graph.add_edge(block_fid)

        output_fid = fileIdGen.next_file_id()
        output_fid.make_ephemeral()
        output_ids.append(output_fid.get_ident())
        graph.add_edge(output_fid)

        split_reader_node = dfs_split_reader.make_dfs_split_reader_node([block_fid.get_ident()], output_fid.get_ident(), split_num, config.HDFS_PREFIX)
        graph.add_node(split_reader_node)

    # Remove the HDFS Cat command as it's not used anymore
    graph.remove_node(hdfs_cat.get_id())

    ## input of next command is output of new merger.
    input_id = next_node.get_standard_inputs()[0]
    new_merger = make_cat_node(output_ids, input_id)
    graph.add_node(new_merger)

    return new_merger


## TODO: There needs to be some state to keep track of open r-split sessions
##       (that either end at r-merge or at r_unwrap before a commutative command).
##
## TODO: At the moment we greedily try to add r-splits if possible, so we need to have a better procedure of deciding whether to put them or not.
##       For example for non-commutative pure commands.

## This function takes a node (id) and parallelizes it
def parallelize_node(curr_id, graph, fileIdGen, fan_out,
                     batch_size, no_cat_split_vanish, r_split_flag, r_split_batch_size):
    assert(False)
    curr = graph.get_node(curr_id)
    new_nodes_for_workset = []

    # TODO: this whole fragment could be moved to the graph after picking a parallelizer
    option_parallelizer_rr = curr.get_option_implemented_round_robin_parallelizer()
    # for now, we use the `r_split_flag` here again:
    if r_split_flag and option_parallelizer_rr is not None:
        parallelizer_rr = option_parallelizer_rr
        aggregator_spec = parallelizer_rr.get_aggregator_spec()
        aggregator_kind = aggregator_spec.get_kind()
        if aggregator_kind == AggregatorKindEnum.CONCATENATE: # is turned into an r_merge
            streaming_inputs = curr.get_streaming_inputs()
            assert(len(streaming_inputs) == 1)
            streaming_input = streaming_inputs[0]
            configuration_inputs = curr.get_configuration_inputs()
            assert(len(configuration_inputs) == 0)
            streaming_outputs = curr.get_output_list()
            assert(len(streaming_outputs) == 1)
            streaming_output = streaming_outputs[0]
            original_cmd_invocation_with_io_vars = curr.cmd_invocation_with_io_vars

            graph.remove_node(curr_id) # remove it here already as as we need to remove edge end points ow. to avoid disconnecting graph to avoid disconnecting graph

            out_split_ids = graph.generate_ephemeral_edges(fileIdGen, fan_out)
            splitter = r_split.make_r_split(streaming_input, out_split_ids, r_split_batch_size)
            graph.set_edge_to(streaming_input, splitter.get_id())
            for out_split_id in out_split_ids:
                graph.set_edge_from(out_split_id, splitter.get_id())
            graph.add_node(splitter)

            in_mapper_ids = out_split_ids
            out_mapper_ids = graph.generate_ephemeral_edges(fileIdGen, fan_out)
            zip_mapper_in_out_ids = zip(in_mapper_ids, out_mapper_ids)

            all_mappers = []
            for (in_id, out_id) in zip_mapper_in_out_ids:
                # BEGIN: these 4 lines could be refactored to be a function in graph such that
                # creating end point of edges and the creation of edges is not decoupled
                mapper_cmd_inv = parallelizer_rr.get_actual_mapper(original_cmd_invocation_with_io_vars, in_id, out_id)
                mapper = DFGNode.make_simple_dfg_node_from_cmd_inv_with_io_vars(mapper_cmd_inv)
                # add r_wrap here:
                mapper_r_wrapped = r_wrap.wrap_node(mapper, graph.edges)
                graph.set_edge_to(in_id, mapper_r_wrapped.get_id())
                graph.set_edge_from(out_id, mapper_r_wrapped.get_id())
                # END
                all_mappers.append(mapper_r_wrapped)
            for new_node in all_mappers:
                graph.add_node(new_node)

            in_aggregator_ids = out_mapper_ids
            out_aggregator_id = streaming_output
            aggregator = r_merge.make_r_merge_node(in_aggregator_ids, out_aggregator_id)
            for in_aggregator_id in in_aggregator_ids:
                graph.set_edge_to(in_aggregator_id, aggregator.get_id())
            graph.set_edge_from(streaming_output, aggregator.get_id())
            all_aggregators = [aggregator]
            ## Add the merge commands in the graph
            for new_node in all_aggregators:
                graph.add_node(new_node)
    elif option_parallelizer_rr is not None: # do consecutive chunks
        # TODO: we do consecutive chunks here but from a rr splitter
        parallelizer_rr = option_parallelizer_rr
        streaming_inputs = curr.get_streaming_inputs()
        assert(len(streaming_inputs) == 1)
        streaming_input = streaming_inputs[0]
        configuration_inputs = curr.get_configuration_inputs()
        assert(len(configuration_inputs) == 0)
        streaming_outputs = curr.get_output_list()
        assert(len(streaming_outputs) == 1)
        streaming_output = streaming_outputs[0]
        original_cmd_invocation_with_io_vars = curr.cmd_invocation_with_io_vars

        graph.remove_node(curr_id) # remove it here already as as we need to remove edge end points ow. to avoid disconnecting graph to avoid disconnecting graph

        out_split_ids = graph.generate_ephemeral_edges(fileIdGen, fan_out)
        splitter = pash_split.make_split_file(streaming_input, out_split_ids)
        graph.set_edge_to(streaming_input, splitter.get_id())
        for out_split_id in out_split_ids:
            graph.set_edge_from(out_split_id, splitter.get_id())
        graph.add_node(splitter)

        in_mapper_ids = out_split_ids
        out_mapper_ids = graph.generate_ephemeral_edges(fileIdGen, fan_out)
        zip_mapper_in_out_ids = zip(in_mapper_ids, out_mapper_ids)

        all_mappers = []
        for (in_id, out_id) in zip_mapper_in_out_ids:
            # BEGIN: these 4 lines could be refactored to be a function in graph such that
            # creating end point of edges and the creation of edges is not decoupled
            mapper_cmd_inv = parallelizer_rr.get_actual_mapper(original_cmd_invocation_with_io_vars, in_id, out_id)
            mapper = DFGNode.make_simple_dfg_node_from_cmd_inv_with_io_vars(mapper_cmd_inv)
            graph.set_edge_to(in_id, mapper.get_id())
            graph.set_edge_from(out_id, mapper.get_id())
            # END
            all_mappers.append(mapper)
        for new_node in all_mappers:
            graph.add_node(new_node)

        in_aggregator_ids = out_mapper_ids
        out_aggregator_id = streaming_output
        aggregator_spec = parallelizer_rr.get_aggregator_spec()
        aggregator_kind = aggregator_spec.get_kind()
        if aggregator_kind == AggregatorKindEnum.CONCATENATE or aggregator_kind == AggregatorKindEnum.CUSTOM_N_ARY:
            aggregator_cmd_inv = parallelizer_rr.get_actual_aggregator(original_cmd_invocation_with_io_vars, in_aggregator_ids, out_aggregator_id)
            aggregator = DFGNode.make_simple_dfg_node_from_cmd_inv_with_io_vars(aggregator_cmd_inv)
            for in_aggregator_id in in_aggregator_ids:
                graph.set_edge_to(in_aggregator_id, aggregator.get_id())
            graph.set_edge_from(streaming_output, aggregator.get_id())
            all_aggregators = [aggregator]
            ## Add the merge commands in the graph
            for new_node in all_aggregators:
                graph.add_node(new_node)
        elif aggregator_kind == AggregatorKindEnum.CUSTOM_2_ARY:
            # TODO: we simplify and assume that every mapper produces a single output for now:
            map_in_aggregator_ids = [[id] for id in in_aggregator_ids]
            graph.create_generic_aggregator_tree(original_cmd_invocation_with_io_vars, parallelizer_rr, map_in_aggregator_ids, out_aggregator_id, fileIdGen)
        else:
            raise Exception("aggregator kind not yet implemented")

    return new_nodes_for_workset

## TODO: Instead of moving a cat after a node, we need to parallelize cat,
##       then remove cat (since it takes a single input to a single output),
##       then parallelize the next node. This will allow us to handle `comm -23 p1 p2`
##
## TODO: A nice interface would be (check/apply transformation)
##
## TODO: This could be a method of IR.
##
## TODO: We need to check if the previous node is a cat or a merge
def check_parallelize_dfg_node(merger_id, node_id, graph, fileIdGen):
    assert(False)

    ## Get merger inputs (cat or r_merge).
    merger_input_edge_ids = graph.get_node_input_ids(merger_id)

    ## If the merger has more than one input, then the next node could be parallelized
    new_nodes = []
    if (len(merger_input_edge_ids) > 1):
        ## If the merger is r-merge, then the next node needs to either be stateless, or commutative parallelizable.
        merger = graph.get_node(merger_id)
        node = graph.get_node(node_id)
        if((isinstance(merger, Cat)
            and node.is_parallelizable())
           or (isinstance(merger, r_merge.RMerge)
               and (node.is_stateless()
                    or node.is_commutative()))):
            new_nodes = parallelize_dfg_node(merger_id, node_id, graph, fileIdGen)

    return new_nodes

def parallelize_dfg_node(old_merger_id, node_id, graph, fileIdGen):
    assert(False)
    node = graph.get_node(node_id)
    assert(node.is_parallelizable())

    ## TODO: Delete this
    ## Get cat inputs and output. Note that there is only one output.
    # old_merger_input_edge_ids = graph.get_node_input_ids(old_merger_id)
    # old_merger_output_edge_ids = graph.get_node_output_ids(old_merger_id)
    # assert(len(old_merger_output_edge_ids) == 1)
    # old_merger_output_edge_id = old_merger_output_edge_ids[0]

    new_nodes = []
    ## We assume that every stateless and pure parallelizable command has one output_file_id for now.
    ##
    ## TODO: Check if this can be lifted.
    node_output_edge_ids = graph.get_node_output_ids(node_id)
    assert(len(node_output_edge_ids) == 1)
    node_output_edge_id = node_output_edge_ids[0]

    ## TODO: Add a commutativity check before actually applying this transformation if the current node is pure parallelizable.
    new_parallel_nodes, map_output_ids = graph.parallelize_node(node_id, fileIdGen)
    new_nodes += new_parallel_nodes

    # log("after duplicate graph nodes:", graph.nodes)
    # log("after duplicate graph edges:", graph.edges)

    ## Make a merge command that joins the results of all the duplicated commands
    ##
    ## TODO: We need to figure out what to do with r_merge when commands are not commutative
    if(node.is_pure_parallelizable()):
        merge_commands, new_edges, final_output_id = create_merge_commands(node,
                                                                           map_output_ids,
                                                                           fileIdGen)
        graph.add_edges(new_edges)

        ## Add the merge commands in the graph
        for merge_command in merge_commands:
            graph.add_node(merge_command)

        ## Replace the previous final_output_id with the previous id
        final_merge_node_id = graph.edges[final_output_id][1]
        final_merge_node = graph.get_node(final_merge_node_id)
        final_merge_node.replace_edge(final_output_id, node_output_edge_id)
        graph.set_edge_from(node_output_edge_id, final_merge_node_id)
        graph.set_edge_from(final_output_id, None)

        ## Only add the final node to the new_nodes
        new_nodes.append(final_merge_node)


    # log("after merge graph nodes:", graph.nodes)
    # log("after merge graph edges:", graph.edges)

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
    return new_nodes

## Creates a merge command for all pure commands that can be
## parallelized using a map and a reduce/merge step
##
## Currently adding an aggregator can be done by adding another branch to this function
##
## TODO: Make that generic to work through annotations
def create_merge_commands(curr, new_output_ids, fileIdGen):
    assert(False)
    if(str(curr.com_name) == "uniq"):
        return create_uniq_merge_commands(curr, new_output_ids, fileIdGen)
    else:
        return create_generic_aggregator_tree(curr, new_output_ids, fileIdGen)

## TODO: These must be generated using some file information
##
## TODO: Find a better place to put these functions
def create_sort_merge_commands(curr, new_output_ids, fileIdGen):
    assert(False)
    output = create_reduce_tree(lambda ids: SortGReduce(curr, ids),
                                new_output_ids, fileIdGen)
    return output

## Instead of creating a tree, we just create a single level reducer for uniq
def create_uniq_merge_commands(curr, new_output_ids, fileIdGen):
    assert(False)
    ## Make an intermediate cat node
    intermediate_fid = fileIdGen.next_ephemeral_file_id()
    intermediate_id = intermediate_fid.get_ident()
    new_cat = make_cat_node(flatten_list(new_output_ids), intermediate_id)

    ## Make the new uniq output
    new_out_fid = fileIdGen.next_ephemeral_file_id()
    new_out_id = new_out_fid.get_ident()

    ## TODO: Pass the options of `curr` correctly

    ## Make the uniq merge node
    uniq_com_name = Arg(string_to_argument("uniq"))
    com_category = "pure_parallelizable"
    node = DFGNode([intermediate_id],
                   [new_out_id],
                   uniq_com_name,
                   com_category)

    return ([new_cat, node], [intermediate_fid, new_out_fid], new_out_id)


## This functions adds an eager on a given edge.
def add_eager(eager_input_id, graph, fileIdGen, intermediateFileIdGen, use_dgsh_tee):
    new_fid = fileIdGen.next_ephemeral_file_id()
    new_id = new_fid.get_ident()

    if use_dgsh_tee:
        ## TODO: seperate to better use dgsh-tee params and maybe deprecate eager
        eager_node = dgsh_tee.make_dgsh_tee_node(eager_input_id, new_id)
    else:
        ## TODO: Remove the line below if eager creates its intermediate file
        ##       on its own.
        # TODO: find a better solution to make unique numbers, currently: set to max-value + 1
        intermediateFileIdGen.bump_counter_to_value_of(fileIdGen)
        intermediate_fid = intermediateFileIdGen.next_temporary_file_id()
        # TODO: this edge will never have to since eager is set to output even though it reads from it
        graph.add_edge(intermediate_fid)
        fileIdGen.bump_counter_to_value_of(intermediateFileIdGen)

        eager_exec_path = '{}/{}'.format(config.PASH_TOP, runtime_config['eager_executable_path'])

        eager_node = make_eager_node(eager_input_id, new_id, intermediate_fid, eager_exec_path)

    ## Add the edges and the nodes to the graph
    graph.add_edge(new_fid)

    ## Modify the next node inputs to be the new inputs
    next_node_id = graph.edges[eager_input_id][2]
    if(not next_node_id is None):
        next_node = graph.get_node(next_node_id)
        next_node.replace_edge(eager_input_id, new_id)
        graph.set_edge_to(new_id, next_node_id)

    graph.add_node(eager_node)


## This function adds eager nodes wherever the width of graph is
## becoming smaller.
def add_eager_nodes(graph, use_dgsh_tee):
    source_node_ids = graph.source_nodes()

    ## Generate a fileIdGen that doesnt clash with graph fids.
    fileIdGen = graph.get_file_id_gen()
    intermediateFileIdGen = FileIdGen(0, runtime_config['eager_intermediate_prefix'])

    ## Get the next nodes
    workset = [node for source_node_id in source_node_ids for node in graph.get_next_nodes(source_node_id)]
    visited = set()
    while (len(workset) > 0):
        curr_id = workset.pop(0)
        curr = graph.get_node(curr_id)
        if (not curr_id in visited):
            visited.add(curr_id)
            next_node_ids = graph.get_next_nodes(curr_id)
            workset += next_node_ids

            ## TODO: Make sure that we don't add duplicate eager nodes

            ## Add eager nodes if the node has more than one input
            curr_input_ids = graph.get_node_input_ids(curr_id)
            if (len(curr_input_ids) > 1):
                ## TODO: If we know that a command reads its inputs in a list,
                ##       then we might not need to put an eager on its first input.
                ## Note: This cannot be done for `sort -m` so we need to know in the
                ##       annotations whether input consumption is in order or not.

                for curr_input_id in curr_input_ids:
                    _fid, from_node, to_node = graph.edges[curr_input_id]
                    assert(to_node == curr_id)
                    ## If the edge is an input edge, then we don't want to put eager.
                    if(not from_node is None):
                        add_eager(curr_input_id, graph, fileIdGen, intermediateFileIdGen, use_dgsh_tee)

            if(isinstance(curr, Split)):
                eager_input_ids = curr.get_output_list()[:-1]
                for edge_id in eager_input_ids:
                    add_eager(edge_id, graph, fileIdGen, intermediateFileIdGen, use_dgsh_tee)

            ## Add an eager after r_unwrap            
            if(isinstance(curr, r_unwrap.RUnwrap)):
                eager_input_id = curr.get_output_list()[0]
                add_eager(eager_input_id, graph, fileIdGen, intermediateFileIdGen, use_dgsh_tee)

            ## Add an eager after r_split
            if(isinstance(curr, r_split.RSplit)):
                eager_input_ids = curr.get_output_list()
                for edge_id in eager_input_ids:
                    add_eager(edge_id, graph, fileIdGen, intermediateFileIdGen, use_dgsh_tee)

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
