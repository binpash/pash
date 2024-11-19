import sys
import pickle
import traceback
from datetime import datetime

from sh_expand import env_vars_util
from sh_expand.expand import ExpansionError 

import config
from ir import *
from ast_to_ir import compile_asts
from ir_to_ast import to_shell
from pash_graphviz import maybe_generate_graphviz
from util import *
from custom_error import *

from definitions.ir.aggregator_node import *

from definitions.ir.nodes.eager import *
from definitions.ir.nodes.pash_split import *

import definitions.ir.nodes.r_split as r_split
import definitions.ir.nodes.r_unwrap as r_unwrap
import definitions.ir.nodes.dgsh_tee as dgsh_tee
import definitions.ir.nodes.dfs_split_reader as dfs_split_reader

# Distirbuted Exec
import dspash.hdfs_utils as hdfs_utils

from cli import CompilerParser

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
    config.set_config_globals_from_pash_args(args)

    ## Load the configuration
    if not config.config:
        config.load_config(args.config_path)

    runtime_config = config.config["distr_planner"]

    ## Read any shell variables files if present
    vars_dict = env_vars_util.read_vars_file(args.var_file, config.BASH_VERSION)
    config.set_vars_file(args.var_file, vars_dict)

    log("Input:", args.input_ir, "Compiled file:", args.compiled_script_file)

    ## Call the main procedure
    compiler_config = CompilerConfig(args.width)
    ast_or_ir = compile_optimize_output_script(
        args.input_ir, args.compiled_script_file, args, compiler_config
    )
    maybe_generate_graphviz(ast_or_ir, args)


def parse_args():
    parser = CompilerParser()
    args, _ = parser.parse_known_args()
    return args


## TODO: Add more fields from args in this
class CompilerConfig:
    def __init__(self, width):
        self.width = width

    def __repr__(self):
        return f"CompilerConfig(Width:{self.width})"


def compile_ir(ir_filename, compiled_script_file, args, compiler_config):
    """
    Return IR object for compilation success. None otherwise.
    """
    ret = None
    try:
        ret = compile_optimize_output_script(
            ir_filename, compiled_script_file, args, compiler_config
        )
    except ExpansionError as e: 
        log("WARNING: Exception caught because some region(s) are not expandable and therefore unparallelizable:", e) 
        raise NotAllRegionParallelizableError()
    except UnparallelizableError as e: 
        log("WARNING: Exception caught because some region(s) are unparallelizable:", e) 
        raise NotAllRegionParallelizableError()
        # log(traceback.format_exc()) # uncomment for exact trace report (PaSh user should see informative messages for unparellizable regions) 
    except (AdjLineNotImplementedError, NotImplementedError) as e: 
        log("WARNING: Exception caught because some part is not implemented:", e)
        log(traceback.format_exc())
    except Exception as e:
        log("WARNING: Exception caught:", e)
        log(traceback.format_exc())

    return ret


def compile_optimize_output_script(
    ir_filename, compiled_script_file, args, compiler_config
):
    global runtime_config

    ret = None

    ## Load the df_region from a file
    candidate_df_region = load_df_region(ir_filename)

    ## Compile it
    optimized_ast_or_ir = compile_optimize_df_region(
        candidate_df_region, args, compiler_config
    )

    ## Call the backend that executes the optimized dataflow graph
    ## TODO: Should never be the case for now. This is obsolete.
    assert not runtime_config["distr_backend"]

    ## If the candidate DF region was indeed a DF region then we have an IR
    ## which should be translated to a parallel script.
    if isinstance(optimized_ast_or_ir, IR):
        if args.distributed_exec:
            ir_filename = ptempfile()
            script_to_execute = (
                f"$PASH_TOP/compiler/dspash/remote_exec_graph.sh {ir_filename}\n"
            )
            ## This might not be needed anymore (since the output script is output anyway)
            ## TODO: This is probably useless, remove
            maybe_log_optimized_script(script_to_execute, args)

            with open(ir_filename, "wb") as f:
                obj = (optimized_ast_or_ir, config.config["shell_variables"])
                pickle.dump(obj, f)
        else:
            script_to_execute = to_shell(optimized_ast_or_ir, args)

        log("Optimized script saved in:", compiled_script_file)
        with open(compiled_script_file, "w") as f:
            f.write(script_to_execute)

        ret = optimized_ast_or_ir
    else:
        raise UnparallelizableError("Script failed to compile!")

    return ret


def load_df_region(ir_filename):
    log("Retrieving candidate DF region: {} ... ".format(ir_filename), end="")
    with open(ir_filename, "rb") as ir_file:
        candidate_df_region = pickle.load(ir_file)
    log("Done!")
    return candidate_df_region


def compile_optimize_df_region(df_region, args, compiler_config):
    ## Compile the candidate DF regions
    compilation_start_time = datetime.now()
    asts_and_irs = compile_candidate_df_region(df_region, config.config)
    compilation_end_time = datetime.now()
    print_time_delta("Compilation", compilation_start_time, compilation_end_time)

    ## Optimize all the IRs that can be optimized
    if args.no_optimize:
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
    assert len(optimized_asts_and_irs) == 1
    optimized_ast_or_ir = optimized_asts_and_irs[0]

    return optimized_ast_or_ir


def maybe_log_optimized_script(script_to_execute, args):
    ## TODO: Merge this write with the one below. Maybe even move this logic in `pash_runtime.sh`
    ## Output the optimized shell script for inspection
    if args.output_optimized:
        output_script_path = runtime_config["optimized_script_filename"]
        with open(output_script_path, "w") as output_script_file:
            log("Optimized script:")
            log(script_to_execute)
            output_script_file.write(script_to_execute)


def compile_candidate_df_region(candidate_df_region, config):
    ## This is for the files in the IR
    fileIdGen = FileIdGen()

    ## If the candidate DF region is not from the top level then
    ## it won't be a list and thus we need to make it into a list to compile it.
    if not isinstance(candidate_df_region, list):
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
        if isinstance(ast_or_ir, IR):
            ## Assert that the graph that was returned from compilation is valid
            assert ast_or_ir.valid()

            # log(ir_node)
            # with cProfile.Profile() as pr:
            distributed_graph = choose_and_apply_parallelizing_transformations(
                ast_or_ir,
                compiler_config.width,
                runtime_config["batch_size"],
                args.r_split_batch_size,
            )
            # pr.print_stats()

            # Eagers are added in remote notes when using distributed exec
            if not args.no_eager and not args.distributed_exec:
                eager_distributed_graph = add_eager_nodes(distributed_graph)
            else:
                eager_distributed_graph = distributed_graph

            ## Assert that the graph stayed valid after all transformations
            assert eager_distributed_graph.valid()

            ## Print statistics of output nodes
            print_graph_statistics(eager_distributed_graph)

            optimized_asts_and_irs.append(eager_distributed_graph)
        else:
            optimized_asts_and_irs.append(ast_or_ir)

    optimization_end_time = datetime.now()
    print_time_delta("Optimization", optimization_start_time, optimization_end_time)

    return optimized_asts_and_irs


def print_graph_statistics(graph):
    total_nodes = graph.nodes
    eager_nodes = [node for node in total_nodes.values() if isinstance(node, Eager)]
    log("Total nodes after optimization:", len(total_nodes))
    log(" -- out of which:")
    log("Eager nodes:", len(eager_nodes))


def choose_and_apply_parallelizing_transformations(
    graph, fan_out, batch_size, r_split_batch_size
):
    parallelizer_map = choose_parallelizing_transformations(graph)
    apply_parallelizing_transformations(
        graph, parallelizer_map, fan_out, batch_size, r_split_batch_size
    )
    return graph


def choose_parallelizing_transformations(graph):  # shall return map
    source_node_ids = graph.source_nodes()
    parallelizer_map = {}
    workset = source_node_ids
    visited = set()
    # We apply a modified BFS such that we ensure that we know which parallelizer was chosen for all previous nodes
    # and assume that the decision for any subsequent node will exploit any potential synergy effects
    while len(workset) > 0:
        curr_id = workset.pop(0)
        assert isinstance(curr_id, int)
        all_previous_nodes_visited = all(
            prev in visited for prev in graph.get_previous_nodes(curr_id)
        )
        if not all_previous_nodes_visited:
            workset.append(curr_id)
        elif not curr_id in visited:
            next_node_ids = graph.get_next_nodes(curr_id)
            workset += next_node_ids
            parallelizer_map[curr_id] = choose_parallelizing_transformation(
                curr_id, graph
            )
            visited.add(curr_id)
    return parallelizer_map


## This currently chooses the best parallelization based on priority:
## 1. The round robin
## 2. The round robin after having performed unwrap (not sure why this is the second priority)
## 3. The consecutive chunks
##
## TODO: In the future, we could develop more complex strategies
def choose_parallelizing_transformation(curr_id, graph):  # shall return map entry
    curr = graph.get_node(curr_id)
    list_all_parallelizers_in_priority = [
        curr.get_option_implemented_round_robin_parallelizer(),
        curr.get_option_implemented_round_robin_with_unwrap_parallelizer(),
        curr.get_option_implemented_consecutive_chunks_parallelizer(),
    ]
    return next(
        (item for item in list_all_parallelizers_in_priority if item is not None), None
    )


def apply_parallelizing_transformations(
    graph, parallelizer_map, fan_out, batch_size, r_split_batch_size
):
    fileIdGen = graph.get_file_id_gen()
    node_id_non_none_parallelizer_list = [
        (node_id, parallelizer)
        for (node_id, parallelizer) in parallelizer_map.items()
        if parallelizer is not None
    ]
    for node_id, parallelizer in node_id_non_none_parallelizer_list:
        graph.apply_parallelization_to_node(
            node_id, parallelizer, fileIdGen, fan_out, r_split_batch_size
        )


def split_hdfs_cat_input(hdfs_cat, next_node, graph, fileIdGen):
    """
    Replaces hdfs cat with a cat per block, each cat uses has an HDFSResource input fid
    Returns: A normal Cat that merges the blocks (will be removed when parallizing next_node)
    """
    assert isinstance(hdfs_cat, HDFSCat)

    ## At the moment this only works for nodes that have one standard input.
    if len(next_node.get_standard_inputs()) != 1:
        return

    hdfscat_input_id = hdfs_cat.get_standard_inputs()[0]
    hdfs_fid = graph.get_edge_fid(hdfscat_input_id)
    hdfs_filepath = str(hdfs_fid.get_resource())
    output_ids = []

    # Create a cat command per file block
    file_config = hdfs_utils.get_file_config(hdfs_filepath)
    dummy_config_path = ptempfile()  # Dummy config file, should be updated by workers
    for split_num, block in enumerate(file_config.blocks):
        resource = DFSSplitResource(
            file_config.dumps(), dummy_config_path, split_num, block.hosts
        )
        block_fid = fileIdGen.next_file_id()
        block_fid.set_resource(resource)
        graph.add_edge(block_fid)

        output_fid = fileIdGen.next_file_id()
        output_fid.make_ephemeral()
        output_ids.append(output_fid.get_ident())
        graph.add_edge(output_fid)

        split_reader_node = dfs_split_reader.make_dfs_split_reader_node(
            [block_fid.get_ident()],
            output_fid.get_ident(),
            split_num,
            config.HDFS_PREFIX,
        )
        graph.add_node(split_reader_node)

    # Remove the HDFS Cat command as it's not used anymore
    graph.remove_node(hdfs_cat.get_id())

    ## input of next command is output of new merger.
    input_id = next_node.get_standard_inputs()[0]
    new_merger = make_cat_node(output_ids, input_id)
    graph.add_node(new_merger)

    return new_merger


## This functions adds an eager on a given edge.
def add_eager(eager_input_id, graph, fileIdGen):
    new_fid = fileIdGen.next_ephemeral_file_id()
    new_id = new_fid.get_ident()

    ## TODO: seperate to better use dgsh-tee params and maybe deprecate eager
    eager_node = dgsh_tee.make_dgsh_tee_node(eager_input_id, new_id)

    ## Add the edges and the nodes to the graph
    graph.add_edge(new_fid)

    ## Modify the next node inputs to be the new inputs
    next_node_id = graph.edges[eager_input_id][2]
    if not next_node_id is None:
        next_node = graph.get_node(next_node_id)
        next_node.replace_edge(eager_input_id, new_id)
        graph.set_edge_to(new_id, next_node_id)

    graph.add_node(eager_node)


## This function adds eager nodes wherever the width of graph is
## becoming smaller.
def add_eager_nodes(graph):
    source_node_ids = graph.source_nodes()

    ## Generate a fileIdGen that doesnt clash with graph fids.
    fileIdGen = graph.get_file_id_gen()

    ## Get the next nodes
    workset = [
        node
        for source_node_id in source_node_ids
        for node in graph.get_next_nodes(source_node_id)
    ]
    visited = set()
    while len(workset) > 0:
        curr_id = workset.pop(0)
        curr = graph.get_node(curr_id)
        if not curr_id in visited:
            visited.add(curr_id)
            next_node_ids = graph.get_next_nodes(curr_id)
            workset += next_node_ids

            ## TODO: Make sure that we don't add duplicate eager nodes

            ## Add eager nodes if the node has more than one input
            curr_input_ids = graph.get_node_input_ids(curr_id)
            if len(curr_input_ids) > 1:
                ## TODO: If we know that a command reads its inputs in a list,
                ##       then we might not need to put an eager on its first input.
                ## Note: This cannot be done for `sort -m` so we need to know in the
                ##       annotations whether input consumption is in order or not.

                for curr_input_id in curr_input_ids:
                    _fid, from_node, to_node = graph.edges[curr_input_id]
                    assert to_node == curr_id
                    ## If the edge is an input edge, then we don't want to put eager.
                    if not from_node is None:
                        add_eager(curr_input_id, graph, fileIdGen)

            if isinstance(curr, Split):
                eager_input_ids = curr.get_output_list()[:-1]
                for edge_id in eager_input_ids:
                    add_eager(edge_id, graph, fileIdGen)

            ## Add an eager after r_unwrap
            if isinstance(curr, r_unwrap.RUnwrap):
                eager_input_id = curr.get_output_list()[0]
                add_eager(eager_input_id, graph, fileIdGen)

            ## Add an eager after r_split
            if isinstance(curr, r_split.RSplit):
                eager_input_ids = curr.get_output_list()
                for edge_id in eager_input_ids:
                    add_eager(edge_id, graph, fileIdGen)

    return graph


if __name__ == "__main__":
    main()
