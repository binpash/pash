# TODO: this file can properly be deleted

# imports from annotation framework
# for typing
# for use
from pash_annotations.annotation_generation.datatypes.parallelizability.Mapper import Mapper

from definitions.ir.dfg_node import DFGNode
from annotations_utils.util_cmd_invocations import get_command_invocation_prefix_from_dfg_node
from util import log

def get_actual_mapper_from_node(node, parallelizer) -> Mapper:
    assert(False)
    cmd_inv_pref = get_command_invocation_prefix_from_dfg_node(node)
    return parallelizer.get_actual_mapper(cmd_inv_pref)

def get_mapper_as_dfg_node_from_node(node, parallelizer, inputs, outputs) -> DFGNode:
    assert(False)
    mapper = get_actual_mapper_from_node(node, parallelizer)
    log(f'mapper for cmd_name: {node.com_name}')
    log(f'here mapper: {mapper}')
    return DFGNode(inputs=inputs,
                outputs=outputs,
                com_name=mapper.cmd_name,
                # com_options=node.com_options,
                com_redirs=node.com_redirs,
                com_assignments=node.com_assignments,
                flag_option_list=mapper.flag_option_list,
                positional_config_list=mapper.positional_config_list,
                positional_input_list=None,     # TODO: somehow from inputs, future shift
                positional_output_list=None    # TODO: somehow from outputs, future shift
            # TODO:
            # implicit_use_of_stdin = False,
            # implicit_use_of_stdout = False,
            # omitted for now since we do not consider nested parallelization
            # parallelizer_list = None,
            # cmd_related_properties = None,
    )

## MOVED from dfg_node
## Get the file names of the outputs of the map commands. This
## differs if the command is stateless, pure that can be
## written as a map and a reduce, and a pure that can be
## written as a generalized map and reduce.
# BEGIN ANNO
# OLD
# def get_map_output_files(node, input_edge_ids, fileIdGen):
# NEW
def get_map_output_files(node, input_edge_ids, fileIdGen, parallelizer):
    assert(False)
    assert (node.is_parallelizable())
    # TODO ANNO: How to substitute? @KK
    if (node.com_category == "stateless"):
        map_output_fids = [fileIdGen.next_ephemeral_file_id() for in_fid in input_edge_ids]
    elif (node.is_pure_parallelizable()):
        # BEGIN ANNO
        # OLD
        # map_output_fids = node.pure_get_map_output_files(input_edge_ids, fileIdGen)
        # NEW
        map_output_fids = pure_get_map_output_files(node, input_edge_ids, fileIdGen, parallelizer)
        # END ANNO
    else:
        log("Unreachable code reached :(")
        assert (False)
        ## This should be unreachable

    return map_output_fids

## TODO: Fix this somewhere in the annotations and not in the code
# BEGIN ANNO
# OLD
# def pure_get_map_output_files(node, input_edge_ids, fileIdGen):
# NEW
def pure_get_map_output_files(node, input_edge_ids, fileIdGen, parallelizer):
    assert(False)
    assert (node.is_pure_parallelizable())
    # BEGIN ANNO
    # OLD
    ## The number of the mapper outputs defaults to 1
    # if(node.com_mapper is None):
    #     number_outputs = 1
    # else:
    #     number_outputs = node.com_mapper.num_outputs
    # NEW
    # TODO: which parallelizer did we choose?
    actual_mapper = get_actual_mapper_from_node(node, parallelizer)
    number_outputs = actual_mapper.num_outputs  # defaults to 1 in class Mapper
    # END ANNO

    new_output_fids = [[fileIdGen.next_ephemeral_file_id() for i in range(number_outputs)]
                       for in_fid in input_edge_ids]
    return new_output_fids


