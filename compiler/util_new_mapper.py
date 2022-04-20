# imports from annotation framework
import sys
sys.path.insert(1, "/home/felix/git-repos/MIT/annotations")
# for typing
# for use
from annotation_generation_new.datatypes.parallelizability.Mapper import Mapper

from definitions.ir.dfg_node import DFGNode
from util_new_cmd_invocations import get_command_invocation_prefix_from_dfg_node
from util import log

def get_mapper_as_dfg_node_from_node(node, parallelizer, inputs, outputs) -> DFGNode:
    cmd_inv_pref = get_command_invocation_prefix_from_dfg_node(node)
    mapper = parallelizer.get_actual_mapper(cmd_inv_pref)
    log(f'node.com_options: {node.com_options}')
    return DFGNode(inputs=inputs,
                outputs=outputs,
                com_name=mapper.cmd_name,
                com_options=node.com_options,
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





