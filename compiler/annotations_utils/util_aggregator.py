# TODO: this file can properly be deleted

from definitions.ir.dfg_node import DFGNode
from definitions.ir.nodes.cat import Cat
from annotations_utils.util_cmd_invocations import get_command_invocation_prefix_from_dfg_node
from util import log
from ir_utils import string_to_argument
from definitions.ir.arg import Arg

def get_aggregator_as_dfg_node_from_node(node, parallelizer, inputs, outputs) -> DFGNode:
    assert(False)
    cmd_inv_pref = get_command_invocation_prefix_from_dfg_node(node)
    log(f'cmdinvpref for agg: {cmd_inv_pref}')
    aggregator = parallelizer.get_actual_aggregator(cmd_inv_pref)
    log(f'here agg: {aggregator}')
    # TODO: this could be simplified once we use the new attributes
    if aggregator.cmd_name == 'cat':
        return Cat(inputs=inputs,
                    outputs=outputs,
                    com_name=Arg(string_to_argument(aggregator.cmd_name)),
                    com_options=[], # empty and not taking over from other one
                    com_category="stateless",
                    com_redirs=node.com_redirs,
                    com_assignments=node.com_assignments,
                    flag_option_list=aggregator.flag_option_list,
                    positional_config_list=aggregator.positional_config_list,
                    positional_input_list=None,     # TODO: somehow from inputs, future shift
                    positional_output_list=None    # TODO: somehow from outputs, future shift
                # TODO:
                # implicit_use_of_stdin = False,
                # implicit_use_of_stdout = False,
                # omitted for now since we do not consider nested parallelization
                # parallelizer_list = None,
                # cmd_related_properties = None,
        )
    else:
        log(f'agg_com_name: {aggregator.cmd_name}')
        log(f'agg_flag_option_list: {aggregator.flag_option_list}')
        return DFGNode(inputs=inputs,
                       outputs=outputs,
                       com_name=Arg(string_to_argument(aggregator.cmd_name)),
                       com_options=node.com_options,
                       com_redirs=node.com_redirs,
                       com_assignments=node.com_assignments,
                       flag_option_list=aggregator.flag_option_list,
                       positional_config_list=aggregator.positional_config_list,
                       positional_input_list=None,     # TODO: somehow from inputs, future shift
                       positional_output_list=None    # TODO: somehow from outputs, future shift
                       # TODO:
                       # implicit_use_of_stdin = False,
                       # implicit_use_of_stdout = False,
                       # omitted for now since we do not consider nested parallelization
                       # parallelizer_list = None,
                       # cmd_related_properties = None,
                       )
