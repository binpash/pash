from definitions.ir.dfg_node import *

# from definitions.ir.nodes.arg import Arg
from annotations_utils.util_cmd_invocations import (
    get_command_invocation_prefix_from_dfg_node,
)


## This class corresponds to a generic n-ary aggregator
##
## TODO: Do we need to do anything special for binary aggregators?
class MapperAggregatorNode(DFGNode):
    def __init__(
        self,
        old_node,
        input_ids,
        output_ids,
        name_string,
        new_options,
        flag_option_list,
    ):
        ## The name of the aggregator command
        name = Arg.string_to_arg(name_string)

        ## TODO: The category should also be acquired through annotations (and maybe should be asserted to be at most pure)
        com_category = "pure"

        ## TODO: Not sure if redirections need to be copied to new function.
        com_redirs = [redir.to_ast() for redir in old_node.com_redirs]
        super().__init__(
            input_ids,
            output_ids,
            name,
            com_category,
            com_options=new_options,  # changed that all are already in there and not appended
            flag_option_list=flag_option_list,
            com_redirs=com_redirs,
            com_assignments=old_node.com_assignments,
        )


class AggregatorNode(MapperAggregatorNode):
    def __init__(self, old_node, input_ids, output_ids):
        used_parallelizer = old_node.get_used_parallelizer()
        cmd_inv_pref = get_command_invocation_prefix_from_dfg_node(old_node)
        used_aggregator = used_parallelizer.get_actual_aggregator(cmd_inv_pref)
        log(f"used_agg: {used_aggregator}")
        log(f"old_node: {old_node}")

        ## Check if an aggregator can be instantiated from the node
        if used_aggregator is None:
            log(
                "Error: Node:",
                old_node,
                "does not contain information to instantiate an aggregator!",
            )
            raise Exception("No information to instantiate aggregator")

        ## The name of the aggregator command
        agg_name_string = used_aggregator.cmd_name
        all_options_incl_new = [
            Arg.string_to_arg(el.get_name())
            for el in used_aggregator.flag_option_list
            + used_aggregator.positional_config_list
        ]
        # TODO: zip is nicer
        all_options_incl_new_right_format = [
            (i, all_options_incl_new[i]) for i in range(len(all_options_incl_new))
        ]

        super().__init__(
            old_node,
            input_ids,
            output_ids,
            agg_name_string,
            all_options_incl_new_right_format,
            flag_option_list=used_aggregator.flag_option_list,
        )

        log("Generic Aggregator Created:", self)
