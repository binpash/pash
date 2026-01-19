from pash_annotations.datatypes.AccessKind import make_stream_input, make_stream_output
from pash_annotations.datatypes.CommandInvocationWithIOVars import (
    CommandInvocationWithIOVars,
)

from definitions.ir.dfg_node import *


class RMerge(DFGNode):
    def __init__(
        self,
        cmd_invocation_with_io_vars,
        com_redirs=None,
        com_assignments=None,
        parallelizer_list=None,
        cmd_related_properties=None,
    ):
        com_redirs = [] if com_redirs is None else com_redirs
        com_assignments = [] if com_assignments is None else com_assignments
        super().__init__(
            cmd_invocation_with_io_vars=cmd_invocation_with_io_vars,
            com_redirs=com_redirs,
            com_assignments=com_assignments,
            parallelizer_list=parallelizer_list,
            cmd_related_properties=cmd_related_properties,
        )


def make_r_merge_node(inputs, output):
    r_merge_bin = os.path.join(
        config.PASH_TOP, config.config["runtime"]["r_merge_binary"]
    )
    # TODO: assume that the inputs and output is provided as operands
    access_map = {input_id: make_stream_input() for input_id in inputs}
    access_map[output] = make_stream_output()
    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=r_merge_bin,
        flag_option_list=[],
        operand_list=inputs,
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=output,
        access_map=access_map,
    )
    return RMerge(cmd_inv_with_io_vars)
