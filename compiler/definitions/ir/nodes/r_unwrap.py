from pash_annotations.datatypes.AccessKind import make_stream_input, make_stream_output
from pash_annotations.datatypes.CommandInvocationWithIOVars import (
    CommandInvocationWithIOVars,
)

from definitions.ir.dfg_node import *


class RUnwrap(DFGNode):
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
            cmd_invocation_with_io_vars,
            com_redirs=com_redirs,
            com_assignments=com_assignments,
            parallelizer_list=parallelizer_list,
            cmd_related_properties=cmd_related_properties,
        )


def make_unwrap_node(inputs, output):
    assert len(inputs) == 1
    input_id = inputs[0]
    access_map = {input_id: make_stream_input(), output: make_stream_output()}
    r_unwrap_bin = os.path.join(
        config.PASH_TOP, config.config["runtime"]["r_unwrap_binary"]
    )
    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=r_unwrap_bin,
        flag_option_list=[],
        operand_list=[],
        implicit_use_of_streaming_input=input_id,
        implicit_use_of_streaming_output=output,
        access_map=access_map,
    )
    return RUnwrap(cmd_inv_with_io_vars)
