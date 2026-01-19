from pash_annotations.datatypes.AccessKind import (
    AccessKind,
    make_stream_output,
    make_stream_input,
    make_other_output,
)
from pash_annotations.datatypes.CommandInvocationWithIOVars import (
    CommandInvocationWithIOVars,
)

from definitions.ir.dfg_node import *


class Eager(DFGNode):
    def __init__(self, cmd_invocation_with_io_vars, com_redirs=None, com_assignments=None):
        com_redirs = [] if com_redirs is None else com_redirs
        com_assignments = [] if com_assignments is None else com_assignments

        super().__init__(
            cmd_invocation_with_io_vars,
            com_redirs=com_redirs,
            com_assignments=com_assignments,
        )


def make_eager_node(input_id, output_id, intermediate_file_id, eager_exec_path):
    eager_name = eager_exec_path
    intermediate_file_id_id = intermediate_file_id.get_ident()
    operand_list = [input_id, output_id, intermediate_file_id_id]
    access_map = {
        output_id: make_stream_output(),
        input_id: make_stream_input(),
        intermediate_file_id_id: make_other_output(),
    }
    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=eager_name,
        flag_option_list=[],
        operand_list=operand_list,
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=None,
        access_map=access_map,
    )
    return Eager(cmd_inv_with_io_vars)
