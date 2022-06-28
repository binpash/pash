from datatypes_new.AccessKind import AccessKind
from datatypes_new.CommandInvocationWithIOVars import CommandInvocationWithIOVars

from definitions.ir.dfg_node import *
from ir_utils import *

class RUnwrap(DFGNode):
    def __init__(self,
                 cmd_invocation_with_io_vars,
                 com_redirs=[],
                 com_assignments=[],
                 parallelizer_list=None,
                 cmd_related_properties=None):
        # TODO []: default
        super().__init__(cmd_invocation_with_io_vars,
                         com_redirs=com_redirs,
                         com_assignments=com_assignments,
                         parallelizer_list=parallelizer_list,
                         cmd_related_properties=cmd_related_properties)

def make_unwrap_node(inputs, output):
    assert(len(inputs) == 1)
    input_id = inputs[0]
    access_map = {input_id: AccessKind.make_stream_input(), output: AccessKind.make_stream_output()}
    r_unwrap_bin = os.path.join(config.PASH_TOP, config.config['runtime']['r_unwrap_binary'])
    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=r_unwrap_bin,
        flag_option_list=[],
        operand_list=[],
        implicit_use_of_streaming_input=input_id,
        implicit_use_of_streaming_output=output,
        access_map=access_map)
    return RUnwrap(cmd_inv_with_io_vars)
