from datatypes_new.AccessKind import AccessKind
from datatypes_new.CommandInvocationWithIOVars import CommandInvocationWithIOVars

from definitions.ir.file_id import *
from definitions.ir.dfg_node import *
from ir_utils import string_to_argument

import config
import os

class Split(DFGNode):
    def __init__(self,
                 cmd_invocation_with_io_vars,
                 com_redirs=[],
                 com_assignments=[],
                 parallelizer_list=None,
                 cmd_related_properties=None):
        # TODO []: default arguments!
        super().__init__(cmd_invocation_with_io_vars=cmd_invocation_with_io_vars,
                         com_redirs=com_redirs,
                         com_assignments=com_assignments,
                         parallelizer_list=parallelizer_list,
                         cmd_related_properties=cmd_related_properties)

def make_split_file(input_id, out_ids):
    auto_split_bin = os.path.join(config.PASH_TOP, config.config['runtime']['auto_split_binary'])
    operand_list = [input_id]
    operand_list.extend(out_ids)
    access_map = {output_id: AccessKind.make_stream_output() for output_id in out_ids}
    access_map[input_id] = AccessKind.make_stream_input()
    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=auto_split_bin,
        flag_option_list=[],
        operand_list=operand_list,
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=None,
        access_map=access_map)
    return Split(cmd_inv_with_io_vars)
