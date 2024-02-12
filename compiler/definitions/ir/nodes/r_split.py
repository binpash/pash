import os

from pash_annotations.datatypes.AccessKind import (
    AccessKind,
    make_stream_input,
    make_stream_output,
)
from pash_annotations.datatypes.BasicDatatypes import Operand, Flag
from pash_annotations.datatypes.CommandInvocationWithIOVars import (
    CommandInvocationWithIOVars,
)

import config

from definitions.ir.dfg_node import *
from definitions.ir.file_id import *
from shell_ast.ast_util import string_to_argument


class RSplit(DFGNode):
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

    def add_r_flag(self):
        self.cmd_invocation_with_io_vars.flag_option_list.append(Flag("-r"))


def make_r_split(input_id, out_ids, r_split_batch_size):
    r_split_bin = os.path.join(
        config.PASH_TOP, config.config["runtime"]["r_split_binary"]
    )
    operand_list = [input_id, Operand(Arg.string_to_arg(str(r_split_batch_size)))]
    operand_list.extend(out_ids)
    access_map = {output_id: make_stream_output() for output_id in out_ids}
    access_map[input_id] = make_stream_input()
    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=r_split_bin,
        flag_option_list=[],
        operand_list=operand_list,
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=None,
        access_map=access_map,
    )
    return RSplit(cmd_inv_with_io_vars)


def make_r_split_with_unwrap_flag(input_id, out_ids, r_split_batch_size):
    standard_r_split = make_r_split(input_id, out_ids, r_split_batch_size)
    standard_r_split.add_r_flag()
    return standard_r_split
