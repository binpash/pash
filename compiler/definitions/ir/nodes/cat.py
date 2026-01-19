from pash_annotations.datatypes.CommandInvocationWithIOVars import (
    CommandInvocationWithIOVars,
)
from definitions.ir.dfg_node import DFGNode


def make_cat_node(inputs, output):
    cmd_inv_cat = CommandInvocationWithIOVars.make_cat_command_invocation_with_io_vars(
        inputs, output
    )
    return DFGNode.make_simple_dfg_node_from_cmd_inv_with_io_vars(cmd_inv_cat)
