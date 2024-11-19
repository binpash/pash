from pash_annotations.datatypes.AccessKind import make_stream_output, make_stream_input
from pash_annotations.datatypes.BasicDatatypes import ArgStringType
from pash_annotations.datatypes.CommandInvocationWithIOVars import (
    CommandInvocationWithIOVars,
)

from annotations_utils.util_cmd_invocations import (
    to_arg_from_cmd_inv_with_io_vars_without_streaming_inputs_or_outputs_for_wrapping,
)
from definitions.ir.dfg_node import *
from shell_ast.ast_util import *


class RWrap(DFGNode):
    def __init__(
        self,
        cmd_invocation_with_io_vars,
        com_redirs=None,
        com_assignments=None,
        parallelizer_list=None,
        cmd_related_properties=None,
        wrapped_node_name=None,
    ):
        com_redirs = [] if com_redirs is None else com_redirs
        com_assignments = [] if com_assignments is None else com_assignments
        self.wrapped_node_name = wrapped_node_name
        super().__init__(
            cmd_invocation_with_io_vars,
            com_redirs=com_redirs,
            com_assignments=com_assignments,
            parallelizer_list=parallelizer_list,
            cmd_related_properties=cmd_related_properties,
        )

    ## Get the label of the node. By default, it is simply the name
    def get_dot_label(self) -> str:
        ## The name could be a full path
        name = self.cmd_invocation_with_io_vars.cmd_name
        basename = os.path.basename(str(name))

        wrapped_node_name = self.wrapped_node_name
        return f"{basename}({wrapped_node_name})"


def wrap_node(node: DFGNode, edges):
    r_wrap_bin = os.path.join(
        config.PASH_TOP, config.config["runtime"]["r_wrap_binary"]
    )

    ## At the moment we can only wrap a node that takes its input from stdin
    ## and outputs to stdout. Therefore the node needs to have only one input and one output.
    ## TO CHECK: with the remodelling also other cases should be handled
    inputs = node.get_input_list()
    assert len(inputs) == 1
    input_id = inputs[0]
    outputs = node.get_output_list()
    ## TODO: Would it make sense for outputs to be less than one?
    ## TODO: changed this from <= to == 1 to simplify reasoning later for now
    assert len(outputs) == 1
    output_id = outputs[0]
    access_map = {input_id: make_stream_input(), output_id: make_stream_output()}

    # create bash -c argument
    cmd_inv_with_io_vars: CommandInvocationWithIOVars = node.cmd_invocation_with_io_vars
    # do we need to copy here? currently, it seems fine
    cmd_inv_with_io_vars.remove_streaming_inputs()
    cmd_inv_with_io_vars.remove_streaming_outputs()
    # any non-streaming inputs or outputs are converted here already!
    cmd = to_arg_from_cmd_inv_with_io_vars_without_streaming_inputs_or_outputs_for_wrapping(
        cmd_inv_with_io_vars, edges
    )

    bash_command_arg = [Arg.string_to_arg("bash -c")]
    operand_list = bash_command_arg + [cmd]

    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=r_wrap_bin,
        flag_option_list=[],
        operand_list=operand_list,
        implicit_use_of_streaming_input=input_id,
        implicit_use_of_streaming_output=output_id,
        access_map=access_map,
    )

    ## TODO: It is not clear if it is safe to just pass redirections and assignments down the line as is
    redirs = node.com_redirs
    assignments = node.com_assignments

    return RWrap(
        cmd_inv_with_io_vars,
        com_redirs=redirs,
        com_assignments=assignments,
        wrapped_node_name=node.cmd_invocation_with_io_vars.cmd_name,
    )
