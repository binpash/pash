from pash_annotations.datatypes.BasicDatatypes import Flag, ArgStringType, Operand
from pash_annotations.datatypes.BasicDatatypesWithIO import OptionWithIO
from pash_annotations.datatypes.CommandInvocationInitial import CommandInvocationInitial
from pash_annotations.annotation_generation.datatypes.InputOutputInfo import (
    InputOutputInfo,
)
from pash_annotations.annotation_generation.datatypes.ParallelizabilityInfo import (
    ParallelizabilityInfo,
)
from pash_annotations.annotation_generation.datatypes.CommandProperties import (
    CommandProperties,
)
from pash_annotations.annotation_generation.AnnotationGeneration import (
    get_input_output_info_from_cmd_invocation,
    get_parallelizability_info_from_cmd_invocation,
)
from pash_annotations.datatypes.CommandInvocationWithIOVars import (
    CommandInvocationWithIOVars,
)

from definitions.ir.arg import Arg

# for typing
from pash_annotations.datatypes.CommandInvocationPrefix import CommandInvocationPrefix

from shell_ast.ast_util import (
    string_to_argument,
    redir_stdout_to_file,
    redir_file_to_stdin,
    make_command,
)


def get_command_invocation_prefix_from_dfg_node(dfg_node):
    return CommandInvocationPrefix(
        cmd_name=dfg_node.com_name,
        flag_option_list=dfg_node.flag_option_list,
        positional_config_list=dfg_node.positional_config_list,
    )


# TODO: ideally methods in the respective classes but requires refactoring of parsing infrastructure
# TODO: isn't this `to_ast`?
def to_node_cmd_inv_with_io_vars(cmd_inv, edges, redirs, assignments):
    ast_cmd_name = string_to_argument(cmd_inv.cmd_name)
    ast_flagoptions = []
    for flagoption in cmd_inv.flag_option_list:
        ast_flagoptions += to_ast_flagoption(flagoption, edges)
    ast_operands = [to_ast_operand(operand, edges) for operand in cmd_inv.operand_list]
    cmd_asts = [ast_cmd_name] + ast_flagoptions + ast_operands

    # TODO: check for actual stdin
    stdin_redir = []
    if cmd_inv.implicit_use_of_streaming_input is not None:
        fid, _, _ = edges[cmd_inv.implicit_use_of_streaming_input]
        if not (fid.has_file_descriptor_resource() and fid.resource.is_stdin()):
            stdin_redir = [redir_file_to_stdin(fid.to_ast())]

    # TODO: check for actual stdout
    stdout_redir = []
    if cmd_inv.implicit_use_of_streaming_output is not None:
        fid, _, _ = edges[cmd_inv.implicit_use_of_streaming_output]
        if not (fid.has_file_descriptor_resource() and fid.resource.is_stdout()):
            stdout_redir = [redir_stdout_to_file(fid.to_ast())]

    new_redirs = redirs + stdin_redir + stdout_redir
    node = make_command(cmd_asts, redirections=new_redirs, assignments=assignments)
    return node


def to_ast_flagoption(flagoption, edges):
    if isinstance(flagoption, Flag):
        return [string_to_argument(flagoption.get_name())]
    elif isinstance(flagoption, OptionWithIO):  # retype to IOVar
        opt_name_ast = string_to_argument(flagoption.get_name())
        opt_arg_ast = translate_io_var_if_applicable(flagoption.get_arg(), edges)
        return [opt_name_ast, opt_arg_ast]


def to_ast_operand(operand, edges):
    if isinstance(operand, Operand):
        return translate_io_var_if_applicable(operand.get_name(), edges)
    return translate_io_var_if_applicable(operand, edges)


def translate_io_var_if_applicable(pot_io_var, edges):
    # TODO: this is currently a hack but eventually every possible type gets their own to_ast-function
    if isinstance(pot_io_var, int):
        return dereference_io_var(pot_io_var, edges)
    elif isinstance(pot_io_var, ArgStringType):
        return to_ast_arg_string_type(pot_io_var)
    elif isinstance(pot_io_var, CommandInvocationWithIOVars):
        assert False
        # only happens as r-wrapped node
        return to_node_cmd_inv_with_io_vars(pot_io_var, edges, [], [])
    elif isinstance(pot_io_var, Arg):
        return pot_io_var.to_ast()
    else:
        raise Exception("Unhandled type for operand in to_ast!")


def to_ast_arg_string_type(arg_string_type):
    return arg_string_type.get_name().arg_char_list  # is of type Arg


# assumes io_var is an edge id
def dereference_io_var(io_var, edges):
    fid, _, _ = edges[io_var]
    return fid.to_ast()


def get_input_output_info_from_cmd_invocation_util(
    cmd_invocationInitial: CommandInvocationInitial,
) -> InputOutputInfo:
    return get_input_output_info_from_cmd_invocation(cmd_invocationInitial)


def get_parallelizability_info_from_cmd_invocation_util(
    cmd_invocationInitial: CommandInvocationInitial,
) -> ParallelizabilityInfo:
    return get_parallelizability_info_from_cmd_invocation(cmd_invocationInitial)


def construct_property_container_from_list_of_properties(list_properties):
    return CommandProperties(dict(list_properties))


# this function is needed to wrap a node in `r_wrap`
def to_arg_from_cmd_inv_with_io_vars_without_streaming_inputs_or_outputs_for_wrapping(
    cmd_inv, edges
):
    # we already expand here
    whole_cmd = Arg.string_to_arg("'")
    arg_cmd_name = Arg.string_to_arg(cmd_inv.cmd_name)
    arg_flagoptions = []
    for flagoption in cmd_inv.flag_option_list:
        arg_flagoptions += to_arg_flagoption(flagoption, edges)
    arg_operands = [to_arg_operand(operand, edges) for operand in cmd_inv.operand_list]
    all_cmd_parts_arg = [arg_cmd_name]
    all_cmd_parts_arg.extend(arg_flagoptions)
    all_cmd_parts_arg.extend(arg_operands)
    for part in all_cmd_parts_arg:
        whole_cmd.concatenate(part)
    whole_cmd.concatenate(Arg.string_to_arg("'"))
    return whole_cmd


def to_arg_flagoption(flagoption, edges):
    if isinstance(flagoption, Flag):
        return [Arg.string_to_arg(flagoption.get_name())]
    elif isinstance(flagoption, OptionWithIO):
        opt_name_arg = Arg.string_to_arg(flagoption.get_name())
        opt_arg_arg = translate_io_var_to_arg_if_applicable(flagoption.get_arg(), edges)
        return [opt_name_arg, opt_arg_arg]


def to_arg_operand(operand, edges):
    if isinstance(operand, Operand):
        return translate_io_var_to_arg_if_applicable(operand.get_name(), edges)
    return translate_io_var_to_arg_if_applicable(operand, edges)


def translate_io_var_to_arg_if_applicable(pot_io_var, edges):
    if isinstance(pot_io_var, int):
        return Arg(dereference_io_var(pot_io_var, edges))
    elif isinstance(pot_io_var, ArgStringType):
        result = pot_io_var.get_name()  # is of type Arg
        return result
    else:
        raise Exception("Unhandled type for operand in to_arg!")
