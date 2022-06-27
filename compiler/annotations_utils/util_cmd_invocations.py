import sys

from datatypes_new.BasicDatatypes import Flag
from datatypes_new.BasicDatatypesWithIO import OptionWithIO
from datatypes_new.CommandInvocationInitial import CommandInvocationInitial
from annotation_generation_new.datatypes.InputOutputInfo import InputOutputInfo
from annotation_generation_new.datatypes.ParallelizabilityInfo import ParallelizabilityInfo
from annotation_generation_new.datatypes.CommandProperties import CommandProperties
from annotation_generation_new.AnnotationGeneration import get_input_output_info_from_cmd_invocation, \
    get_parallelizability_info_from_cmd_invocation

from util import log

from config import get_path_annotation_repo
sys.path.insert(1, get_path_annotation_repo())

# for typing
from datatypes_new.CommandInvocationPrefix import CommandInvocationPrefix

from ir_utils import string_to_argument, redir_stdout_to_file, redir_file_to_stdin, make_command

def get_command_invocation_prefix_from_dfg_node(dfg_node):
    return CommandInvocationPrefix(cmd_name = dfg_node.com_name,
                                   flag_option_list = dfg_node.flag_option_list,
                                   positional_config_list = dfg_node.positional_config_list)

# TODO: ideally methods in the respective classes but requires refactoring of parsing infrastructure
def to_node_cmd_inv_with_io_vars(cmd_inv, edges, redirs, assignments):
    log("edges", edges)
    ast_cmd_name = string_to_argument(cmd_inv.cmd_name)
    log("ast_cmd_name", ast_cmd_name)
    ast_flagoptions = []
    for flagoption in cmd_inv.flag_option_list:
        ast_flagoptions += to_ast_flagoption(flagoption, edges)
    log("flagoptions", cmd_inv.flag_option_list)
    log("ast_flagoptions", ast_flagoptions)
    ast_operands = [to_ast_operand(operand, edges) for operand in cmd_inv.operand_list]
    log("operands", cmd_inv.operand_list)
    log("ast_operands", ast_operands)
    # log("type of ast_operands [0]", type(ast_operands[0])) # can only be used if there are operands
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
    log("node", node)
    return node

def to_ast_flagoption(flagoption, _edges):
    if isinstance(flagoption, Flag):
        return [string_to_argument(flagoption.get_name())]
    elif isinstance(flagoption, OptionWithIO): # retype to IOVar
        opt_name_ast = string_to_argument(flagoption.get_name())
        opt_arg_ast = translate_io_var_if_applicable(flagoption.get_arg())
        return [opt_name_ast, opt_arg_ast]

def to_ast_operand(operand, edges):
    return translate_io_var_if_applicable(operand, edges)

def translate_io_var_if_applicable(pot_io_var, edges):
    if isinstance(pot_io_var, int):
        return dereference_io_var(pot_io_var, edges)
    else:
        return to_ast_arg_string_type(pot_io_var)

def to_ast_arg_string_type(arg_string_type):
    return arg_string_type.get_name().arg_char_list # is of type Arg

# assumes io_var is an edge id
def dereference_io_var(io_var, edges):
    fid, _, _ = edges[io_var]
    log(fid)
    return fid.to_ast()

def get_input_output_info_from_cmd_invocation_util(cmd_invocationInitial : CommandInvocationInitial) -> InputOutputInfo:
    return get_input_output_info_from_cmd_invocation(cmd_invocationInitial)

def get_parallelizability_info_from_cmd_invocation_util(cmd_invocationInitial : CommandInvocationInitial) -> ParallelizabilityInfo:
    return get_parallelizability_info_from_cmd_invocation(cmd_invocationInitial)

def construct_property_container_from_list_of_properties(list_properties):
    return CommandProperties(dict(list_properties))

