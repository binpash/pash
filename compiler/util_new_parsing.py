# imports from annotation framework
import sys
sys.path.insert(1, "/home/felix/git-repos/MIT/annotations")
# for typing
from datatypes_new.CommandInvocation import CommandInvocation
# for use
from datatypes_new.BasicDatatypes import Option, ArgStringType
from parser_new.parser import parse

from ir_utils import format_arg_chars, string_to_argument, log


def merge_to_single_string_with_space(list_str):
    if len(list_str) == 1:
        return list_str[0]
    else:
        return " ".join(list_str)

def get_command_invocation(command, options) -> CommandInvocation:
    command_as_string: str = format_arg_chars(command)
    options_and_operands_as_string: str = merge_to_single_string_with_space([format_arg_chars(option) for option in options])
    command_invocation_as_string: str = f'{command_as_string} {options_and_operands_as_string}'
    command_invocation: CommandInvocation = parse(command_invocation_as_string)
    return command_invocation

def get_ast_for_flagoption(flagoption):
    result = string_to_argument(flagoption.get_name())
    if isinstance(flagoption, Option):
        # TODO: add argument here as well but eventually also fid
        assert False
    return result

def get_ast_for_argstringtype(arg):
    return string_to_argument(arg.get_name())

# TODO: this is a hack to fix the wrong parsing of "
def fix_parsing_newline(arg):
    if arg.get_name() == '\\n':
        return ArgStringType(r'"\n"')
    else:
        return arg


