# imports from annotation framework
import sys
sys.path.insert(1, "/home/felix/git-repos/MIT/annotations")
# for typing
from datatypes_new.CommandInvocation import CommandInvocation
# for use
from parser_new.parser import parse

from ir_utils import format_arg_chars


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
