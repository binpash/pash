from typing import Set, List, Any

from definitions.ir.arg import Arg

from pash_annotations.datatypes.CommandInvocationInitial import CommandInvocationInitial
from pash_annotations.datatypes.BasicDatatypes import (
    Option,
    ArgStringType,
    Flag,
    Operand,
)
from pash_annotations.parser.parser import (
    parse,
    get_set_of_all_flags,
    get_dict_flag_to_primary_repr,
    get_set_of_all_options,
    get_dict_option_to_primary_repr,
    are_all_individually_flags,
)
from pash_annotations.parser.util_parser import get_json_data


from shell_ast.ast_util import format_arg_chars, string_to_argument


def merge_to_single_string_with_space(list_str):
    if len(list_str) == 1:
        return list_str[0]
    else:
        return " ".join(list_str)


def get_command_invocation(command, options) -> CommandInvocationInitial:
    command_as_string: str = format_arg_chars(command)
    options_and_operands_as_string: str = merge_to_single_string_with_space(
        [format_arg_chars(option) for option in options]
    )
    command_invocation_as_string: str = (
        f"{command_as_string} {options_and_operands_as_string}"
    )
    command_invocation: CommandInvocationInitial = parse(command_invocation_as_string)
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
    if arg.get_name() == "\\n":
        return ArgStringType(r'"\n"')
    else:
        return arg


def parse_arg_list_to_command_invocation(
    command, flags_options_operands
) -> CommandInvocationInitial:
    cmd_name = format_arg_chars(command)
    json_data = get_json_data(cmd_name)

    set_of_all_flags: Set[str] = get_set_of_all_flags(json_data)
    dict_flag_to_primary_repr: dict[str, str] = get_dict_flag_to_primary_repr(json_data)
    set_of_all_options: Set[str] = get_set_of_all_options(json_data)
    dict_option_to_primary_repr: dict[str, str] = get_dict_option_to_primary_repr(
        json_data
    )
    # we keep the Arg for everything but flag and option names

    # parse list of command invocation terms
    flag_option_list: List[Any] = []
    i = 0
    while i < len(flags_options_operands):
        potential_flag_or_option_arg = flags_options_operands[i]
        potential_flag_or_option_name = format_arg_chars(potential_flag_or_option_arg)
        if potential_flag_or_option_name in set_of_all_flags:
            flag_name_as_string: str = dict_flag_to_primary_repr.get(
                potential_flag_or_option_name, potential_flag_or_option_name
            )
            flag: Flag = Flag(flag_name_as_string)
            flag_option_list.append(flag)
        elif (potential_flag_or_option_name in set_of_all_options) and (
            (i + 1) < len(flags_options_operands)
        ):
            option_name_as_string: str = dict_option_to_primary_repr.get(
                potential_flag_or_option_name, potential_flag_or_option_name
            )
            option_arg_as_arg: Arg = Arg(flags_options_operands[i + 1])
            option = Option(option_name_as_string, option_arg_as_arg)
            flag_option_list.append(option)
            i += 1  # since we consumed another term for the argument
        elif (
            potential_flag_or_option_name == "-"
        ):  # switch to operand mode (interpreted as hyphen-stdin)
            break
        elif are_all_individually_flags(
            potential_flag_or_option_name, set_of_all_flags
        ):
            for split_el in list(potential_flag_or_option_name[1:]):
                flag: Flag = Flag(f"-{split_el}")
                flag_option_list.append(flag)
        else:
            break  # next one is Operand, and we keep these in separate list
        i += 1

    # we would probably want to skip '--' but then the unparsed command could have a different meaning so we'd need to keep it
    # for now, omitted
    # if parsed_elements_list[i] == '--':
    #     i += 1

    operand_list = [
        Operand(Arg(operand_arg)) for operand_arg in flags_options_operands[i:]
    ]
    # log("type of operand_list[0].get_name()", type(operand_list[0].get_name()))   can only be used if there are operands

    return CommandInvocationInitial(cmd_name, flag_option_list, operand_list)
