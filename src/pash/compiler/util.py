"""
Utility functions for the PaSh compiler.
"""

from datetime import timedelta
import logging
from typing import Optional, TypeVar, Union, List, Any

TType = TypeVar("TType")
import os
import tempfile

import config
from shasta.ast_node import (
        CArgChar,
        VArgChar,
        QArgChar,
        CommandNode,
        BackgroundNode,
        SubshellNode,
        SemiNode,
        DefunNode,
        AssignNode,
        FileRedirNode,
        )


# === List utilities ===


def flatten_list(lst):
    return [item for sublist in lst for item in sublist]


# === Logging utilities ===


def print_time_delta(prefix, start_time, end_time):
    """Log time delta between start and end time."""
    time_difference = (end_time - start_time) / timedelta(milliseconds=1)
    log("{} time:".format(prefix), time_difference, " ms")

def log(*args, end="\n", level=2):
    """Wrapper for logging."""
    if level == 1:
        concatted_args = " ".join([str(a) for a in list(args)])
        logging.warning(f"{config.LOGGING_PREFIX} {concatted_args}")
    elif level >= 2:
        concatted_args = " ".join([str(a) for a in list(args)])
        logging.info(f"{config.LOGGING_PREFIX} {concatted_args}")


# === File utilities ===


def ptempfile():
    fd, name = tempfile.mkstemp(dir=config.PASH_TMP_PREFIX)
    ## TODO: Get a name without opening the fd too if possible
    os.close(fd)
    return name


# === Optional/default utilities ===


def return_empty_list_if_none_else_itself(
    arg: Optional[TType],
) -> Union[TType, List[Any]]:  # list always empty
    if arg is None:
        return []
    else:
        return arg


def return_default_if_none_else_itself(arg: Optional[TType], default: TType) -> TType:
    if arg is None:
        return default
    else:
        return arg


# === AST key-value utilities ===


## This function gets a key and a value from the ast json format
def get_kv(dic):
    return (dic[0], dic[1])


def make_kv(key, val):
    return [key, val]


# === Custom compiler error codes ===


class UnparallelizableError(Exception):
    pass


class AdjLineNotImplementedError(Exception):
    pass


# to be raised in pash_compiler if a UnparallelizableError is caught at any point running the compiler
#   primarily to differentiate
#       --assert_compiler_success (exit with error only under general exceptions caught)
#       --assert_all_regions_parallelizable (exit with error when regions are found not parallelizable + general exceptions)
class NotAllRegionParallelizableError(Exception):
    pass


# === AST classes ===


class UnparsedScript:
    """
    Represents text that was not modified at all by preprocessing,
    and therefore does not need to be unparsed.
    """

    def __init__(self, text):
        self.text = text


# === AST pattern matching and construction ===


def format_arg_chars(arg_chars):
    chars = [arg_char.format() for arg_char in arg_chars]
    return "".join(chars)


def string_to_carg_char_list(string: str) -> "list[CArgChar]":
    ret = [CArgChar(ord(char)) for char in string]
    return ret


def string_to_arguments(string):
    return [string_to_argument(word) for word in string.split(" ")]


def string_to_argument(string):
    ret = [char_to_arg_char(char) for char in string]
    return ret


def char_to_arg_char(char):
    return CArgChar(ord(char))


def standard_var_ast(string):
    return VArgChar("Normal", False, string, [])


def make_quoted_variable(string):
    return QArgChar(arg=[standard_var_ast(string)])


def quote_arg(arg):
    return QArgChar(arg=arg)


def redir_append_stderr_to_string_file(string):
    return FileRedirNode("Append", ("fixed", 2), string_to_argument(string))

def redir_stdout_to_file(arg):
    return FileRedirNode("To", ("fixed", 1), arg)


def redir_file_to_stdin(arg):
    return FileRedirNode("From", ("fixed", 0), arg)

def make_background(body, redirections=None):
    redirections = [] if redirections is None else redirections
    lineno = 0
    node = BackgroundNode(lineno, body, redirections)
    return node


def make_subshell(body, redirections=None):
    redirections = [] if redirections is None else redirections
    lineno = 0
    node = SubshellNode(lineno, body, redirections)
    return node


def make_command(arguments, redirections=None, assignments=None):
    redirections = [] if redirections is None else redirections
    assignments = [] if assignments is None else assignments
    lineno = 0
    node = CommandNode(lineno, assignments, arguments, redirections)
    return node


def make_assignment(var, value):
    lineno = 0
    assignment = AssignNode(var, value)
    return CommandNode(lineno, [assignment],[],[])



def make_semi_sequence(asts):
    if len(asts) == 0:
        return make_command([string_to_argument(":")])

    if len(asts) == 1:
        return asts[0]
    else:
        acc = asts[-1]
        # Remove the last ast
        iter_asts = asts[:-1]
        for ast in iter_asts[::-1]:
            acc = SemiNode(ast, acc)
        return acc


def make_defun(name, body):
    lineno = 0
    node = DefunNode(lineno, string_to_argument(name), body )
    return node
