import os
import config
import subprocess
import sys

from shell_ast.ast_util import UnparsedScript
from shasta.ast_node import ast_node_to_untyped_deep
from shasta.json_to_ast import to_ast_node
from shasta.ast_node import string_of_arg

from util import *

import libdash.parser


## Parses straight a shell script to an AST
## through python without calling it as an executable
INITIALIZE_LIBDASH = True
def parse_shell_to_asts(input_script_path):
    global INITIALIZE_LIBDASH
    try:
        new_ast_objects = libdash.parser.parse(input_script_path, INITIALIZE_LIBDASH)
        INITIALIZE_LIBDASH = False
        ## Transform the untyped ast objects to typed ones
        typed_ast_objects = []
        for (
            untyped_ast,
            original_text,
            linno_before,
            linno_after,
        ) in new_ast_objects:
            typed_ast = to_ast_node(untyped_ast)
            typed_ast_objects.append(
                (typed_ast, original_text, linno_before, linno_after)
            )

        return typed_ast_objects
    except libdash.parser.ParsingException as e:
        log("Parsing error!", e)
        sys.exit(1)


def parse_shell_to_asts_interactive(input_script_path: str):
    return libdash.parser.parse(input_script_path)


def from_ast_objects_to_shell(asts):
    shell_list = []
    for ast in asts:
        # log("Ast:", ast)
        if isinstance(ast, UnparsedScript):
            shell_list.append(ast.text)
        else:
            shell_list.append(ast.pretty())
    return "\n".join(shell_list) + "\n"


def from_ast_objects_to_shell_file(asts, new_shell_filename):
    script = from_ast_objects_to_shell(asts)
    with open(new_shell_filename, "w") as new_shell_file:
        new_shell_file.write(script)


## Simply wraps the string_of_arg
def pash_string_of_arg(arg, quoted=False):
    return string_of_arg(arg, quoted)
