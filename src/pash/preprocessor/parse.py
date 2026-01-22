"""
Shell script parsing and unparsing utilities.
"""

import sys

from shell_ast.ast_util import UnparsedScript
from shasta.json_to_ast import to_ast_node
from shasta.bash_to_shasta_ast import to_ast_node as bash_to_shasta_ast

from util import log

import libdash.parser
import libbash


def parse_shell_to_asts(input_script_path, bash_mode=False):
    """Parse a shell script to a list of AST objects."""
    if bash_mode:
        return parse_shell_to_asts_bash(input_script_path)
    else:
        return parse_shell_to_asts_dash(input_script_path)


INITIALIZE_LIBDASH = True


def parse_shell_to_asts_dash(input_script_path):
    """Parse a POSIX shell script using libdash."""
    global INITIALIZE_LIBDASH
    try:
        new_ast_objects = libdash.parser.parse(input_script_path, INITIALIZE_LIBDASH)
        INITIALIZE_LIBDASH = False
        # Transform the untyped ast objects to typed ones
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


def parse_shell_to_asts_bash(input_script_path):
    """Parse a bash script using libbash."""
    try:
        new_ast_objects = libbash.bash_to_ast(input_script_path, with_linno_info=True)

        # Convert the libbash AST to a shasta AST
        typed_ast_objects = []
        for (
            untyped_ast,
            original_text,
            linno_before,
            linno_after,
        ) in new_ast_objects:
            typed_ast = bash_to_shasta_ast(untyped_ast)
            typed_ast_objects.append(
                (
                    typed_ast,
                    original_text.decode("utf-8", errors="replace"),
                    linno_before,
                    linno_after,
                )
            )
        return typed_ast_objects
    except RuntimeError as e:
        log("Parsing error!", e)
        sys.exit(1)


def from_ast_objects_to_shell(asts):
    """Convert AST objects back to shell script text."""
    shell_list = []
    for ast in asts:
        if isinstance(ast, UnparsedScript):
            shell_list.append(ast.text)
        else:
            shell_list.append(ast.pretty())

    shell_list = [
        x.decode("utf-8", errors="replace") if isinstance(x, bytes) else x
        for x in shell_list
    ]
    return "\n".join(shell_list) + "\n"
