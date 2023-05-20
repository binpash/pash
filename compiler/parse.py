import os
import config
import subprocess
import sys

from shell_ast.ast_util import UnparsedScript
from shell_ast.ast_node import AstNode, ast_node_to_untyped_deep
from shell_ast.untyped_to_ast import to_ast_node

from util import *

import libdash.parser
import libdash.printer

## Parses straight a shell script to an AST
## through python without calling it as an executable
def parse_shell_to_asts(input_script_path):
    try:
        new_ast_objects = libdash.parser.parse(input_script_path)

        ## Transform the untyped ast objects to typed ones
        typed_ast_objects = []
        for untyped_ast, original_text, linno_before, linno_after, in new_ast_objects:
             typed_ast = to_ast_node(untyped_ast)
             typed_ast_objects.append((typed_ast, original_text, linno_before, linno_after))

        return typed_ast_objects
    except libdash.parser.ParsingException as e:
        log("Parsing error!", e)
        sys.exit(1)

def parse_shell_to_asts_interactive(input_script_path: str):
    return libdash.parser.parse(input_script_path)

def from_ast_objects_to_shell(asts):
    shell_list = []
    for ast in asts:
        if(isinstance(ast, UnparsedScript)):
            shell_list.append(ast.text)
        else:
            ## We are working with two different abstractions for ASTs, one is the class and the other
            ## is its JSON object form. Due to Python's _disgusting_ lack of types (and our bad code)
            ## you can sometimes end up here with both.
            ##
            ## TODO: At some point this should be fixed and we should only work with the AstNode abstraction
            ##       and only serialize at the end. There is more info on that in ast_node.py
            if(isinstance(ast, AstNode)):
                serialized_ast = ast_node_to_untyped_deep(ast)
            else:
                serialized_ast = ast

            shell_list.append(libdash.printer.to_string(serialized_ast))
    return "\n".join(shell_list) + "\n"

def from_ast_objects_to_shell_file(asts, new_shell_filename):
    script = from_ast_objects_to_shell(asts)
    with open(new_shell_filename, 'w') as new_shell_file:
        new_shell_file.write(script)

def parse_shell(input_script_path):
    if(not os.path.isfile(input_script_path)):
        log("Error! File:", input_script_path, "does not exist.", level=0)
        sys.exit(1)
    parser_output = subprocess.run([config.PARSER_BINARY, input_script_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if (not parser_output.returncode == 0):
        log(parser_output.stderr)
        parser_output.check_returncode()
    return parser_output.stdout


## Simply wraps the ceda string_of_arg
def pash_string_of_arg(arg, quoted=False):
    return libdash.printer.string_of_arg(arg, quoted)

### Legacy

def from_ir_to_shell_legacy(ir_filename):
    printer_output = subprocess.run([config.PRINTER_BINARY, ir_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    printer_output.check_returncode()
    preprocessed_script = printer_output.stdout
    return preprocessed_script
