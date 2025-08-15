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
        # log("Ast:", ast)
        if(isinstance(ast, UnparsedScript)):
            shell_list.append(ast.text)
        else:
            if os.environ.get("SERVERLESS_PASH", "0") == "0":
                # Original one
                shell_list.append(ast.pretty())
            else:
                # Shortcut for cat and output redirection
                line = ast.pretty()
                if "cat" in line and "inputs" in line:
                    cat_list = line.split()
                    obj_key = cat_list[2]
                    s3_cmd = f"python3 {os.path.join(config.PASH_TOP, 'runtime/serverless/aws/s3-get-object.py')} {cat_list[2]} /dev/stdout"
                    new_cmd = line.replace(f"cat {obj_key}", s3_cmd)
                    shell_list.append(new_cmd)
                elif "stdout" in line:
                    stdout_str = line.split()
                    obj_key = stdout_str[-3].replace(">", "")
                    s3_cmd = f"| python3 {os.path.join(config.PASH_TOP, 'runtime/serverless/aws/s3-put-object.py')} {obj_key} /dev/stdin dummy"
                    new_cmd = line.replace(f">{obj_key}", s3_cmd)
                    shell_list.append(new_cmd)
                elif "auto-split" in line:
                    shell_list.append(line.replace("runtime", config.PASH_TOP + "/runtime"))
                else:
                    shell_list.append(line)
    return "\n".join(shell_list) + "\n"

def from_ast_objects_to_shell_file(asts, new_shell_filename):
    script = from_ast_objects_to_shell(asts)
    with open(new_shell_filename, 'w') as new_shell_file:
        new_shell_file.write(script)

## Simply wraps the string_of_arg
def pash_string_of_arg(arg, quoted=False):
    return string_of_arg(arg, quoted)
