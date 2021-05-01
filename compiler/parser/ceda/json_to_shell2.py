#!/usr/bin/env python3

import os
import sys

if 'PASH_TOP' in os.environ:
    PASH_TOP = os.environ['PASH_TOP']
else:
    GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
    PASH_TOP = subprocess.run(GIT_TOP_CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.rstrip()

sys.path.append(os.path.join(PASH_TOP, "compiler"))

from json_ast import parse_json_ast, parse_json_ast_string
from ast2shell import *


def main ():
    if (len (sys.argv) == 1):
        asts = parse_json_ast ("/dev/stdin")
    else:
        asts = parse_json_ast (sys.argv [1])

    for ast in asts:
        print (to_string (ast))

def json_string_to_shell_string(json_string):
    asts = parse_json_ast_string(json_string)
    shell_list = []
    for ast in asts:
        shell_list.append(to_string(ast))
    return "\n".join(shell_list) + "\n"

def json_to_shell_string(input_filename):
    with open(input_filename) as json_file:
        json_string = json_file.read()
        return json_string_to_shell_string(json_string)

if __name__ == "__main__":
    main ()
