#!/usr/bin/env python3

import os
import sys
from parse_to_ast2 import parse_to_ast

if 'PASH_TOP' in os.environ:
    PASH_TOP = os.environ['PASH_TOP']
else:
    GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
    PASH_TOP = subprocess.run(GIT_TOP_CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.rstrip()

sys.path.append(os.path.join(PASH_TOP, "compiler"))
from json_ast import serialize_asts_to_json


def main ():
    if (len (sys.argv) == 1):
        new_asts = parse_to_ast ("-", True)
    else:
        new_asts = parse_to_ast (sys.argv [1], True)

    json = serialize_asts_to_json (new_asts)

    print (json)


if __name__ == "__main__":
    main ()
