#!/usr/bin/env python3

import os
import sys
import subprocess
from parse_to_ast2 import parse_to_ast

if 'PASH_TOP' in os.environ:
    PASH_TOP = os.environ['PASH_TOP']
else:
    GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
    PASH_TOP = subprocess.run(GIT_TOP_CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.rstrip()

sys.path.append(os.path.join(PASH_TOP, "compiler"))
from json_ast import serialize_asts_to_json


def main ():
    sys.setrecursionlimit (90001)

    if (len (sys.argv) == 1):
        inputPath = "-"
    else:
        inputPath = sys.argv [1]

    new_asts = []
    for output in (parse_to_ast (inputPath, True)):
        (new_ast, verbatim, linno_before, linno_after) = output;
        new_asts.append (new_ast);

        # Debugging
        if (False):
            print ("### Parsed lines [%d, %d)" % (linno_before, linno_after));
            print ("--------------------");
            print (verbatim, end='');
            print ("--------------------");

    json = serialize_asts_to_json (new_asts)
    print (json)

if __name__ == "__main__":
    main ()
