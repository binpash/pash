#!/usr/bin/python3


import sys;

sys.path.append("/pash/compiler");
from json_ast import parse_json_ast;


from ast2shell import *;


def main ():
    if (len (sys.argv) == 1):
        asts = parse_json_ast ("/dev/stdin");
    else:
        asts = parse_json_ast (sys.argv [1]);

    for ast in asts:
        print (to_string (ast));


if __name__ == "__main__":
    main ();
