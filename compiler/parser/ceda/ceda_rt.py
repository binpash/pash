import sys;

from parse_to_ast2 import parse_to_ast;
from ast2shell import to_string;

import cProfile;


sys.setrecursionlimit (9001);


def print_asts (new_asts):
    for ast in new_asts:
        print (to_string (ast));
#        to_string (ast);


init = True

#cProfile.runctx ("parse_to_ast (sys.argv [1], init)", globals (), locals ());
#sys.exit (0);

for f in range (10):
    if (len (sys.argv) == 1):
        new_asts = parse_to_ast ("-", init);
    else:
        new_asts = parse_to_ast (sys.argv [1], init);

    init = False;

#print ("NEW IMPLEMENTATION:");
#print (new_asts);
#print ("");

    print_asts (new_asts);

#    cProfile.runctx ("print_asts (new_asts)", globals (), locals ());
