import sys;

from parse_to_ast2 import parse_to_ast;
from ast2shell import to_string;


sys.setrecursionlimit (30001)


if (len (sys.argv) == 1):
    new_asts = parse_to_ast ("-");
else:
    new_asts = parse_to_ast (sys.argv [1]);


#print ("NEW IMPLEMENTATION:");
#print (new_asts);
#print ("");

for ast in new_asts:
    print (to_string (ast));
