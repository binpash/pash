#!/usr/bin/python3

import sys
import codecs

# export PYTHONIOENCODING=charmap


# sys.stdout = codecs.getwriter('charmap')(sys.stdout)


sys.path.append("/pash/compiler")

from parse import parse_shell, from_ir_to_shell, from_ir_to_shell_file
from json_ast import parse_json_ast_string, serialize_asts_to_json, json_to_shell
from ast2shell import to_string


if (len (sys.argv) != 2):
    print ("Usage: rt.py shell.sh");
    exit (1);

inputFile = sys.argv [1];

# json_ast_string = parse_shell (...) 
# json = parse_json_ast (inputFile);
# from_ir_to_shell

json = parse_shell (inputFile);
#print ("Shell script -> JSON:");
#print (json);

if (len (json) == 0):
    sys.exit (0);

asts = parse_json_ast_string (json);
#print ("JSON -> Pash AST:");
#print (asts);
#print ();

#print ("TODO: directly convert AST to shell script\n");

#json_rt = serialize_asts_to_json (asts)
#print ("JSON round-trip: %s" % json_rt);
#print ();


#shell_rt = json_to_shell (json_rt);
#print ("Shell round-trip: %s" % shell_rt);

#print ("to_string");
for ast in asts:
    str1 = to_string (ast);
    shell_direct = str1;
#    shell_direct = str.encode ("utf-8") + str1;

    # Some shell scripts have characters with ASCII value >= 128,
    # which disagrees with the default Python print.
    print (shell_direct);

#    sys.stdout.buffer.write (shell_direct.encode ('utf-8'));

#    newline = (10).to_bytes (1, byteorder='little');
#    sys.stdout.buffer.write (newline);

#    sys.stdout.buffer.write (

#    shell_direct_bytes = shell_direct.encode ('utf-8').strip ();
#
#    print (len (shell_direct_bytes));
#    for b in range (len (shell_direct_bytes)):
#        print ("%s" % chr (shell_direct_bytes [b]), end ="");

#    sys.stdout.write (shell_direct);
#    print (shell_direct);
#    print (shell_direct.encode ('utf8'));
#    sys.stdout.buffer.write (shell_direct);
