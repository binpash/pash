import os
import sys
from ctypes import *

from ast2a import of_node
from dash2 import *

if 'PASH_TOP' in os.environ:
    PASH_TOP = os.environ['PASH_TOP']
else:
    GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
    PASH_TOP = subprocess.run(GIT_TOP_CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.rstrip()


# TODO: use Pash root directory
LIBDASH_LIBRARY_PATH = os.path.join(PASH_TOP, "compiler/parser/libdash/src/.libs/libdash.so")


# This is a mix of dash.ml:parse_next and parse_to_json.ml.
def parse_to_ast (inputPath, init=True):
    libdash = CDLL (LIBDASH_LIBRARY_PATH);

    if (init):
        initialize (libdash);

    if (inputPath == "-"):
        setinputtostdin (libdash);
    else:
        setinputfile (libdash, inputPath);

    smark = init_stack (libdash);

    NEOF = addressof (c_int.in_dll (libdash, "tokpushback"));
    NERR = addressof (c_int.in_dll (libdash, "lasttoken"));

    new_asts = [];

    while (True):
        n_ptr_C = parsecmd_safe (libdash, False);

        if (n_ptr_C == None): # Dash.Null
            pass;
        elif (n_ptr_C == NEOF): # Dash.Done
            break;
        elif (n_ptr_C == NERR): # Dash.Error
            break;
        else:
            n_ptr = cast (n_ptr_C, POINTER (union_node));
            new_ast = of_node (n_ptr);

            new_asts.append (new_ast);

            pop_stack (libdash, smark);

    return (new_asts);
