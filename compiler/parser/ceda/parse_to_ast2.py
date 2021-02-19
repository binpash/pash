import sys;
from ctypes import *

from ast2a import of_node;
from dash2 import *;

# For testing
sys.path.append("/pash/compiler")

from ast2shell import to_string


LIBDASH_LIBRARY_PATH = "../libdash/src/.libs/libdash.so";
LIBDASH2_LIBRARY_PATH = "./libdash2.so"; # Uses shims in dash2.c


def parse_to_ast (inputPath, init):
    libdash2 = CDLL (LIBDASH2_LIBRARY_PATH);

    libdash = CDLL (LIBDASH_LIBRARY_PATH);

    if (init):
        initialize (libdash);

    if (inputPath == "-"):
        setinputtostdin (libdash);
    else:
        setinputfile (libdash, sys.argv [1]);

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
