import sys;
from ctypes import *

from ast2a import of_node;
from dash2 import *;

# For testing
sys.path.append("/pash/compiler")

from ast2shell import to_string


#LIBDASH_LIBRARY_PATH = "../libdash/src/.libs/libdash.so";
LIBDASH2_LIBRARY_PATH = "./libdash2.so"; # Uses shims in dash2.c


def parse_to_ast (inputPath, init):
    libdash = CDLL (LIBDASH2_LIBRARY_PATH);

    if (init):
        libdash.Dash_initialize.argtypes = [];
        libdash.Dash_initialize.restypes = None;
        libdash.Dash_initialize ();

    if (inputPath == "-"):
        libdash.Dash_setinputtostdin ();
    else:
        libdash.Dash_setinputfile.argtypes = [c_char_p, c_int];
        libdash.Dash_setinputfile.restypes = c_int;
        libdash.Dash_setinputfile (sys.argv [1].encode ('utf-8'), 0);

    libdash.Dash_init_stack.argtypes = [];
    # libdash.Dash_init_stack.restype = POINTER (stackmark);
    libdash.Dash_init_stack.restype = c_void_p; # We want the opaque C pointer
    smark = libdash.Dash_init_stack ();

    libdash.Dash_parse_next.argtypes = [c_int];
    # libdash.Dash_parse_next.restype = POINTER (union_node);
    # We need the unadultered C pointer to test for NEOF and NERR
    # We cast to POINTER (union_node) later
    libdash.Dash_parse_next.restype = c_void_p;

    # libdash.Dash_pop_stack.argtypes = [POINTER (stackmark)];
    libdash.Dash_pop_stack.argtypes = [c_void_p];
    libdash.Dash_pop_stack.restype = None;

    NEOF = addressof (c_int.in_dll (libdash, "tokpushback"));
    NERR = addressof (c_int.in_dll (libdash, "lasttoken"));

    new_asts = [];

    while (True):
        n_ptr_C = libdash.Dash_parse_next (0);

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

            libdash.Dash_pop_stack (smark);

    return (new_asts);
