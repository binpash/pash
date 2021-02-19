import sys;

from parse_to_ast2 import *;


def main ():
    if (len (sys.argv) == 1):
        new_asts = parse_to_ast ("-", True);
    else:
        new_asts = parse_to_ast (sys.argv [1], True);

    json = serialize_asts_to_json (new_asts);

    print (json);


if __name__ == "__main__":
    main ();
