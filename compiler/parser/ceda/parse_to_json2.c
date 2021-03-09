#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>

#include "nodes.h"
#include "parser.h"

#include "dash2.h"
#include "ast2a.h"
#include "ast2json.h"

#include "json.h"


#define TRUE  1
#define FALSE 0


void parse_args (void) {
}



void set_input_src (char* inputPath) {
    // TODO: allow file input
    if (inputPath == NULL) {
        Dash_setinputtostdin ();
    } else {
        Dash_setinputfile (inputPath, 0);
    }
}


void parse_all (void) {
    struct stackmark* smark = Dash_init_stack ();

    while (1) { // Dijkstra would not approve
        union node* n = Dash_parse_next (0); // not interactive

        if (n == NEOF) { // Dash.Done
            break;

            // []
        } else if (n == NERR) { // Dash.Error
            break;

            // raise Parse_error
        } else if (n == NULL) { // Dash.Null
//            printf ("null\n");
        } else { // Dash.Parsed
            struct t_TYPE* t = of_node (n);
            pour_the_t (t);

            printf ("\n");

            Dash_pop_stack (smark);
        }
    }
}


int main (int argc, char* argv []) {
    Dash_initialize ();

//    parse_args ();

    char* inputPath = NULL;

    if (argc == 2) {
        if (strcmp (argv [1], "-") != 0) {
            inputPath = argv [1];
        }
    }

    set_input_src (inputPath);

    parse_all ();

    // TODO: print_ast (JSON output)

    return 0;
}
