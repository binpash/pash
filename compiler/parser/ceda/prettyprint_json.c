#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <string.h>

#include "ast2b.h"


#include "json_tokener.h"


// 640KB ought to be enough for anybody
#define MAX_LINE_LENGTH (640 * 1024)


json_object* json_text_to_jobj (char* json_text) {
    // Nesting can be very deep e.g., scripts/intermediary/web-index_p2_1_funs.sh.json
    struct json_tokener* tok = json_tokener_new_ex (JSON_MAX_DEPTH);

    // Based on http://json-c.github.io/json-c/json-c-current-release/doc/html/json__tokener_8h.html
    json_object* jobj = json_tokener_parse_ex(tok, json_text, strlen (json_text));
    enum json_tokener_error jerr = json_tokener_get_error(tok);
    if (jerr != json_tokener_success) {
        fprintf (stderr, "Error: %s\n", json_tokener_error_desc (jerr));
        abort ();
    }

    assert (jobj != NULL);
    return (jobj);
}


int main (int argc, char* argv []) {
    FILE* fp = stdin;

    if (argc == 2) {
        if (strcmp (argv [1], "-") != 0) {
            fp = fopen (argv [1], "r");
            assert (fp != NULL);
        }
    }

    char* buf = malloc (MAX_LINE_LENGTH * sizeof (char));
    assert (buf != NULL);

    while (fgets (buf, 640 * 1024, fp) != NULL) {
        json_object* root = json_text_to_jobj (buf);
        const char* text = json_object_to_json_string_ext (root, JSON_C_TO_STRING_PRETTY);
        printf ("%s\n", text);
        putchar ('\n');
    }

    free (buf); // Pointless since we exit immediately

    return 0;
}
