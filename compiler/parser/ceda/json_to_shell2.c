#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <string.h>

#include "ast2b.h"


// 640KB ought to be enough for anybody
#define MAX_LINE_LENGTH (640 * 1024)


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
        json_text_to_string (buf);
        putchar ('\n');
    }

    free (buf); // Pointless since we exit immediately

    return 0;
}
