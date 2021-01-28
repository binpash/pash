/*


        ast2b.c : Exports 'to_string' which, analogously to the function in
                  libdash's ast.ml, converts a JSON representation of the
                  libdash AST into an ordinary shell script.

                  The functions here intentionally closely match ast.ml;
                  relevant OCaml snippets are included as comments.


*/


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "ast2b.h"

#include "json_tokener.h"


#define TRUE  1
#define FALSE 0


// Helper functions for debugging JSON
static void debug_jsonI (json_object* jobj, int depth);

// Helper functions for pattern-matching JSON
static int JSONMatchesStr (json_object* obj, char* magicWord);
static int JSONArrayStartsWithStr (json_object* obj, char* magicWord);
static int isJSONArrayOfLength (json_object* obj, int len);
static int isEmptyJSONArray (json_object* obj);

// Helper functions that move out complicated cases from 'to_string'
static void string_of_commandI (json_object* assigns, json_object* cmds, json_object* redirs);
static void string_of_pipeI (json_object* bg, json_object* ps);

// Helper functions that correspond fairly closing to ast.ml.
// Note that some simple functions (e.g., separated/background)
// have been manually inlined.
static void string_of_var_type (json_object* vt);
static void show_unless (int expected, int actual);
static void string_of_if (json_object* c, json_object* t, json_object* e);
static void string_of_arg_char (json_object* first, json_object* second);
static void string_of_arg (json_object* arg1);
static void string_of_assign (json_object* va);
static void string_of_case (json_object* c);
static void string_of_redir (json_object* r);
static void string_of_redirs (json_object* rs);


// ----------------------------------------------------------------------------


void debug_json (json_object* jobj) {
    printf ("\n");
    printf ("-----\n");
    debug_jsonI (jobj, 0);
    printf ("-----\n");
    printf ("\n");
}


// ----------------------------------------------------------------------------


static void debug_jsonI (json_object* jobj, int depth) {
    for (int i = 0; i < depth; i++) {
        printf ("    ");
    }

    switch (json_object_get_type (jobj)) {
        case json_type_null:
            printf ("null\n");
            break;

        case json_type_boolean:
            printf ("boolean\n");
            break;

        case json_type_double:
            printf ("double\n");
            break;

        case json_type_int:
            printf ("int: %d\n", json_object_get_int (jobj));
            break;

        case json_type_object:
            printf ("object\n");
            break;

        case json_type_array:
            printf ("array\n");

            for (int i = 0; i < json_object_array_length (jobj); i++) {
                debug_jsonI (json_object_array_get_idx (jobj, i), depth + 1);
            }
            break;

        case json_type_string:
            (void) 0; // Avoid compiler 'error: a label can only be part of a statement'
            const char* str = json_object_get_string (jobj);

            printf ("string: %s\n", str);
            break;

        default:
            assert ("! Unexpected case");
            break;
    }
}


// string: <magicWord>
static int JSONMatchesStr (json_object* obj, char* magicWord) {
    return    (json_object_get_type (obj) == json_type_string)
           && (strcmp (json_object_get_string (obj), magicWord) == 0);
}


/*
   array
       string: <magicWord>
       ...
*/
static int JSONArrayStartsWithStr (json_object* obj, char* magicWord) {
    return    (json_object_get_type (obj) == json_type_array)
           && (json_object_array_length (obj) >= 1)
           && (JSONMatchesStr (json_object_array_get_idx (obj, 0), magicWord));
}


static int isJSONArrayOfLength (json_object* obj, int len) {
    return    (json_object_get_type (obj) == json_type_array)
           && (json_object_array_length (obj) == len);
}


// array
static int isEmptyJSONArray (json_object* obj) {
    return isJSONArrayOfLength (obj, 0);
}


// ----------------------------------------------------------------------------


int json_object_get_flex_boolean (json_object* bg) {
    int truth = FALSE;

    if (json_object_get_type (bg) == json_type_boolean) {
        truth = json_object_get_boolean (bg);
    } else if (json_object_get_type (bg) == json_type_string) {
        if (strcmp (json_object_get_string (bg), "false") == 0) {
            truth = FALSE;
        } else if (strcmp (json_object_get_string (bg), "true") == 0) {
            truth = TRUE;
        } else {
            abort ();
        }
    }

    return (truth);
}


// string_of_commandI contains code moved out of the OCaml 'to_string' for clarity
/*
    | Command (_,assigns,cmds,redirs) ->
       separated string_of_assign assigns ^
       (if List.length assigns = 0 || List.length cmds = 0 then "" else " ") ^
       separated string_of_arg cmds ^ string_of_redirs redirs
*/
void string_of_commandI (json_object* assigns, json_object* cmds, json_object* redirs) {
    assert (assigns != NULL);
    assert (cmds != NULL);
    assert (redirs != NULL);

    assert (json_object_get_type (assigns) == json_type_array);
    assert (json_object_get_type (cmds) == json_type_array);
    assert (json_object_get_type (redirs) == json_type_array);

    for (int i = 0; i < json_object_array_length (assigns); i++) {
        if (i > 0) {
            putchar (' '); // Inline 'separated'
        }
        string_of_assign (json_object_array_get_idx (assigns, i));
    }

    if (   (json_object_array_length (assigns) != 0)
        && (json_object_array_length (cmds) != 0)) { // De Morgan's law
        putchar (' ');
    }

    for (int i = 0; i < json_object_array_length (cmds); i++) {
        if (i > 0) {
            putchar (' '); // Inline 'separated'
        }
        string_of_arg (json_object_array_get_idx (cmds, i));
    }

    string_of_redirs (redirs);
}


// string_of_pipeI contains code moved out of the OCaml 'to_string' for clarity
/*
    background s = "{ " ^ s ^ " & }"

    | Pipe (bg,ps) ->
       let p = intercalate " | " (List.map to_string ps) in
       if bg then background p else p
*/
void string_of_pipeI (json_object* bg, json_object* ps) {
    assert (bg != NULL);
    assert (ps != NULL);

    assert (json_object_get_type (ps) == json_type_array);

    if (json_object_get_flex_boolean (bg)) { // Inline 'background'
        printf ("{ ");
    }

    for (int i = 0; i < json_object_array_length (ps); i++) {
        if (i > 0) {
            printf (" | "); // Inline 'intercalate " | "'
        }
        to_string (json_object_array_get_idx (ps, i));
    }

    if (json_object_get_flex_boolean (bg)) { // Inline 'background'
        printf (" & }");
    }
}


// ----------------------------------------------------------------------------


/*
    let string_of_var_type = function
     | Normal -> ""
     | Minus -> "-"
     | Plus -> "+"
     | Question -> "?"
     | Assign -> "="
     | TrimR -> "%"
     | TrimRMax -> "%%"
     | TrimL -> "#"
     | TrimLMax -> "##"
     | Length -> "#"
*/
void string_of_var_type (json_object* vt) {
    assert (vt != NULL);

    assert (json_object_get_type (vt) == json_type_string);

    const char* vt_str = json_object_get_string (vt);
    if (strcmp (vt_str, "Normal") == 0) {
        // No-op
    } else if (strcmp (vt_str, "Minus") == 0) {
        putchar ('-');
    } else if (strcmp (vt_str, "Plus") == 0) {
        putchar ('+');
    } else if (strcmp (vt_str, "Question") == 0) {
        putchar ('?');
    } else if (strcmp (vt_str, "Assign") == 0) {
        putchar ('=');
    } else if (strcmp (vt_str, "TrimR") == 0) {
        putchar ('%');
    } else if (strcmp (vt_str, "TrimRMax") == 0) {
        putchar ('%');
        putchar ('%');
    } else if (strcmp (vt_str, "TrimL") == 0) {
        putchar ('#');
    } else if (strcmp (vt_str, "TrimLMax") == 0) {
        putchar ('#');
        putchar ('#');
    } else if (strcmp (vt_str, "Length") == 0) {
        putchar ('#');
    } else {
        debug_json (vt);

        assert (! "Unexpected type for string_of_var_type\n");
    }
}


/*
void separated (void) {
    assert (! "This function has been manually inlined\n");
}
*/


/*
    let show_unless expected actual =
      if expected = actual
      then ""
      else string_of_int actual
*/
void show_unless (int expected, int actual) {
    if (expected != actual) {
        printf ("%d", actual);
    }
}


/*
void background (void) {
    assert (! "This function has been manually inlined\n");
}
*/


void json_text_to_string (char* json_text) {
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
    to_string (jobj);
}


// OCaml function body is marvelous but this comment is too narrow to contain it
void to_string (json_object* jobj) {
    switch (json_object_get_type (jobj)) {
        case json_type_array:
            assert (json_object_array_length (jobj) == 2);

            json_object* obj1 = json_object_array_get_idx (jobj, 0);
            assert (obj1 != NULL);
            assert (json_object_get_type (obj1) == json_type_string);

            const char* str1 = json_object_get_string (obj1);
            assert (str1 != NULL);

            json_object* obj2 = json_object_array_get_idx (jobj, 1);
            assert (obj2 != NULL);
            assert (json_object_get_type (obj2) == json_type_array);

            if (strcmp (str1, "Command") == 0) {
                assert (json_object_array_length (obj2) == 4);

                string_of_commandI (// First param (lineno) is ignored!
                                    json_object_array_get_idx (obj2, 1),   // assigns
                                    json_object_array_get_idx (obj2, 2),   // cmds
                                    json_object_array_get_idx (obj2, 3));  // redirs
            } else if (strcmp (str1, "Pipe") == 0) {
                assert (json_object_array_length (obj2) == 2);

                string_of_pipeI (json_object_array_get_idx (obj2, 0),
                                 json_object_array_get_idx (obj2, 1));
            } else if (strcmp (str1, "Redir") == 0) {
                /*
                    | Redir (_,a,redirs) ->
                       to_string a ^ string_of_redirs redirs
                */

                assert (json_object_array_length (obj2) == 3);
                // First param is ignored
                json_object* a      = json_object_array_get_idx (obj2, 1);
                json_object* redirs = json_object_array_get_idx (obj2, 2);

                to_string (a);
                string_of_redirs (redirs);
            } else if (strcmp (str1, "Background") == 0) {
                // | Background (_,a,redirs) ->
                //    background (to_string a ^ string_of_redirs redirs)

                assert (json_object_array_length (obj2) == 3);
                // First param is ignored
                json_object* a      = json_object_array_get_idx (obj2, 1);
                json_object* redirs = json_object_array_get_idx (obj2, 2);

                printf ("{ "); // Inline 'background'
                to_string (a);
                string_of_redirs (redirs);
                printf (" & }"); // Inline 'background'
            } else if (strcmp (str1, "Subshell") == 0) {
                //   | Subshell (_,a,redirs) ->
                //      parens (to_string a ^ string_of_redirs redirs)
                //
                // let parens s = "( " ^ s ^ " )"

                assert (json_object_array_length (obj2) == 3);
                // First param is ignored
                json_object* a      = json_object_array_get_idx (obj2, 1);
                json_object* redirs = json_object_array_get_idx (obj2, 2);

                printf ("( "); // Inline 'parens'
                to_string (a);
                string_of_redirs (redirs);
                printf (" )"); // Inline 'parens'
            } else if (strcmp (str1, "And") == 0) {
                // | And (a1,a2) -> to_string  a1 ^ " && " ^ to_string a2

                assert (json_object_array_length (obj2) == 2);
                json_object* a1 = json_object_array_get_idx (obj2, 0);
                json_object* a2 = json_object_array_get_idx (obj2, 1);

                to_string (a1);
                printf (" && ");
                to_string (a2);
            } else if (strcmp (str1, "Or") == 0) {
                // | Or (a1,a2) -> to_string a1 ^ " || " ^ to_string a2

                assert (json_object_array_length (obj2) == 2);
                json_object* a1 = json_object_array_get_idx (obj2, 0);
                json_object* a2 = json_object_array_get_idx (obj2, 1);

                to_string (a1);
                printf (" || ");
                to_string (a2);
            } else if (strcmp (str1, "Not") == 0) {
                // | Not a -> "! " ^ braces (to_string a)
                //
                // let braces s = "{ " ^ s ^ " ; }"

                printf ("! ");
                printf ("{ "); // Inline 'braces'
                to_string (obj2);
                printf (" ; }"); // Inline 'braces'
            } else if (strcmp (str1, "Semi") == 0) {
                //   | Semi (a1,a2) -> to_string a1 ^ " ; " ^ to_string a2
                assert (json_object_array_length (obj2) == 2);
                json_object* a1 = json_object_array_get_idx (obj2, 0);
                json_object* a2 = json_object_array_get_idx (obj2, 1);

                to_string (a1);
                printf (" ; ");
                to_string (a2);
            } else if (strcmp (str1, "If") == 0) {
                assert (json_object_array_length (obj2) == 3);

                string_of_if (json_object_array_get_idx (obj2, 0),
                              json_object_array_get_idx (obj2, 1),
                              json_object_array_get_idx (obj2, 2));
            } else if (strcmp (str1, "While") == 0) {
                /*
                    | While (Not t,b) ->
                       "until " ^ to_string t ^ "; do " ^ to_string b ^ "; done "
                    | While (t,b) ->
                       "while " ^ to_string t ^ "; do " ^ to_string b ^ "; done "
                */
                assert (json_object_array_length (obj2) == 2);
                json_object* obj2_0 = json_object_array_get_idx (obj2, 0);
                json_object* obj2_1 = json_object_array_get_idx (obj2, 1);

                if (JSONArrayStartsWithStr (obj2_0, "Not")) {
                    //  | While (Not t,b) ->
                    //    "until " ^ to_string t ^ "; do " ^ to_string b ^ "; done "
                    assert (json_object_array_length (obj2_0) == 2);

                    // obj2_0_0 is "Not"
                    json_object* t = json_object_array_get_idx (obj2_0, 1);
                    json_object* b = obj2_1;

                    printf ("until ");
                    to_string (t);
                    printf ("; do ");
                    to_string (b);
                    printf ("; done ");
                } else {
                    // | While (t,b) ->
                    //   "while " ^ to_string t ^ "; do " ^ to_string b ^ "; done "
                    json_object* t = obj2_0;
                    json_object* b = obj2_1;

                    printf ("while ");
                    to_string (t);
                    printf ("; do ");
                    to_string (b);
                    printf ("; done ");
                    // Trailing space can look ugly but we keep it for consistency
                    // with ast.ml (which probably has a good reason for it) for:
                    // - scripts/circular/sine.sh'
                    // - scripts/circular/incr.sh'
                    // - scripts/usecases/shellcheck/distrotest_funs.sh
                    // - scripts/intermediary/web-index_p2_1_funs.sh
                }
            } else if (strcmp (str1, "For") == 0) {
                /*
                    | For (_,a,body,var) ->
                       "for " ^ var ^ " in " ^ string_of_arg a ^ "; do " ^
                       to_string body ^ "; done"
                */

                assert (json_object_array_length (obj2) == 4);
                // First param is ignored
                json_object* a    = json_object_array_get_idx (obj2, 1);
                json_object* body = json_object_array_get_idx (obj2, 2);
                json_object* var  = json_object_array_get_idx (obj2, 3);

                assert (json_object_get_type (var) == json_type_string);

                printf ("for ");
                printf ("%s", json_object_get_string (var));
                printf (" in ");
                string_of_arg (a);
                printf ("; do ");
                to_string (body);
                printf ("; done");
            } else if (strcmp (str1, "Case") == 0) {
                /*
                    | Case (_,a,cs) ->
                       "case " ^ string_of_arg a ^ " in " ^
                       separated string_of_case cs ^ " esac"
                */

                assert (json_object_array_length (obj2) == 3);
                // First param is ignored
                json_object* a  = json_object_array_get_idx (obj2, 1);
                json_object* cs = json_object_array_get_idx (obj2, 2);

                printf ("case ");
                string_of_arg (a);
                printf (" in ");

                assert (json_object_get_type (cs) == json_type_array);
                for (int i = 0; i < json_object_array_length (cs); i++) {
                    if (i > 0) {
                        putchar (' '); // Inline 'separated'
                    }
                    string_of_case (json_object_array_get_idx (cs, i));
                }
                printf (" esac");
            } else if (strcmp (str1, "Defun") == 0) {
                // | Defun (_,name,body) -> name ^ "() {\n" ^ to_string body ^ "\n}"

                assert (json_object_array_length (obj2) == 3);
                // First param is ignored
                json_object* name = json_object_array_get_idx (obj2, 1);
                json_object* body = json_object_array_get_idx (obj2, 2);

                assert (json_object_get_type (name) == json_type_string);

                printf ("%s", json_object_get_string (name));
                printf ("() {\n");
                to_string (body);
                printf ("\n}");
            } else {
                debug_json (obj2);
                printf ("Type: %s\n", str1);
                assert (! "Unexpected case");
            }

            break;

        case json_type_null:
        case json_type_boolean:
        case json_type_double:
        case json_type_int:
        case json_type_object:
        case json_type_string:
            assert (! "Unexpected type at top-level of to_string\n");
            break;

        default:
            assert ("! Not a valid json_type");
            break;
    }
}


/*
    string_of_if c t e =
      "if " ^ to_string c ^
      "; then " ^ to_string t ^
      (match e with
       | Command (-1,[],[],[]) -> "; fi" (* one-armed if *)
       | If (c,t,e) -> "; el" ^ string_of_if c t e
       | _ -> "; else " ^ to_string e ^ "; fi")
*/
void string_of_if (json_object* c, json_object* t, json_object* e) {
    assert (c != NULL);
    assert (t != NULL);
    assert (e != NULL);

    printf ("if ");
    to_string (c);
    printf ("; then ");
    to_string (t);

    assert (json_object_get_type (e) == json_type_array);
    assert (json_object_array_length (e) >= 1);

    json_object* e1 = json_object_array_get_idx (e, 0);
    assert (e1 != NULL);
    assert (json_object_get_type (e1) == json_type_string);

    const char* e1_str = json_object_get_string (e1);

    // | Command (-1,[],[],[]) -> "; fi" (* one-armed if *)
    int done = FALSE;
    if (json_object_array_length (e) == 2) {
        json_object* e2 = json_object_array_get_idx (e, 1);

        if (   (strcmp (e1_str, "Command") == 0)
            && (json_object_get_type (e2) == json_type_array)

            && (json_object_get_type (json_object_array_get_idx (e2, 0)) == json_type_int)
            && (json_object_get_int (json_object_array_get_idx (e2, 0)) == -1)

            && isEmptyJSONArray (json_object_array_get_idx (e2, 1))
            && isEmptyJSONArray (json_object_array_get_idx (e2, 2))
            && isEmptyJSONArray (json_object_array_get_idx (e2, 3))) {
            printf ("; fi");
            done = TRUE;
        }
    }

    if (strcmp (e1_str, "If") == 0) {
        // | If (c,t,e) -> "; el" ^ string_of_if c t e
        assert (json_object_array_length (e) == 2);
        json_object* e2 = json_object_array_get_idx (e, 1);

        assert (isJSONArrayOfLength (e2, 3));

        // Rename c/t/e to avoid confusion (even though lexical scoping allows it)
        json_object* cI = json_object_array_get_idx (e2, 0);
        json_object* tI = json_object_array_get_idx (e2, 1);
        json_object* eI = json_object_array_get_idx (e2, 2);

        printf ("; el");
        string_of_if (cI, tI, eI);
    } else if (! done) {
        printf ("; else ");
        to_string (e);
        printf ("; fi");
    }
}


/*
    string_of_arg_char = function
      | E '\'' -> "\\'"
      | E '\"' -> "\\\""
      | E '(' -> "\\("
      | E ')' -> "\\)"
      | E '{' -> "\\{"
      | E '}' -> "\\}"
      | E '$' -> "\\$"
      | E '!' -> "\\!"
      | E '&' -> "\\&"
      | E '|' -> "\\|"
      | E ';' -> "\\;"
      | C c -> String.make 1 c
      | E c -> Char.escaped c
      | T None -> "~"
      | T (Some u) -> "~" ^ u
      | A a -> "$((" ^ string_of_arg a ^ "))"
      | V (Length,_,name,_) -> "${#" ^ name ^ "}"
      | V (vt,nul,name,a) ->
         "${" ^ name ^ (if nul then ":" else "") ^ string_of_var_type vt ^ string_of_arg a ^ "}"
      | Q a -> "\"" ^ string_of_arg a ^ "\""
      | B t -> "$(" ^ to_string t ^ ")"
*/
void string_of_arg_char (json_object* first, json_object* second) {
    assert (first != NULL);
    assert (second != NULL);

    assert (json_object_get_type (first) == json_type_string);
    const char* firstStr = json_object_get_string (first);
    if (strcmp (firstStr, "E") == 0) {
        assert (json_object_get_type (second) == json_type_int);

        int second_int = json_object_get_int (second);
        switch (second_int) {
            case '\'':
                printf ("\\'");
                break;
            case '"':
                printf ("\\\"");
                break;
            case '(':
                printf ("\\(");
                break;
            case ')':
                printf ("\\)");
                break;
            case '{':
                printf ("\\{");
                break;
            case '}':
                printf ("\\}");
                break;
            case '$':
                printf ("\\$");
                break;
            case '!':
                printf ("\\!");
                break;
            case '&':
                printf ("\\&");
                break;
            case '|':
                printf ("\\|");
                break;
            case ';':
                printf ("\\;");
                break;

             // | E c -> Char.escaped c
             //
             // "All characters outside the ASCII printable range (32..126)
             //  are escaped, as well as backslash, double-quote, and
             //  single-quote."
             //     -- https://ocaml.org/releases/4.07/htmlman/libref/Char.html
            case '\\': // Special case
                printf ("\\\\");
                break;
            case '\t':
                printf ("\\t");
                break;

            default:
                if (second_int < 32 || second_int > 126) {
                    putchar ('\\');
                    printf ("%d", second_int);
                } else {
                    printf ("%c", second_int);
                }

                break;
        }
    } else if (strcmp (firstStr, "C") == 0) {
        assert (json_object_get_type (second) == json_type_int);
        putchar (json_object_get_int (second));
    } else if (strcmp (firstStr, "T") == 0) {
        /*
            | T None -> "~"
            | T (Some u) -> "~" ^ u
        */

        if (json_object_get_type (second) == json_type_string) {
            assert (strcmp (json_object_get_string (second), "None") == 0); // Left beef

            putchar ('~');
        } else if (json_object_get_type (second) == json_type_array) {
            assert (json_object_array_length (second) == 2);

            json_object* Some = json_object_array_get_idx (second, 0);
            json_object* u    = json_object_array_get_idx (second, 1);

            if (! JSONMatchesStr (Some, "Some")) {
                assert (! "Was not None or Some");
                exit (1);
            }

            assert (json_object_get_type (u) == json_type_string);

            putchar ('~');
            printf ("%s", json_object_get_string (u));
        } else {
            debug_json (second);

            assert (! "Unexpected pattern for T");
        }
    } else if (strcmp (firstStr, "A") == 0) {
        // | A a -> "$((" ^ string_of_arg a ^ "))"
        assert (json_object_get_type (second) == json_type_array);
        json_object* a = second;

        printf ("$((");
        string_of_arg (a);
        printf ("))");
    } else if (strcmp (firstStr, "V") == 0) {
        assert (isJSONArrayOfLength (second, 4));

        json_object* vt   = json_object_array_get_idx (second, 0);
        json_object* nul  = json_object_array_get_idx (second, 1);
        json_object* name = json_object_array_get_idx (second, 2);
        json_object* a    = json_object_array_get_idx (second, 3);

        assert (json_object_get_type (vt) == json_type_string);
        assert (json_object_get_type (name) == json_type_string);

        if (strcmp (json_object_get_string (vt), "Length") == 0) {
            // | V (Length,_,name,_) -> "${#" ^ name ^ "}"

            printf ("${#");
            printf ("%s", json_object_get_string (name));
            printf ("}");

            // assert (! "TODO: This case has been implemented but never tested");
        } else {
            // | V (vt,nul,name,a) ->
            //    "${" ^ name ^ (if nul then ":" else "") ^ string_of_var_type vt ^ string_of_arg a ^ "}"

            printf ("${");
            printf ("%s", json_object_get_string (name));

            if (json_object_get_flex_boolean (nul)) {
                printf (":");
            }

            string_of_var_type (vt);
            string_of_arg (a);
            printf ("}");
        }
    } else if (strcmp (firstStr, "Q") == 0) {
        putchar ('"');
        string_of_arg (second);
        putchar ('"');
    } else if (strcmp (firstStr, "B") == 0) {
        printf ("$(");
        to_string (second);
        printf (")");
    } else {
        assert (! "Unexpected arg_char");
    }
}


/*
    string_of_arg = function
      | [] -> ""
      | c :: a -> string_of_arg_char c ^ string_of_arg a
*/
void string_of_arg (json_object* arg1) {
    assert (arg1 != NULL);

    assert (json_object_get_type (arg1) == json_type_array);
    for (int i = 0; i < json_object_array_length (arg1); i++) {
        json_object* ch = json_object_array_get_idx (arg1, i);

        assert (ch != NULL);
        assert (isJSONArrayOfLength (ch, 2));

        json_object* ch1 = json_object_array_get_idx (ch, 0);
        assert (ch1 != NULL);
        assert (json_object_get_type (ch1) == json_type_string);

        json_object* ch2 = json_object_array_get_idx (ch, 1);
        assert (ch2 != NULL);

        string_of_arg_char (ch1, ch2);
    }
}


// string_of_assign (v,a) = v ^ "=" ^ string_of_arg a
void string_of_assign (json_object* va) {
    assert (va != NULL);

    assert (isJSONArrayOfLength (va, 2));

    json_object* v = json_object_array_get_idx (va, 0);
    assert (v != NULL);
    assert (json_object_get_type (v) == json_type_string);

    printf ("%s", json_object_get_string (v));
    putchar ('=');

    json_object* a = json_object_array_get_idx (va, 1);
    assert (a != NULL);
    assert (json_object_get_type (a) == json_type_array);
    string_of_arg (a);
}


/*
    string_of_case c =
      let pats = List.map string_of_arg c.cpattern in
      intercalate "|" pats ^ ") " ^ to_string c.cbody ^ ";;"
*/
// This case is unique because it uses JSON key/value.
void string_of_case (json_object* c) {
    assert (json_object_get_type (c) == json_type_object);

    json_object* cpattern;
    json_object* cbody;

    assert (json_object_object_get_ex (c, "cpattern", &cpattern));
    assert (json_object_object_get_ex (c, "cbody", &cbody));

    assert (json_object_get_type (cpattern) == json_type_array);
    assert (json_object_get_type (cbody) == json_type_array);

    for (int i = 0; i < json_object_array_length (cpattern); i++) {
        if (i > 0) {
            putchar ('|'); // Inline 'intercalate'
        }
        string_of_arg (json_object_array_get_idx (cpattern, i));
    }

    printf (") ");
    to_string (cbody);
    printf (";;");
}


/*
    string_of_redir = function
      | File (To,fd,a)      -> show_unless 1 fd ^ ">" ^ string_of_arg a
      | File (Clobber,fd,a) -> show_unless 1 fd ^ ">|" ^ string_of_arg a
      | File (From,fd,a)    -> show_unless 0 fd ^ "<" ^ string_of_arg a
      | File (FromTo,fd,a)  -> show_unless 0 fd ^ "<>" ^ string_of_arg a
      | File (Append,fd,a)  -> show_unless 1 fd ^ ">>" ^ string_of_arg a
      | Dup (ToFD,fd,tgt)   -> show_unless 1 fd ^ ">&" ^ string_of_arg tgt
      | Dup (FromFD,fd,tgt) -> show_unless 0 fd ^ "<&" ^ string_of_arg tgt
      | Heredoc (t,fd,a) ->
         let heredoc = string_of_arg a in
         let marker = fresh_marker (lines heredoc) "EOF" in
         show_unless 0 fd ^ "<<" ^
         (if t = XHere then marker else "'" ^ marker ^ "'") ^ "\n" ^ heredoc ^ marker ^ "\n"
*/
void string_of_redir (json_object* r) {
    assert (isJSONArrayOfLength (r, 2));

    // File (To, fd, a)
    // -r1- --- r2 ----

    json_object* r1 = json_object_array_get_idx (r, 0);
    assert (r1 != NULL);
    assert (json_object_get_type (r1) == json_type_string);

    json_object* r2 = json_object_array_get_idx (r, 1);
    assert (r2 != NULL);
    assert (json_object_get_type (r2) == json_type_array);

    const char* r1_str = json_object_get_string (r1);
    if (   (strcmp (r1_str, "File") == 0)
        || (strcmp (r1_str, "Dup") == 0)
        || (strcmp (r1_str, "Heredoc") == 0)) {
        assert (json_object_array_length (r2) == 3);

        json_object* direction = json_object_array_get_idx (r2, 0);
        json_object* fd        = json_object_array_get_idx (r2, 1);
        json_object* a         = json_object_array_get_idx (r2, 2);

        assert (direction != NULL);
        assert (json_object_get_type (direction) == json_type_string);

        assert (fd != NULL);
        assert (json_object_get_type (fd) == json_type_int);

        assert (a != NULL);
        assert (json_object_get_type (a) == json_type_array);

        const char* direction_str = json_object_get_string (direction);

        if (strcmp (r1_str, "File") == 0) {
            if (strcmp (direction_str, "To") == 0) {
                show_unless (1, json_object_get_int (fd));
                printf (">");
                string_of_arg (a);
            } else if (strcmp (direction_str, "Clobber") == 0) {
                show_unless (1, json_object_get_int (fd));
                printf (">|");
                string_of_arg (a);
            } else if (strcmp (direction_str, "From") == 0) {
                show_unless (0, json_object_get_int (fd));
                printf ("<");
                string_of_arg (a);
            } else if (strcmp (direction_str, "FromTo") == 0) {
                show_unless (0, json_object_get_int (fd));
                printf ("<>");
                string_of_arg (a);
            } else if (strcmp (direction_str, "Append") == 0) {
                show_unless (1, json_object_get_int (fd));
                printf (">>");
                string_of_arg (a);
            } else {
                assert (! "Invalid File case");
            }
        } else if (strcmp (r1_str, "Dup") == 0) {
            if (strcmp (direction_str, "ToFD") == 0) {
                show_unless (1, json_object_get_int (fd));
                printf (">&");
                string_of_arg (a);
            } else if (strcmp (direction_str, "FromFD") == 0) {
                show_unless (0, json_object_get_int (fd));
                printf ("<&");
                string_of_arg (a);
            } else {
                debug_json (r1);
                debug_json (r2);

                assert (! "Invalid Dup case");
            }
        } else if (strcmp (r1_str, "Heredoc") == 0) {
            /*
                | Heredoc (t,fd,a) ->
                   let heredoc = string_of_arg a in
                   let marker = fresh_marker (lines heredoc) "EOF" in
                   show_unless 0 fd ^ "<<" ^
                   (if t = XHere then marker else "'" ^ marker ^ "'") ^ "\n" ^ heredoc ^ marker ^ "\n"

                let lines = Str.split (Str.regexp "[\n\r]+")

                // If one of the lines contains the marker "EOF", keep appending 'F' until it's a unique marker
                let rec fresh_marker ls s =
                  if List.mem s ls
                  then fresh_marker ls (s ^ (String.sub s (String.length s - 1) 1))
                  else s
            */

            show_unless (0, json_object_get_int (fd));
            printf ("<<");

            // TODO: implement fresh_marker. Kludge would be to simply use a long, obscure piece of text.
            char* marker = "EOF";

            if (strcmp (direction_str, "XHere") == 0) {
                printf ("%s", marker);
            } else {
                putchar ('\'');
                printf ("%s", marker);
                putchar ('\'');
            }
            printf ("\n");
            string_of_arg (a);
            printf ("%s", marker);
            printf ("\n");
        } else {
            assert (! "This case shouldn't happen");
        }
    } else {
        debug_json (r1);
        debug_json (r2);

        assert (! "Invalid redir case");
    }
}


/*
    string_of_redirs rs =
      let ss = List.map string_of_redir rs in
      (if List.length ss > 0 then " " else "") ^ intercalate " " ss
*/
void string_of_redirs (json_object* rs) {
    assert (rs != NULL);

    assert (json_object_get_type (rs) == json_type_array);
    for (int i = 0; i < json_object_array_length (rs); i++) {
        // Inlined 'intercalate " "' and
        // optimized out special case whitespace
        putchar (' ');
        string_of_redir (json_object_array_get_idx (rs, i));
    }
}
