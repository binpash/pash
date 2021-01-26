/*


        ast2a.c : Exports 'of_node' which, analogously to the function in
                  libdash's ast.ml, converts an ordinary shell script into
                  a C-approximation of the libdash AST.

                  The functions here intentionally closely match ast.ml;
                  relevant OCaml snippets are included as comments.


        TODO: serialize to JSON fix, memory leaks


*/


#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <string.h>

#include "shell.h"

#include "init.h"
#include "input.h"
#include "main.h"
#include "memalloc.h"
#include "nodes.h"
#include "parser.h"

#include "dash2.h"
#include "ast2a.h"


#include "json_object.h"

#include "arg_char.h"
#include "CharList.h"
#include "Stack.h"


#define TRUE  1
#define FALSE 0


//------------------------------------------------------------------------------------------


#define STACK_CTLVar 100
#define STACK_CTLAri 101
#define STACK_CTLQuo 102


//------------------------------------------------------------------------------------------


#define TYPE_T_COMMAND    0
#define TYPE_T_PIPE       1
#define TYPE_T_REDIR      2
#define TYPE_T_BACKGROUND 3
#define TYPE_T_SUBSHELL   4
#define TYPE_T_AND        5
#define TYPE_T_OR         6
#define TYPE_T_NOT        7
#define TYPE_T_SEMI       8
#define TYPE_T_IF         9
#define TYPE_T_WHILE      10
#define TYPE_T_FOR        11
#define TYPE_T_CASE       12
#define TYPE_T_DEFUN      13

const char* SERIALIZE_TYPE_T []
    = {
       "Command", "Pipe", "Redir", "Background", "Subshell",
       "And", "Or", "Not", "Semi",
       "If", "While", "For", "Case", "Defun"
      };


/*
    var_type =
       | Normal
       | Minus
       | Plus
       | Question
       | Assign
       | TrimR
       | TrimRMax
       | TrimL
       | TrimLMax
       | Length
*/
#define VAR_TYPE_NORMAL    0x0
#define VAR_TYPE_MINUS     0x2
#define VAR_TYPE_PLUS      0x3
#define VAR_TYPE_QUESTION  0x4
#define VAR_TYPE_ASSIGN    0x5
#define VAR_TYPE_TRIMR     0x6
#define VAR_TYPE_TRIMRMAX  0x7
#define VAR_TYPE_TRIML     0x8
#define VAR_TYPE_TRIMLMAX  0x9
#define VAR_TYPE_LENGTH    0xA

const char* SERIALIZE_VAR_TYPE []
    = {
       "Normal",   // 0x0
       "UNUSED",
       "Minus",    // 0x2
       "Plus",     // 0x3
       "Question", // 0x4
       "Assign",   // 0x5
       "TrimR",    // 0x6
       "TrimRMax", // 0x7
       "TrimL",    // 0x8
       "TrimLMax", // 0x9
       "Length"    // 0xA
      };


/*
   redirection =
   | File of redir_type * int * arg
   | Dup of dup_type * int * arg
   | Heredoc of heredoc_type * int * arg
*/
// Don't mix up with REDIR_TYPE_*!
#define REDIRECTION_TYPE_FILE    0x0
#define REDIRECTION_TYPE_DUP     0x1
#define REDIRECTION_TYPE_HEREDOC 0x2

const char* SERIALIZE_REDIRECTION_TYPE [] = {"File", "Dup", "Heredoc"};


// redir_type = To | Clobber | From | FromTo | Append
#define REDIR_TYPE_TO      0x0
#define REDIR_TYPE_CLOBBER 0x1
#define REDIR_TYPE_FROM    0x2
#define REDIR_TYPE_FROMTO  0x3
#define REDIR_TYPE_APPEND  0x4

const char* SERIALIZE_REDIR_TYPE [] = {"To", "Clobber", "From", "FromTo", "Append"};


// dup_type = ToFD | FromFD
#define DUP_TYPE_TOFD   0x0
#define DUP_TYPE_FROMFD 0x1

const char* SERIALIZE_DUP_TYPE [] = {"ToFD", "FromFD"};


// heredoc_type = Here | XHere (* for when in a quote... not sure when this comes up *)
#define HEREDOC_TYPE_HERE  0x0
#define HEREDOC_TYPE_XHERE 0x1

const char* SERIALIZE_HEREDOC_TYPE [] = {"Here", "XHere"};


typedef CharList arg_TYPE; // arg = arg_char list


//------------------------------------------------------------------------------------------
// SIMPLE LIST TYPES


struct t_list {
    struct t_TYPE* t;
    struct t_list* next;
};

struct assign_list {
    struct assign_TYPE* assign;
    struct assign_list* next;
};

struct redirectionList {
    struct redirection_TYPE* redir;
    struct redirectionList* next;
};

// args = arg list
// Note that 'arg_TYPE' is typedef'ed above as a CharList
struct args_TYPE {
    arg_TYPE arg;
    struct args_TYPE* next;
};

struct case_list {
    struct case_TYPE* casey;
    struct case_list* next;
};


//------------------------------------------------------------------------------------------


// assign = string * arg
struct assign_TYPE {
    char* string;
    arg_TYPE arg;
};


// | File of redir_type * int * arg
//
// Not to be mistaken with <stdio.h>'s FILE
struct file_TYPE {
    int redir_type;
    int fd;
    arg_TYPE a;
};


// | Dup of dup_type * int * arg
struct dup_TYPE {
    int dup_type;
    int fd;
    arg_TYPE tgt;
};


// | Heredoc of heredoc_type * int * arg
struct heredoc_TYPE {
    int heredoc_type;
    int fd;
    arg_TYPE a;
};


/*
   redirection =
   | File of redir_type * int * arg
   | Dup of dup_type * int * arg
   | Heredoc of heredoc_type * int * arg
*/
struct redirection_TYPE {
    int type;

    union {
        struct file_TYPE file;
        struct dup_TYPE dup;
        struct heredoc_TYPE heredoc;
    };
};


// case = { cpattern : arg list; cbody : t }
//
// Note: the hash table is useful only for JSON serialization/deserialization
// purposes; we don't store it that way internally.
struct case_TYPE {
    struct args_TYPE cpattern;
    struct t* cbody;
};


//------------------------------------------------------------------------------------------
// type t = ...


//  | Command of linno * assign list * args * redirection list (* assign, args, redir *)
struct Command_TYPE {
    int linno;
    struct assign_list* assign;
    struct args_TYPE* args;
    struct redirectionList* redirect;
};


//  | Pipe of bool * t list (* background?, commands *)
struct Pipe_TYPE {
    int background;
    struct t_list* spill;
};


//  | Redir of linno * t * redirection list
struct Redir_TYPE {
    int linno;
    struct t_TYPE* t;
    struct redirectionList* redirect;
};


//  | Background of linno * t * redirection list 
struct Background_TYPE {
    int linno;
    struct t_TYPE* t;
    struct redirectionList* redirect;
};


//  | Subshell of linno * t * redirection list
struct Subshell_TYPE {
    int linno;
    struct t_TYPE* t;
    struct redirectionList* redirect;
};


//  | And of t * t
struct And_TYPE {
    struct t_TYPE* left;
    struct t_TYPE* right;
};


//  | Or of t * t
struct Or_TYPE {
    struct t_TYPE* left;
    struct t_TYPE* right;
};


//  | Not of t
struct Not_TYPE {
    struct t_TYPE* t;
};


//  | Semi of t * t
struct Semi_TYPE {
    struct t_TYPE* left;
    struct t_TYPE* right;
};


//  | If of t * t * t (* cond, then, else *)
struct If_TYPE {
    struct t_TYPE* test;
    struct t_TYPE* ifpart;
    struct t_TYPE* elsepart;
};


//  | While of t * t (* test, body *) (* until encoded as a While . Not *)
struct While_TYPE {
    struct t_TYPE* test;
    struct t_TYPE* body;
};


//  | For of linno * arg * t * string (* args, body, var *)
struct For_TYPE {
    int linno;
    arg_TYPE arg;
    struct t_TYPE* body;
    char* var;
};


// | Case of linno * arg * case list
struct Case_TYPE {
    int linno;
    arg_TYPE arg;
    struct case_list* cases;
};


//  | Defun of linno * string * t (* name, body *)
struct Defun_TYPE {
    int linno;
    char* name;
    struct t_TYPE* body;
};


struct t_TYPE {
    int type;

    union {
        struct Command_TYPE Command;
        struct Pipe_TYPE Pipe;
        struct Redir_TYPE Redir;
        struct Background_TYPE Background;
        struct Subshell_TYPE Subshell;
        struct And_TYPE And;
        struct Or_TYPE Or;
        struct Not_TYPE Not;
        struct Semi_TYPE Semi;
        struct If_TYPE If;
        struct While_TYPE While;
        struct For_TYPE For;
        struct Case_TYPE Case;
        struct Defun_TYPE Defun;
    };
};


//----------------------------------------------------------------------------------------------------------------------

struct json_object* json_arg_char_TYPE (struct arg_char_TYPE* head);
struct json_object* json_assign_list (struct assign_list* assign);
struct json_object* json_arg_TYPE (arg_TYPE arg);
struct json_object* json_args_TYPE (struct args_TYPE* args);
struct json_object* json_redirectionList (struct redirectionList* redirect);
struct json_object* json_t_TYPE (struct t_TYPE* t);

void debug_arg_char_TYPE (struct arg_char_TYPE* head);
void debug_assign_list (struct assign_list* assign);
void debug_arg_TYPE (arg_TYPE arg);
void debug_args_TYPE (struct args_TYPE* args);
void debug_redirectionList (struct redirectionList* redirect);
void debug_t_TYPE (struct t_TYPE* t);
void pour_the_t (struct t_TYPE* t);

void both (union node* n);

int var_type (int vstype);

int needs_escaping (char c);

struct t_TYPE* of_node (union node* n);
void mk_file (struct redirection_TYPE* redirection, union node* n, int ty);
void mk_dup (struct redirection_TYPE* redirection, union node* n, int ty);
void mk_here (struct redirection_TYPE* redirection, union node* n, int ty);
struct redirectionList* redirs (union node* n);
arg_TYPE to_arg (struct narg* n);
arg_TYPE parse_arg (CharList s, struct nodelist** bqlist, Stack stack);
char* parse_tilde (CharList left, CharList s);
arg_TYPE arg_char_FUNC (struct arg_char_TYPE* c, CharList s, struct nodelist** bqlist, Stack stack);
char* reverseStr (char* str);
struct assign_TYPE* to_assign (CharList left, CharList right);
struct assign_list* to_assigns (union node* assign);
struct args_TYPE* to_args (union node* n);

//----------------------------------------------------------------------------------------------------------------------


struct json_object* json_arg_char_TYPE (struct arg_char_TYPE* head) {
    struct json_object* root = json_object_new_array ();
    assert (root != NULL);

    assert (head->type >= 0);
//    assert (head->type < sizeof (SERIALIZE_TYPE_ARG_CHAR) / sizeof (char*));
    json_object_array_add (root, json_object_new_string (SERIALIZE_TYPE_ARG_CHAR [head->type]));

    switch (head->type) {
        case TYPE_ARG_CHAR_C: {
            json_object_array_add (root, json_object_new_int (head->C.c));
        }
            break;

        case TYPE_ARG_CHAR_E: {
            json_object_array_add (root, json_object_new_int (head->E.c));
        }
            break;

        case TYPE_ARG_CHAR_T: {
            if (head->T.str == NULL) {
                json_object_array_add (root, json_object_new_string ("None"));
            } else {
                struct json_object* nest = json_object_new_array ();
                assert (nest != NULL);

                json_object_array_add (nest, json_object_new_string ("Some"));
                json_object_array_add (nest, json_object_new_string (head->T.str));

                json_object_array_add (root, nest);
            }
        }
            break;

        case TYPE_ARG_CHAR_A: {
            json_object_array_add (root, json_arg_TYPE (head->A.arg));
        }
            break;

        case TYPE_ARG_CHAR_V: {
            struct json_object* nest = json_object_new_array ();
            assert (nest != NULL);

            json_object_array_add (nest, json_object_new_string (SERIALIZE_VAR_TYPE [head->V.var_type]));
            json_object_array_add (nest, json_object_new_boolean (head->V.vsnul));
            json_object_array_add (nest, json_object_new_string (head->V.str));
            json_object_array_add (nest, json_arg_TYPE (head->V.arg));

            json_object_array_add (root, nest);
        }
            break;

        case TYPE_ARG_CHAR_Q: {
            json_object_array_add (root, json_arg_TYPE (head->Q.arg));
        }
            break;

        case TYPE_ARG_CHAR_B: {
            json_object_array_add (root, json_t_TYPE (head->B.t));
        }
            break;

        default:
            break;
    }

    return (root);
}


struct json_object* json_assign_list (struct assign_list* assigns) {
    struct json_object* root = json_object_new_array ();
    assert (root != NULL);

    while (assigns != NULL) {
        struct assign_TYPE* assign = assigns->assign;

        struct json_object* nest = json_object_new_array ();
        assert (nest != NULL);

        json_object_array_add (nest, json_object_new_string (assign->string));
        json_object_array_add (nest, json_arg_TYPE (assign->arg));

        json_object_array_add (root, nest);

        assigns = assigns->next;
    }

    return (root);
}


struct json_object* json_arg_TYPE (arg_TYPE arg) {
    struct json_object* root = json_object_new_array ();
    assert (root != NULL);

    // arg is a list of arg_char, but stored in CharList
    while (! isCharListEmpty (arg)) {
        struct arg_char_TYPE* arg_char = charListHead_arg_char (arg);
        json_object_array_add (root, json_arg_char_TYPE (arg_char));

        charListTail (arg);
    }

    return (root);
}


struct json_object* json_args_TYPE (struct args_TYPE* args) {
    struct json_object* root = json_object_new_array ();
    assert (root != NULL);

    while (args != NULL) {
        json_object_array_add (root, json_arg_TYPE (args->arg));

        args = args->next;
    }

    return (root);
}


struct json_object* json_redirectionListI (struct redirectionList* redirect,
                                           json_object* root) {
    struct json_object* nest = json_object_new_array ();
    assert (nest != NULL);

    if (redirect == NULL) {
    } else {
        // TODO: bounds check
        json_object_array_add (nest, json_object_new_string (SERIALIZE_REDIRECTION_TYPE [redirect->redir->type]));

        switch (redirect->redir->type) {
            case REDIRECTION_TYPE_FILE: {
                struct json_object* nest2 = json_object_new_array ();
                assert (nest2 != NULL);

                json_object_array_add (nest2, json_object_new_string (SERIALIZE_REDIR_TYPE [redirect->redir->file.redir_type]));
                json_object_array_add (nest2, json_object_new_int (redirect->redir->file.fd));
                json_object_array_add (nest2, json_arg_TYPE (redirect->redir->file.a));

                json_object_array_add (nest, nest2);
            }
                break;

            case REDIRECTION_TYPE_DUP: {
                assert (! "TODO");
            }
                break;

            case REDIRECTION_TYPE_HEREDOC: {
                assert (! "TODO");
            }
                break;

            default:
                assert (! "Invalid redirection type");
                break;
        }

        json_object_array_add (root, nest);
        json_redirectionListI (redirect->next, root);
    }

    return (root);
}


// TODO: serialize JSON
struct json_object* json_redirectionList (struct redirectionList* redirect) {
    struct json_object* root = json_object_new_array ();
    assert (root != NULL);

    return (json_redirectionListI (redirect, root));
}


struct json_object* json_t_TYPE (struct t_TYPE* t) {
    struct json_object* root = json_object_new_array ();
    assert (root != NULL);

    if (t == NULL) {
        return root;
    }

    assert (t->type >= 0);
    assert (t->type < sizeof (SERIALIZE_TYPE_T) / sizeof (char*));
    json_object_array_add (root, json_object_new_string (SERIALIZE_TYPE_T [t->type]));

    struct json_object* nest = json_object_new_array ();
    assert (nest != NULL);

    switch (t->type) {
        case TYPE_T_COMMAND: {
            json_object_array_add (nest, json_object_new_int (t->Command.linno));
            json_object_array_add (nest, json_assign_list (t->Command.assign));
            json_object_array_add (nest, json_args_TYPE (t->Command.args));
            json_object_array_add (nest, json_redirectionList (t->Command.redirect));

//            json_object_array_add (root, json_object_new_string ("blah"));
//            json_object_array_add (root, json_object_new_string ("blah2"));
//            json_object_array_add (root, json_object_new_string ("blah3"));
            break;
        }

        case TYPE_T_PIPE: {
            json_object_array_add (nest, json_object_new_boolean (t->Pipe.background));

            struct t_list* spill = t->Pipe.spill;

            struct json_object* nest2 = json_object_new_array ();
            assert (nest2 != NULL);

            while (spill != NULL) {
                json_object_array_add (nest2, json_t_TYPE (spill->t));
                spill = spill->next;
            }

            json_object_array_add (nest, nest2);

            break;
        }

        case TYPE_T_REDIR: {
            json_object_array_add (nest, json_t_TYPE (t->Redir.t));
            json_object_array_add (nest, json_redirectionList (t->Redir.redirect));

            break;
        }

        case TYPE_T_BACKGROUND: {
            json_object_array_add (nest, json_t_TYPE (t->Background.t));
            json_object_array_add (nest, json_redirectionList (t->Background.redirect));

            break;
        }

        case TYPE_T_SUBSHELL: {
            json_object_array_add (nest, json_t_TYPE (t->Subshell.t));
            json_object_array_add (nest, json_redirectionList (t->Subshell.redirect));

            break;
        }

        case TYPE_T_AND: {
            json_object_array_add (nest, json_t_TYPE (t->And.left));
            json_object_array_add (nest, json_t_TYPE (t->And.right));
            break;
        }

        case TYPE_T_OR: {
            json_object_array_add (nest, json_t_TYPE (t->Or.left));
            json_object_array_add (nest, json_t_TYPE (t->Or.right));
            break;
        }

        case TYPE_T_NOT: {
            json_object_array_add (nest, json_t_TYPE (t->Not.t));
            break;
        }

        case TYPE_T_SEMI: {
            json_object_array_add (nest, json_t_TYPE (t->Or.left));
            json_object_array_add (nest, json_t_TYPE (t->Or.right));

            break;
        }

        case TYPE_T_IF: {
            json_object_array_add (nest, json_t_TYPE (t->If.test));
            json_object_array_add (nest, json_t_TYPE (t->If.ifpart));
            json_object_array_add (nest, json_t_TYPE (t->If.elsepart));

            break;
        }

        case TYPE_T_WHILE: {
            json_object_array_add (nest, json_t_TYPE (t->While.test));
            json_object_array_add (nest, json_t_TYPE (t->While.body));

            break;
        }

        case TYPE_T_FOR: {
            assert (! "debug_t_TYPE: TYPE_FOR not implemented");
            break;
        }

        case TYPE_T_CASE: {
            assert (! "debug_t_TYPE: TYPE_CASE not implemented");
            break;
        }

        case TYPE_T_DEFUN: {
            struct json_object* nest2 = json_object_new_array ();
            assert (nest2 != NULL);

            json_object_array_add (nest, json_object_new_int (t->Defun.linno));
            json_object_array_add (nest, json_object_new_string (t->Defun.name));

            json_object_array_add (nest2, json_t_TYPE (t->Defun.body));

            json_object_array_add (nest, nest2);
            break;
        }

        default:
            assert (! "Invalid t type");
            break;
    }

    json_object_array_add (root, nest);

    return (root);
}


//-----------------------------------------------------------------------------


void debug_arg_char_TYPE (struct arg_char_TYPE* head) {
    switch (head->type) {
        case TYPE_ARG_CHAR_C:
            putchar (head->C.c);
            break;
        case TYPE_ARG_CHAR_E:
            if (needs_escaping (head->E.c)) {
                putchar ('\\');
            }
            putchar (head->E.c);
            break;
        case TYPE_ARG_CHAR_T:
            putchar ('~');
            if (head->T.str != NULL) {
                printf ("%s", head->T.str);
            }
            break;
        case TYPE_ARG_CHAR_A:
            printf ("[a]");
            debug_arg_TYPE (head->A.arg);
            printf ("[/a]");
            break;
        case TYPE_ARG_CHAR_V:
            printf ("${%s}", head->V.str);
            debug_arg_TYPE (head->V.arg);
            break;
        case TYPE_ARG_CHAR_Q:
            putchar ('"');
            debug_arg_TYPE (head->Q.arg);
            putchar ('"');
            break;
        case TYPE_ARG_CHAR_B:
            printf ("$(");
            debug_t_TYPE (head->B.t);
            printf (")");
            break;
        default:
            break;
    }
}


void debug_assign_list (struct assign_list* assigns) {
    while (assigns != NULL) {
        struct assign_TYPE* assign = assigns->assign;

        printf ("%s=", assign->string);
        debug_arg_TYPE (assign->arg);

        assigns = assigns->next;
    }
}


void debug_arg_TYPE (arg_TYPE arg) {
    // arg is a list of arg_char, but stored in CharList
    while (! isCharListEmpty (arg)) {
        struct arg_char_TYPE* arg_char = charListHead_arg_char (arg);
        debug_arg_char_TYPE (arg_char);

        charListTail (arg);
    }
}


void debug_args_TYPE (struct args_TYPE* args) {
    while (args != NULL) {
        putchar (' ');
        debug_arg_TYPE (args->arg);

        args = args->next;
    }
}


void debug_redirectionList (struct redirectionList* redirect) {
    if (redirect == NULL) {
    } else {
        switch (redirect->redir->type) {
            case REDIRECTION_TYPE_FILE: {
                printf (" [%s, %d, ", SERIALIZE_REDIR_TYPE [redirect->redir->file.redir_type],
                                  redirect->redir->file.fd);
                debug_arg_TYPE (redirect->redir->file.a);
                putchar (']');
            }
                break;

            case REDIRECTION_TYPE_DUP: {
                assert (! "TODO");
            }
                break;

            case REDIRECTION_TYPE_HEREDOC: {
                assert (! "TODO");
            }
                break;

            default:
                assert (! "Invalid redirection type");
                break;
        }

        debug_redirectionList (redirect->next);
    }
}


void debug_t_TYPE (struct t_TYPE* t) {
    if (t == NULL) {
        return;
    }

    switch (t->type) {
        case TYPE_T_COMMAND: {
            debug_assign_list (t->Command.assign);
            putchar (' ');
            debug_args_TYPE (t->Command.args);
            putchar (' ');
            debug_redirectionList (t->Command.redirect);

            break;
        }

        case TYPE_T_PIPE: {
            if (t->Pipe.background) {
                printf (" {");
            }

            struct t_list* spill = t->Pipe.spill;

            while (spill != NULL) {
                if (spill != t->Pipe.spill) {
                    printf (" | ");
                }

                debug_t_TYPE (spill->t);
                spill = spill->next;
            }

            if (t->Pipe.background) {
                printf (" & } ");
            }

            break;
        }

        case TYPE_T_REDIR: {
            debug_t_TYPE (t->Redir.t);
            debug_redirectionList (t->Redir.redirect);

            break;
        }

        case TYPE_T_BACKGROUND: {
            printf (" {");
            debug_t_TYPE (t->Background.t);
            debug_redirectionList (t->Background.redirect);
            printf (" & };");

            break;
        }

        case TYPE_T_SUBSHELL: {
            printf (" (");
            debug_t_TYPE (t->Subshell.t);
            debug_redirectionList (t->Subshell.redirect);
            printf (" ) ;");

            break;
        }

        case TYPE_T_AND: {
            debug_t_TYPE (t->And.left);
            printf (" && ");
            debug_t_TYPE (t->And.right);
            break;
        }

        case TYPE_T_OR: {
            debug_t_TYPE (t->Or.left);
            printf (" || ");
            debug_t_TYPE (t->Or.right);
            break;
        }

        case TYPE_T_NOT: {
            printf ("!");
            debug_t_TYPE (t->Not.t);
            break;
        }

        case TYPE_T_SEMI: {
            debug_t_TYPE (t->Or.left);
            printf ("; ");
            debug_t_TYPE (t->Or.right);
            break;
        }

        case TYPE_T_IF: {
            printf ("if ");
            debug_t_TYPE (t->If.test);
            printf ("; then ");
            debug_t_TYPE (t->If.ifpart);
            printf ("; ");
            printf ("else ");
            debug_t_TYPE (t->If.elsepart);
            printf ("; fi ");

            break;
        }

        case TYPE_T_WHILE: {
            printf ("while ");
            debug_t_TYPE (t->While.test);
            printf ("; do ");
            debug_t_TYPE (t->While.body);
            printf ("; done\n");

            break;
        }

        case TYPE_T_FOR: {
            assert (! "debug_t_TYPE: TYPE_FOR not implemented");
            break;
        }

        case TYPE_T_CASE: {
            assert (! "debug_t_TYPE: TYPE_CASE not implemented");
            break;
        }

        case TYPE_T_DEFUN: {
            printf ("%s() {\n", t->Defun.name);
            debug_t_TYPE (t->Defun.body);
            printf ("\n");
            printf ("}\n");

            break;
        }

        default:
            assert (! "Invalid t type");
            break;
    }
}


void pour_the_t (struct t_TYPE* t) {
//    debug_t_TYPE (t);
    struct json_object* root = json_t_TYPE (t);

    const char* text = json_object_to_json_string_ext (root, JSON_C_TO_STRING_PRETTY);
    printf ("%s\n", text);
}


// -------------------------------------------------------------------------------------


void both (union node* n) {
    struct t_TYPE* t = of_node (n);
    pour_the_t (t);
}



// -------------------------------------------------------------------------------------



int var_type (int vstype) {
    // We don't have algebraic data types (e.g., 0x0 -> `Normal); instead, we have
    // defined VAR_TYPE_ such that this is the identity function.

    return (vstype);
}


// -------------------------------------------------------------------------------------


// let special_chars : char list = explode "|&;<>()$`\\\"'"


// let needs_escaping c = List.mem c special_chars
int needs_escaping (char c) {

    // Would be more efficient as a lookup table
    switch (c) {
        case '|':
        case '&':
        case ';':
        case '<':
        case '>':
        case '(':
        case ')':
        case '$':
        case '`':
        case '\\':
        case '"':
        case '\'':
            return TRUE;
        default:
            return FALSE;
    }
}


// -------------------------------------------------------------------------------------


// OCaml code excerpts inlined
struct t_TYPE* of_node (union node* n) {
    // printf ("[of_node] %s\n", NODE_NAMES [n->type]);

    struct t_TYPE* t = malloc (sizeof (struct t_TYPE));
    assert (t != NULL);

    /*
        let skip = Command (-1,[],[],[])

        if nullptr n
        then skip
    */
    if (n == NULL) {
        t->type = TYPE_T_COMMAND;

        t->Command.linno = -1;
        t->Command.assign = NULL;
        t->Command.args = NULL;
        t->Command.redirect = NULL;
    } else {
        switch (n->type) {
            case NCMD: {
                (void) 0;

                // struct ncmd {
                //       int type;
                //       int linno;
                //       union node *assign;
                //       union node *args;
                //       union node *redirect;
                // };

                /*
                    let n = n @-> node_ncmd in
                    Command (getf n ncmd_linno,
                             to_assigns (getf n ncmd_assign),
                             to_args (getf n ncmd_args),
                             redirs (getf n ncmd_redirect))
                */

                t->type = TYPE_T_COMMAND;

                t->Command.linno = n->ncmd.linno;
                t->Command.assign = to_assigns (n->ncmd.assign);
                t->Command.args = to_args (n->ncmd.args);
                t->Command.redirect = redirs (n->ncmd.redirect);
            }
                break; // Should be unreachable but we leave it in for clarity

            case NPIPE: {
                // struct npipe {
                //       int type;
                //       int backgnd;
                //       struct nodelist *cmdlist;
                // };

                /*
                    let n = n @-> node_npipe in
                    Pipe (getf n npipe_backgnd <> 0,
                          List.map of_node (nodelist (getf n npipe_cmdlist)))
                */

                t->type = TYPE_T_PIPE;

                t->Pipe.background = (n->npipe.backgnd != 0);
                t->Pipe.spill = NULL;

                struct t_list* last = NULL;

                struct nodelist* cmdlist = n->npipe.cmdlist;
                while (cmdlist != NULL) {
                    struct t_list* newLast = malloc (sizeof (struct t_list));
                    assert (newLast != NULL);

                    newLast->t = of_node (cmdlist->n);
                    newLast->next = NULL;

                    if (last == NULL) {
                        t->Pipe.spill = newLast;
                        last = newLast;
                    } else {
                        last->next = newLast;
                        last = newLast;
                    }

                    cmdlist = cmdlist->next;
                }
            }
                break;

            case NREDIR: {
                // struct nredir {
                //       int type;
                //       int linno;
                //       union node *n;
                //       union node *redirect;
                // };


                /*
                    let (ty,fd,arg) = of_nredir n in Redir (ty,fd,arg)


                    and of_nredir (n : node union ptr) =
                      let n = n @-> node_nredir in
                      (getf n nredir_linno, of_node (getf n nredir_n), redirs (getf n nredir_redirect))
                */

                // Inline 'Redir' function
                t->type = TYPE_T_REDIR;
                t->Redir.linno    = n->nredir.linno;
                t->Redir.t        = of_node (n->nredir.n);
                t->Redir.redirect = redirs (n->nredir.redirect);

                assert (! "of_node: NREDIR not implemented");
            }
                break;

            case NBACKGND: {
                // let (ty,fd,arg) = of_nredir n in Background (ty,fd,arg)

                // Inline 'Redir' function
                t->type = TYPE_T_BACKGROUND;
                t->Background.linno    = n->nredir.linno;
                t->Background.t        = of_node (n->nredir.n);
                t->Background.redirect = redirs (n->nredir.redirect);
            }
                break;

            case NSUBSHELL: {
                // let (ty,fd,arg) = of_nredir n in Subshell (ty,fd,arg)

                // Inline 'Redir' function
               t->type = TYPE_T_SUBSHELL;
               t->Subshell.linno    = n->nredir.linno;
               t->Subshell.t        = of_node (n->nredir.n);
               t->Subshell.redirect = redirs (n->nredir.redirect);
            }
                break;

            case NAND: {
                // struct nbinary {
                //       int type;
                //       union node *ch1;
                //       union node *ch2;
                // };

                // let (l,r) = of_binary n in And (l,r)
                //
                // Manually inlined 'of_binary':
                //    of_binary (n : node union ptr) =
                //       let n = n @-> node_nbinary in
                //       (of_node (getf n nbinary_ch1), of_node (getf n nbinary_ch2))

                t->type = TYPE_T_AND;
                t->And.left  = of_node (n->nbinary.ch1);
                t->And.right = of_node (n->nbinary.ch2);
           }
                break;


            case NOR: {
                // let (l,r) = of_binary n in Or (l,r)

                t->type = TYPE_T_OR;

                t->Or.left = of_node (n->nbinary.ch1);
                t->Or.right = of_node (n->nbinary.ch2);

                return (t);
            }
                break;

            case NSEMI: {
                // let (l,r) = of_binary n in Semi (l,r)

                t->type = TYPE_T_SEMI;

                t->Semi.left = of_node (n->nbinary.ch1);
                t->Semi.right = of_node (n->nbinary.ch2);

                return (t);
            }
                break;

            case NIF: {
                /*
                    let n = n @-> node_nif in
                    If (of_node (getf n nif_test),
                        of_node (getf n nif_ifpart),
                        of_node (getf n nif_elsepart))
                */

                t->type = TYPE_T_IF;
                t->If.test     = of_node (n->nif.test);
                t->If.ifpart   = of_node (n->nif.ifpart);
                t->If.elsepart = of_node (n->nif.elsepart);

                // struct nif {
                //       int type;
                //       union node *test;
                //       union node *ifpart;
                //       union node *elsepart;
                // };
                return (t);
            }
                break;

            case NWHILE: {
                // let (t,b) = of_binary n in While (t,b)

                t->type = TYPE_T_WHILE;
                t->While.test = of_node (n->nbinary.ch1);
                t->While.body = of_node (n->nbinary.ch2);

                return (t);
            }
                break;

            case NUNTIL: {
                // let (t,b) = of_binary n in While (Not t,b)

                assert (! "of_node: NUNTIL not implemented");
            }
                break;

            case NFOR: {
                // struct nfor {
                //       int type;
                //       int linno;
                //       union node *args;
                //       union node *body;
                //       char *var;
                // };

                /*
                    let n = n @-> node_nfor in
                    For (getf n nfor_linno,
                         to_arg (getf n nfor_args @-> node_narg),
                         of_node (getf n nfor_body),
                         getf n nfor_var)
                */
                assert (! "of_node: NFOR not implemented");
            }
                break;

            case NCASE: {
                // struct ncase {
                //       int type;
                //       int linno;
                //       union node *expr;
                //       union node *cases;
                // };

                /*
                   let n = n @-> node_ncase in
                   Case (getf n ncase_linno,
                         to_arg (getf n ncase_expr @-> node_narg),
                         List.map
                           (fun (pattern,body) ->
                             { cpattern = to_args pattern;
                               cbody = of_node body})
                           (caselist (getf n ncase_cases)))
                */
                assert (! "of_node: NCASE not implemented");
            }
                break;

            case NDEFUN: {
                // struct ndefun {
                //       int type;
                //       int linno;
                //       char *text;
                //       union node *body;
                // };

                /*
                    let n = n @-> node_ndefun in
                    Defun (getf n ndefun_linno,
                           getf n ndefun_text,
                           of_node (getf n ndefun_body))
                */
                t->type = TYPE_T_DEFUN;

                t->Defun.linno = n->ndefun.linno;
                t->Defun.name = n->ndefun.text;
                t->Defun.body = of_node (n->ndefun.body);
            }
                break;

            case NNOT: {
                // Not (of_node (getf (n @-> node_nnot) nnot_com))

                struct t_TYPE* t = malloc (sizeof (struct t_TYPE));
                assert (t != NULL);
                t->type = TYPE_T_NOT;

                t->Not.t = of_node (n->nbinary.ch1);

                return (t);
            }
                break;

            default: {
                assert (! "of_node: unexpected node type");
            }
                break;
        }
    }

    return t;
}



/*
    of_nredir (n : node union ptr) =
      let n = n @-> node_nredir in
      (getf n nredir_linno, of_node (getf n nredir_n), redirs (getf n nredir_redirect))
*/
// Manually inlined


/*
        let mk_file ty =
          let n = n @-> node_nfile in
          File (ty,getf n nfile_fd,to_arg (getf n nfile_fname @-> node_narg)) in
*/
void mk_file (struct redirection_TYPE* redirection, union node* n, int ty) {
    assert (redirection != NULL);
    assert (n != NULL);

    redirection->type = REDIRECTION_TYPE_FILE;

    redirection->file.redir_type = ty;
    redirection->file.fd         = n->nfile.fd;
    redirection->file.a          = to_arg (&(n->nfile.fname->narg));
}


/*
        let mk_dup ty =
          let n = n @-> node_ndup in
          let vname = getf n ndup_vname in
          let tgt =
            if nullptr vname
            then let dupfd = getf n ndup_dupfd in
                 if dupfd = -1
                 then [C '-']
                 else List.map (fun c -> C c) (explode (string_of_int dupfd))
            else to_arg (vname @-> node_narg)
          in
          Dup (ty,getf n ndup_fd,tgt) in
*/
void mk_dup (struct redirection_TYPE* redirection, union node* n, int ty) {
    assert (redirection != NULL);
    assert (n != NULL);

/*
    union node* vname = n->ndup.vname;

    arg_TYPE tgt;
    if (vname == NULL) {
        int dupfd = n->ndup.dupfd;

        if (dupfd == -1) {
            tgt = newArgCharC ('-');
        } else {
            // List.map (fun c -> C c) (explode (string_of_int dupfd))
            assert (! "TODO");
        }
    }

    redirection->type = REDIRECTION_TYPE_DUP;
*/
    assert (! "TODO");
}


/*
        let mk_here ty =
          let n = n @-> node_nhere in
          Heredoc (ty,getf n nhere_fd,to_arg (getf n nhere_doc @-> node_narg)) in
*/
void mk_here (struct redirection_TYPE* redirection, union node* n, int ty) {
    assert (redirection != NULL);
    assert (n != NULL);

    redirection->type = REDIRECTION_TYPE_HEREDOC;

    redirection->heredoc.heredoc_type = ty;
    redirection->heredoc.fd           = n->nhere.fd;
    redirection->heredoc.a            = to_arg (&(n->nhere.doc->narg));
}


// TODO
/*
        let mk_file ty =
          let n = n @-> node_nfile in
          File (ty,getf n nfile_fd,to_arg (getf n nfile_fname @-> node_narg)) in
    ...

    redirs (n : node union ptr) =
      if nullptr n
      then []
      else

        let h = match n @-> node_type with
        ...

    in
    h :: redirs (getf (n @-> node_nfile) nfile_next)
*/
struct redirectionList* redirs (union node* n) {
    if (n == NULL) {
        return NULL;
    }

    struct redirection_TYPE* redirection = malloc (sizeof (struct redirection_TYPE));
    assert (redirection != NULL);

    struct redirectionList* rlist = malloc (sizeof (struct redirectionList));
    assert (rlist != NULL);
    rlist->redir = redirection;

    switch (n->type) {
        case NTO: {
            // (* NTO *)
            // | 16 -> mk_file To

            redirection->type = REDIRECTION_TYPE_FILE;
            mk_file (redirection, n, REDIR_TYPE_TO);
        }
            break;

        case NCLOBBER: {
            // (* NCLOBBER *)
            // | 17 -> mk_file Clobber
            assert (! "Not implemented");
        }
            break;

        case NFROM: {
            // (* NFROM *)
            // | 18 -> mk_file From
            assert (! "Not implemented");
        }
            break;

        case NFROMTO: {
            // (* NFROMTO *)
            // | 19 -> mk_file FromTo
            assert (! "Not implemented");
        }
            break;

        case NAPPEND: {
            // (* NAPPEND *)
            // | 20 -> mk_file Append
            assert (! "Not implemented");
        }
            break;

        case NTOFD: {
            // (* NTOFD *)
            // | 21 -> mk_dup ToFD
            assert (! "Not implemented");
        }
            break;

        case NFROMFD: {
            // (* NFROMFD *)
            // | 22 -> mk_dup FromFD
            assert (! "Not implemented");
        }
            break;

        case NHERE: {
            // (* NHERE quoted heredoc---no expansion)*)
            // | 23 -> mk_here Here
            assert (! "Not implemented");
        }
            break;

        case NXHERE: {
            // (* NXHERE unquoted heredoc (param/command/arith expansion) *)
            // | 24 -> mk_here XHere
            assert (! "Not implemented");
        }
            break;

        default: {
            // | nt -> failwith ("unexpected node_type in redirlist: " ^ string_of_int nt)
            assert (! "Invalid redirs type");
        }
            break;
    }

    rlist->next = redirs (n->nfile.next);

    return rlist;
}


/*
    of_binary (n : node union ptr) =
      let n = n @-> node_nbinary in
      (of_node (getf n nbinary_ch1), of_node (getf n nbinary_ch2))
*/
// MANUALLY INLINED
// struct t_TYPE* of_binary (union node* n) {


/*
    to_arg (n : narg structure) : arg =
      let a,s,bqlist,stack = parse_arg (explode (getf n narg_text)) (getf n narg_backquote) [] in
      (* we should have used up the string and have no backquotes left in our list *)
      assert (s = []);
      assert (nullptr bqlist);
      assert (stack = []);
      a  
*/
arg_TYPE to_arg (struct narg* n) {
    CharList s = explode (n->text);
    struct nodelist* bqlist = n->backquote;
    Stack stack = newStack ();

    CharList a = parse_arg (s, &bqlist, stack);

    assert (isCharListEmpty (s));
    assert (bqlist == NULL);
    assert (isStackEmpty (stack));

    return (a);
}


/*
    parse_arg (s : char list) (bqlist : nodelist structure ptr) stack =
       match s,stack with
       ...
*/
arg_TYPE parse_arg (CharList s, struct nodelist** bqlist, Stack stack) {
    if (isCharListEmpty (s) && isStackEmpty (stack)) {
        // | [],[] -> [],[],bqlist,[]
        return newCharList ();
    } else if (isCharListEmpty (s) && topStack (stack) == STACK_CTLVar) {
        // | [],`CTLVar::_ -> failwith "End of string before CTLENDVAR"
        assert (! "End of string before CTLENDVAR");
    } else if (isCharListEmpty (s) && topStack (stack) == STACK_CTLAri) {
        // | [],`CTLAri::_ -> failwith "End of string before CTLENDARI"
        assert (! "End of string before CTLENDARI");
    } else if (isCharListEmpty (s) && topStack (stack) == STACK_CTLQuo) {
        // | [],`CTLQuo::_ -> failwith "End of string before CTLQUOTEMARK"
        assert (! "End of string before CTLENDQUOTEMARK");
    } else if ((charListLength (s) > 1) && (charListHead_char (s) == CTLESC)) {
        // | '\129'::c::s,_ -> arg_char (E c) s bqlist stack
        charListTail (s);

        char c = charListHead_char (s);
        charListTail (s);

        return arg_char_FUNC (newArgCharE (c), s, bqlist, stack);
    } else if ((charListLength (s) > 1) && (charListHead_char (s) == CTLVAR)) {
        /*
            Lightly reformatted

            | '\130'::t::s,_ ->
               let var_name,s = split_at (fun c -> c = '=') s in
               let t = int_of_char t in

               let v,s,bqlist,stack
               = match t land 0x0f, s with
                  (* VSNORMAL and VSLENGTH get special treatment
                  neither ever gets VSNUL
                  VSNORMAL is terminated just with the =, without a CTLENDVAR *)
                  (* VSNORMAL *)
                  | 0x1,'='::s ->
                     V (Normal,false,implode var_name,[]),s,bqlist,stack
                  (* VSLENGTH *)
                  | 0xa,'='::'\131'::s ->
                     V (Length,false,implode var_name,[]),s,bqlist,stack
                  | 0x1,c::_ | 0xa,c::_ ->
                     failwith ("Missing CTLENDVAR for VSNORMAL/VSLENGTH, found " ^ Char.escaped c)
                  (* every other VSTYPE takes mods before CTLENDVAR *)
                  | vstype,'='::s ->
                     let a,s,bqlist,stack' = parse_arg s bqlist (`CTLVar::stack) in
                     V (var_type vstype,t land 0x10 = 0x10,implode var_name,a), s, bqlist, stack'
                  | _,c::_ -> failwith ("Expected '=' terminating variable name, found " ^ Char.escaped c)
                  | _,[] -> failwith "Expected '=' terminating variable name, found EOF"
               in

               arg_char v s bqlist stack
        */

        charListTail (s);
        char t = (charListHead_char (s) & 0x0f); // Inline 'land 0x0f'
        charListTail (s);

        CharList var_name = split_at (s, '=');

        struct arg_char_TYPE* v = NULL;

        if (   (t == 0x1)
            && (charListLength (s) >= 1)
            && (charListHead_char (s) == '=')) {
            charListTail (s);

            v = newArgCharV (VAR_TYPE_NORMAL, FALSE, reverseStr (implode (var_name)), newCharList ());
        } else if (   (t == 0xA)
                   && (charListLength (s) >= 2)
                   && (charListHead_char (s) == '=')
                   && (charListSecond_char (s) == CTLENDVAR)) {
            charListTail (s);
            charListTail (s); // Drop two

            v = newArgCharV (VAR_TYPE_LENGTH, FALSE, reverseStr (implode (var_name)), newCharList ());
        } else if ((t == 0x1 || t == 0xA) && (charListLength (s) >= 0)) {
            assert (! "Missing CTLENDVAR for VSNORMAL/VSLENGTH");
        } else if (   (charListLength (s) >= 1)
                   && (charListHead_char (s) == '=')) {
            /*
                (* every other VSTYPE takes mods before CTLENDVAR *)
                | vstype,'='::s ->
                   let a,s,bqlist,stack' = parse_arg s bqlist (`CTLVar::stack) in
                   V (var_type vstype,t land 0x10 = 0x10,implode var_name,a), s, bqlist, stack'
            */

            charListTail (s);
            char vstype = t;

            pushStack (stack, STACK_CTLVar);
            arg_TYPE a = parse_arg (s, bqlist, stack);

            v = newArgCharV (var_type (vstype), (t == 0x10), reverseStr (implode (var_name)), a);

            assert (! "Implemented but not tested: 1");
        } else if (charListLength (s) >= 1) {
            assert (! "Expected '=' terminating variable name");
        } else {
            assert ("Expected '=' terminating variable name, found EOF");
        }

        destroyCharList (var_name);

        return (arg_char_FUNC (v, s, bqlist, stack));
    } else if (   (! isCharListEmpty (s))
               && (charListHead_char (s) == CTLENDVAR)) {
        /*
            (* CTLENDVAR *)
            | '\131'::s,`CTLVar::stack' -> [],s,bqlist,stack'
            | '\131'::_,`CTLAri::_ -> failwith "Saw CTLENDVAR before CTLENDARI"
            | '\131'::_,`CTLQuo::_ -> failwith "Saw CTLENDVAR before CTLQUOTEMARK"
            | '\131'::_,[] -> failwith "Saw CTLENDVAR outside of CTLVAR"
        */
        if (isStackEmpty (stack)) {
            assert (! "Saw CTLENDVAR outside of CTLVAR");
        } else if (topStack (stack) == STACK_CTLVar) {
            charListTail (s);
            popStack (stack);

            return newCharList ();
        } else if (topStack (stack) == STACK_CTLAri) {
            assert (! "Saw CTLENDVAR before CTLENDARI");
        } else if (topStack (stack) == STACK_CTLQuo) {
            assert (! "Saw CTLENDVAR before CTLQUOTEMARK");
        } else {
            assert (! "Unexpected stack contents");
        }
    } else if (   (! isCharListEmpty (s))
               && (charListHead_char (s) == CTLBACKQ)) {
        /*
            (* CTLBACKQ *)
            | '\132'::s,_ ->
               if nullptr bqlist
               then failwith "Saw CTLBACKQ but bqlist was null"
               else arg_char (B (of_node (bqlist @-> nodelist_n))) s (bqlist @-> nodelist_next) stack
        */

        charListTail (s);

        if (bqlist == NULL) {
            assert (! "Saw CTLBACKQ but bqlist was null");
        } else {
            struct arg_char_TYPE* a = newArgCharB (of_node ((*bqlist)->n));
            *bqlist = (*bqlist)->next; // TODO: check if this should persist

            return (arg_char_FUNC (a, s, bqlist, stack));
        }
    } else if ((! isCharListEmpty (s)) && (charListHead_char (s) == CTLARI)) {
        /*
            (* CTLARI *)
            | '\134'::s,_ ->
               let a,s,bqlist,stack' = parse_arg s bqlist (`CTLAri::stack) in
               assert (stack = stack');
               arg_char (A a) s bqlist stack'
        */

        charListTail (s);

        char* oldStackStr = serializeStack (stack);

        pushStack (stack, STACK_CTLAri);

        struct arg_char_TYPE* a = newArgCharA (parse_arg (s, bqlist, stack));

        char* newStackStr = serializeStack (stack);
        assert (strcmp (oldStackStr, newStackStr) == 0);
        free (oldStackStr);
        free (newStackStr);

        return (arg_char_FUNC (a, s, bqlist, stack));
    } else if (   (! isCharListEmpty (s))
               && (charListHead_char (s) == CTLENDARI)) {
        /*
            | '\135'::s,`CTLAri::stack' -> [],s,bqlist,stack'
            | '\135'::_,`CTLVar::_' -> failwith "Saw CTLENDARI before CTLENDVAR"
            | '\135'::_,`CTLQuo::_' -> failwith "Saw CTLENDARI before CTLQUOTEMARK"
            | '\135'::_,[] -> failwith "Saw CTLENDARI outside of CTLARI"
        */

        if (isStackEmpty (stack)) {
            assert (! "Saw CTLENDARI outside of CTLARI");
        } else if (topStack (stack) == STACK_CTLAri) {
            charListTail (s);
            popStack (stack);

            return newCharList ();
        } else if (topStack (stack) == STACK_CTLVar) {
            assert (! "Saw CTLENDARI before CTLENDVAR");
        } else if (topStack (stack) == STACK_CTLQuo) {
            assert (! "Saw CTLENDARI before CTLQUOTEMARK");
        } else {
            assert (! "Unexpected stack contents");
        }
    } else if (   (! isCharListEmpty (s))
               && (charListHead_char (s) == CTLQUOTEMARK)) {
        /*
            (* CTLQUOTEMARK *)
            | '\136'::s,`CTLQuo::stack' -> [],s,bqlist,stack'
            | '\136'::s,_ ->
               let a,s,bqlist,stack' = parse_arg s bqlist (`CTLQuo::stack) in
               assert (stack' = stack);
               arg_char (Q a) s bqlist stack'
        */

        charListTail (s);

        if ((! isStackEmpty (stack)) && (topStack (stack) == STACK_CTLQuo)) {
            popStack (stack);

            return newCharList ();
        } else {
            char* oldStackStr = serializeStack (stack);

            pushStack (stack, STACK_CTLQuo);

            struct arg_char_TYPE* a = newArgCharQ (parse_arg (s, bqlist, stack));

            char* newStackStr = serializeStack (stack);
            assert (strcmp (oldStackStr, newStackStr) == 0);
            free (oldStackStr);
            free (newStackStr);

            return (arg_char_FUNC (a, s, bqlist, stack));
        }
    } else if ((! isCharListEmpty (s)) && (charListHead_char (s) == '~')) {
        /*
           (* tildes *)
           | '~'::s,stack ->
              if List.exists (fun m -> m = `CTLQuo || m = `CTLAri) stack
              then (* we're in arithmetic or double quotes, so tilde is ignored *)
                 arg_char (C '~') s bqlist stack
              else
                 let uname,s' = parse_tilde [] s in
                 arg_char (T uname) s' bqlist stack
        */

        charListTail (s);
        if (   existsInStack (stack, STACK_CTLQuo)
            || existsInStack (stack, STACK_CTLAri)) {
            // N.B. Would be more efficient to search for both simultaneously

            return (arg_char_FUNC (newArgCharC ('~'), s, bqlist, stack));
        } else {
            char* uname = parse_tilde (newCharList (), s);

            return (arg_char_FUNC (newArgCharT (uname), s, bqlist, stack));
        }

    } else {
        /*
            (* ordinary character *)
            | c::s,_ ->
               arg_char (C c) s bqlist stack
        */
        char c = charListHead_char (s);

        charListTail (s);

        return (arg_char_FUNC (newArgCharC (c), s, bqlist, stack));
    }
}


static char* implodeOrNull (CharList acc, int reverse) {
    if (isCharListEmpty (acc)) {
        return NULL;
    } else if (reverse) {
        return reverseStr (implode (acc));
    } else {
        return implode (acc);
    }
}


/*
    parse_tilde acc =
      let ret = if acc = [] then None else Some (implode acc) in
      function
      ...
*/
char* parse_tilde (CharList acc, CharList s) {
    // OCaml has lazy evaluation but C does not, hence we can't afford
    // to define 'ret' here

    if (isCharListEmpty (s)) {
        // |     [] -> (ret , [])
        return implodeOrNull (acc, TRUE);
    } else if (charListHead_char (s) == CTLESC) {
        // (* CTLESC *)
        // |     '\129'::_ as s -> None, s
        return NULL;
    } else if (charListHead_char (s) == CTLQUOTEMARK) {
        // (* CTLQUOTEMARK *)
        // |     '\136'::_ as s -> None, s
        return NULL;
    } else if (   (charListHead_char (s) == CTLENDVAR)
               || (charListHead_char (s) == ':')
               || (charListHead_char (s) == '/')) {
        // (* terminal: CTLENDVAR, /, : *)
        // |     '\131'::_ as s -> ret, s
        // |     ':'::_ as s -> ret, s
        // |     '/'::_ as s -> ret, s
        return implodeOrNull (acc, TRUE);
    } else {
        // (* ordinary char *)
        // (* TODO 2019-01-03 only characters from the portable character set *)
        // |     c::s' -> parse_tilde (acc @ [c]) s'

        // This is reversed - implode will fix it
        prependCharList_char (acc, charListHead_char (s));
        charListTail (s);

        // Could convert tail recursion to a while loop
        return parse_tilde (acc, s);
    }
}


/*
    arg_char c s bqlist stack =
      let a,s,bqlist,stack = parse_arg s bqlist stack in
      (c::a,s,bqlist,stack)
*/
// Note that, in ast.ml, arg_char is both a type and a function!
arg_TYPE arg_char_FUNC (struct arg_char_TYPE* c, CharList s, struct nodelist** bqlist, Stack stack) {
    arg_TYPE a = parse_arg (s, bqlist, stack);

    if (c != NULL) {
        prependCharList_arg_char (a, (void*) c);
    }

    return (a);
}


// Helper kludge
char* reverseStr (char* str) {
    assert (str != NULL);

    int len = strlen (str);

    char* reverse = malloc (len + 1);
    assert (reverse != NULL);

    int i;
    for (i = 0; i < len; i++) {
        reverse [i] = str [len - 1 - i];
    }
    reverse [len] = '\0';

    return (reverse);
}


/*
    Lightly edited to invalid, but more readable OCaml syntax:

    to_assign v ca =
      | v   [] -> failwith ("Never found an '=' sign in assignment, got " ^ implode v)
      | v   (C '=') :: a    -> (implode v,a)
      | v   (C c  ) :: a    -> to_assign (v @ [c])      a
      | v   _               -> failwith "Unexpected special character in assignment"
*/
struct assign_TYPE* to_assign (CharList v, CharList ca) {
    if (isCharListEmpty (ca)) {
        assert (! "Never found an '=' sign in assignment");
    } else {
        struct arg_char_TYPE* ca_head = charListHead_arg_char (ca);
        charListTail (ca);

        if (ca_head->type != TYPE_ARG_CHAR_C) {
            assert (! "Unexpected special character in assignment");
        }

        if (ca_head->C.c == '=') {
            // Inline implode
            struct assign_TYPE* assign = malloc (sizeof (struct assign_TYPE));
            assert (assign != NULL);

            assign->string = reverseStr (implode (v));
            assign->arg    = ca;

            return (assign);
        } else {
            // Kludge: CharList.c only allows prepending, not appending;
            // implode will fix this.
            prependCharList_char (v, ca_head->C.c);

            return to_assign (v, ca);
        }
    }
}


/*
            // struct narg {
            //       int type;
            //       union node *next;
            //       char *text;
            //       struct nodelist *backquote;
            // };
*/
// to_assigns n = List.map (to_assign []) (to_args n)
struct assign_list* to_assigns (union node* n) {
    struct args_TYPE* args = to_args (n);
    struct assign_list* assigns = NULL;

    while (args != NULL) {
        struct assign_TYPE* assign = to_assign (newCharList (), args->arg);

        struct assign_list* newAssigns = malloc (sizeof (struct assign_list));
        assert (newAssigns != NULL);

        newAssigns->assign = assign;
        newAssigns->next   = assigns;

        assigns = newAssigns;

        args = args->next;
    }

    return assigns;
}


/*
    union node {
          int type;
          ...
          struct narg narg;
          ...
    };

    struct narg {
           int type;
           union node *next;
           char *text;
           struct nodelist *backquote;
    };


    to_args (n : node union ptr) : args =
      if nullptr n
      then []
      else (assert (n @-> node_type = 15);
            let n = n @-> node_narg in
            to_arg n::to_args (getf n narg_next))
*/
struct args_TYPE* to_args (union node* n) {
    if (n == NULL) {
        return NULL;
    } else {
        assert (n->type == NARG);

        struct args_TYPE* args = malloc (sizeof (struct args_TYPE));
        assert (args != NULL);

        args->arg  = to_arg (&(n->narg));
        args->next = to_args ((n->narg).next);

        return (args);
    }
}
