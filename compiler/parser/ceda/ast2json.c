/*


        ast2json.c : Walks through our C-approximation of the libdash AST, and
                     uses JSON-C to generate a JSON representation that closely
                     matches the OCaml output (notably, including quirks such
                     as (None | Some str).


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
#include "ast2json.h"

#include "json_object.h"

#include "arg_char.h"
#include "CharList.h"
#include "Stack.h"


//----------------------------------------------------------------------------------------------------------------------

static struct json_object* json_arg_char_TYPE (struct arg_char_TYPE* head);
static struct json_object* json_assign_list (struct assign_list* assign);
static struct json_object* json_arg_TYPE (arg_TYPE arg);
static struct json_object* json_args_TYPE (struct args_TYPE* args);
static struct json_object* json_redirectionList (struct redirectionList* redirect);
static struct json_object* json_t_TYPE (struct t_TYPE* t);


//----------------------------------------------------------------------------------------------------------------------


static struct json_object* json_arg_char_TYPE (struct arg_char_TYPE* head) {
    struct json_object* root = json_object_new_array ();
    assert (root != NULL);

    assert (head->type >= 0);
    // TODO: bounds check
    // assert (head->type < sizeof (SERIALIZE_TYPE_ARG_CHAR) / sizeof (char*));
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


static struct json_object* json_assign_list (struct assign_list* assigns) {
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


static struct json_object* json_arg_TYPE (arg_TYPE arg) {
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


static struct json_object* json_args_TYPE (struct args_TYPE* args) {
    struct json_object* root = json_object_new_array ();
    assert (root != NULL);

    while (args != NULL) {
        json_object_array_add (root, json_arg_TYPE (args->arg));

        args = args->next;
    }

    return (root);
}


// Explicitly specify root so that the recursive step doesn't add extra nesting
static struct json_object* json_redirectionListI (struct redirectionList* redirect,
                                                  json_object* root) {
    struct json_object* nest = json_object_new_array ();
    assert (nest != NULL);

    if (redirect == NULL) {
    } else {
        // TODO: bounds check
        json_object_array_add (nest, json_object_new_string (SERIALIZE_REDIRECTION_TYPE [redirect->redir->type]));

        struct json_object* nest2 = json_object_new_array ();
        assert (nest2 != NULL);


        switch (redirect->redir->type) {
            case REDIRECTION_TYPE_FILE: {
                json_object_array_add (nest2, json_object_new_string (SERIALIZE_REDIR_TYPE [redirect->redir->file.redir_type]));
                json_object_array_add (nest2, json_object_new_int (redirect->redir->file.fd));
                json_object_array_add (nest2, json_arg_TYPE (redirect->redir->file.a));
            }
                break;

            case REDIRECTION_TYPE_DUP: {
                json_object_array_add (nest2, json_object_new_string (SERIALIZE_DUP_TYPE [redirect->redir->dup.dup_type]));
                json_object_array_add (nest2, json_object_new_int (redirect->redir->dup.fd));
                json_object_array_add (nest2, json_arg_TYPE (redirect->redir->dup.tgt));
            }
                break;

            case REDIRECTION_TYPE_HEREDOC: {
                json_object_array_add (nest2, json_object_new_string (SERIALIZE_HEREDOC_TYPE [redirect->redir->heredoc.heredoc_type]));
                json_object_array_add (nest2, json_object_new_int (redirect->redir->heredoc.fd));
                json_object_array_add (nest2, json_arg_TYPE (redirect->redir->heredoc.a));
            }
                break;

            default:
                assert (! "Invalid redirection type");
                break;
        }

        json_object_array_add (nest, nest2);
        json_object_array_add (root, nest);

        json_redirectionListI (redirect->next, root);
    }

    return (root);
}


static struct json_object* json_redirectionList (struct redirectionList* redirect) {
    struct json_object* root = json_object_new_array ();
    assert (root != NULL);

    return (json_redirectionListI (redirect, root));
}


static struct json_object* json_t_TYPE (struct t_TYPE* t) {
    struct json_object* root = json_object_new_array ();
    assert (root != NULL);

    if (t == NULL) {
        return root;
    }

    assert (t->type >= 0);
    // TODO: bounds check
    // assert (t->type < sizeof (SERIALIZE_TYPE_T) / sizeof (char*));
    json_object_array_add (root, json_object_new_string (SERIALIZE_TYPE_T [t->type]));

    struct json_object* nest = json_object_new_array ();
    assert (nest != NULL);

    switch (t->type) {
        case TYPE_T_COMMAND: {
            json_object_array_add (nest, json_object_new_int (t->Command.linno));
            json_object_array_add (nest, json_assign_list (t->Command.assign));
            json_object_array_add (nest, json_args_TYPE (t->Command.args));
            json_object_array_add (nest, json_redirectionList (t->Command.redirect));

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
            json_object_array_add (nest, json_object_new_int (t->Redir.linno));
            json_object_array_add (nest, json_t_TYPE (t->Redir.t));
            json_object_array_add (nest, json_redirectionList (t->Redir.redirect));

            break;
        }

        case TYPE_T_BACKGROUND: {
            json_object_array_add (nest, json_object_new_int (t->Background.linno));
            json_object_array_add (nest, json_t_TYPE (t->Background.t));
            json_object_array_add (nest, json_redirectionList (t->Background.redirect));

            break;
        }

        case TYPE_T_SUBSHELL: {
            json_object_array_add (nest, json_object_new_int (t->Background.linno));
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
            json_object_array_add (root /* not nested */, json_t_TYPE (t->Not.t));
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
            json_object_array_add (nest, json_object_new_int (t->For.linno));
            json_object_array_add (nest, json_arg_TYPE (t->For.arg));
            json_object_array_add (nest, json_t_TYPE (t->For.body));
            json_object_array_add (nest, json_object_new_string (t->For.var));

            break;
        }

        case TYPE_T_CASE: {
            json_object_array_add (nest, json_object_new_int (t->Case.linno));
            json_object_array_add (nest, json_arg_TYPE (t->Case.arg));

            struct json_object* nest2 = json_object_new_array ();
            assert (nest2 != NULL);

            struct case_list* cases_head = t->Case.cases;
            while (cases_head != NULL) {
                struct json_object* nest3 = json_object_new_object ();
                assert (nest3 != NULL);

                struct case_TYPE* casey = cases_head->casey;

                struct json_object* nest4 = json_object_new_array ();
                assert (nest4 != NULL);

                struct args_TYPE* cpattern_head = casey->cpattern;
                while (cpattern_head != NULL) {
                    json_object_array_add (nest4, json_arg_TYPE (cpattern_head->arg));

                    cpattern_head = cpattern_head->next;
                }
                json_object_object_add_ex (nest3, "cpattern", nest4, 0);
                json_object_object_add_ex (nest3, "cbody", json_t_TYPE (casey->cbody), 0);

                json_object_array_add (nest2, nest3);
                cases_head = cases_head->next;
            }

            json_object_array_add (nest, nest2);

            break;
        }

        case TYPE_T_DEFUN: {
            json_object_array_add (nest, json_object_new_int (t->Defun.linno));
            json_object_array_add (nest, json_object_new_string (t->Defun.name));

            json_object_array_add (nest, json_t_TYPE (t->Defun.body));

            break;
        }

        default:
            assert (! "Invalid t type");
            break;
    }

    if (t->type != TYPE_T_NOT) {
        json_object_array_add (root, nest);
    }

    return (root);
}


void pour_the_t (struct t_TYPE* t) {
    struct json_object* root = json_t_TYPE (t);

    // const char* text = json_object_to_json_string_ext (root, JSON_C_TO_STRING_PRETTY);
    // const char* text = json_object_to_json_string_ext (root, JSON_C_TO_STRING_SPACED);
    const char* text = json_object_to_json_string_ext (root, JSON_C_TO_STRING_PLAIN);
    printf ("%s", text);
}
