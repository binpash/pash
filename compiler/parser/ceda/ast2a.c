/*


        ast2a.c : Exports 'of_node' which, analogously to the function in
                  libdash's ast.ml, converts an ordinary shell script into
                  a C-approximation of the libdash AST.

                  The functions here intentionally closely match ast.ml;
                  relevant OCaml snippets are included as comments.


        N.B. this leaks memory - don't use it for a persistent process.


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

int var_type (int vstype);

int needs_escaping (char c);

void mk_file (struct redirection_TYPE* redirection, union node* n, int ty);
void mk_dup (struct redirection_TYPE* redirection, union node* n, int ty);
void mk_here (struct redirection_TYPE* redirection, union node* n, int ty);
struct redirectionList* redirs (union node* n);
arg_TYPE to_arg (struct narg* n);
arg_TYPE parse_arg (CharList s, struct nodelist** bqlist, Stack stack);
char* parse_tilde (CharList left, CharList s);
arg_TYPE arg_char_FUNC (struct arg_char_TYPE* c, CharList s, struct nodelist** bqlist, Stack stack);
struct assign_TYPE* to_assign (CharList left, CharList right);
struct assign_list* to_assigns (union node* assign);
struct args_TYPE* to_args (union node* n);

//----------------------------------------------------------------------------------------------------------------------


const char* SERIALIZE_TYPE_T []
    = {
       "Command", "Pipe", "Redir", "Background", "Subshell",
       "And", "Or", "Not", "Semi",
       "If", "While", "For", "Case", "Defun"
      };

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

const char* SERIALIZE_REDIRECTION_TYPE [] = {"File", "Dup", "Heredoc"};

const char* SERIALIZE_REDIR_TYPE [] = {"To", "Clobber", "From", "FromTo", "Append"};

const char* SERIALIZE_DUP_TYPE [] = {"ToFD", "FromFD"};

const char* SERIALIZE_HEREDOC_TYPE [] = {"Here", "XHere"};


//----------------------------------------------------------------------------------------------------------------------


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
                /*
                    let n = n @-> node_ncmd in
                    Command (getf n ncmd_linno,
                             to_assigns (getf n ncmd_assign),
                             to_args (getf n ncmd_args),
                             redirs (getf n ncmd_redirect))
                */

                t->type = TYPE_T_COMMAND;
                t->Command.linno    = n->ncmd.linno;
                t->Command.assign   = to_assigns (n->ncmd.assign);
                t->Command.args     = to_args (n->ncmd.args);
                t->Command.redirect = redirs (n->ncmd.redirect);
            }
                break; // Should be unreachable but we leave it in for clarity

            case NPIPE: {
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
                t->Or.left  = of_node (n->nbinary.ch1);
                t->Or.right = of_node (n->nbinary.ch2);
            }
                break;

            case NSEMI: {
                // let (l,r) = of_binary n in Semi (l,r)

                t->type = TYPE_T_SEMI;
                t->Semi.left  = of_node (n->nbinary.ch1);
                t->Semi.right = of_node (n->nbinary.ch2);
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
            }
                break;

            case NWHILE: {
                // let (t,b) = of_binary n in While (t,b)

                t->type = TYPE_T_WHILE;
                t->While.test = of_node (n->nbinary.ch1);
                t->While.body = of_node (n->nbinary.ch2);
            }
                break;

            case NUNTIL: {
                // let (t,b) = of_binary n in While (Not t,b)

                struct t_TYPE* not_node = malloc (sizeof (struct t_TYPE));
                assert (not_node != NULL);
                not_node->type = TYPE_T_NOT;
                not_node->Not.t = of_node (n->nbinary.ch1);

                t->type = TYPE_T_WHILE;
                t->While.test = not_node;
                t->While.body = of_node (n->nbinary.ch2);
            }
                break;

            case NFOR: {
                /*
                    let n = n @-> node_nfor in
                    For (getf n nfor_linno,
                         to_arg (getf n nfor_args @-> node_narg),
                         of_node (getf n nfor_body),
                         getf n nfor_var)
                */

                t->type = TYPE_T_FOR;
                t->For.linno = n->nfor.linno;
                t->For.arg   = to_arg (&(n->nfor.args->narg));
                t->For.body  = of_node (n->nfor.body);
                t->For.var   = n->nfor.var;
            }
                break;

            case NCASE: {
               /*
                   let rec caselist (n : node union ptr) =
                      if nullptr n
                      then []
                      else
                         let n = n @-> node_nclist in
                            assert (getf n nclist_type = 13); (* NCLIST *)
                            (getf n nclist_pattern, getf n nclist_body)::caselist (getf n nclist_next)

                   let n = n @-> node_ncase in
                   Case (getf n ncase_linno,
                         to_arg (getf n ncase_expr @-> node_narg),
                         List.map
                           (fun (pattern,body) ->
                             { cpattern = to_args pattern;
                               cbody = of_node body})
                           (caselist (getf n ncase_cases)))
                */

                t->type = TYPE_T_CASE;
                t->Case.linno = n->ncase.linno;
                t->Case.arg   = to_arg (&(n->ncase.expr->narg));

                struct case_list* cases_ceda = NULL;
                struct case_list* cases_ceda_tail = cases_ceda;

                union node* case_head = n->ncase.cases;
                while (case_head != NULL) {
                    assert (case_head->type == NCLIST);

                    struct case_TYPE* newCase = malloc (sizeof (struct case_TYPE));
                    assert (newCase != NULL);

                    newCase->cpattern = to_args (case_head->nclist.pattern);
                    newCase->cbody    = of_node (case_head->nclist.body);

                    if (cases_ceda == NULL) {
                        cases_ceda = malloc (sizeof (struct case_list));
                        assert (cases_ceda != NULL);

                        cases_ceda_tail = cases_ceda;
                    } else {
                        cases_ceda_tail->next = malloc (sizeof (struct case_list));
                        assert (cases_ceda_tail->next != NULL);

                        cases_ceda_tail = cases_ceda_tail->next;
                    }

                    cases_ceda_tail->casey = newCase;
                    cases_ceda_tail->next = NULL;

                    case_head = case_head->nclist.next;
                }

                t->Case.cases = cases_ceda;
            }
                break;

            case NDEFUN: {
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
                t->type = TYPE_T_NOT;
                t->Not.t = of_node (n->nbinary.ch1);
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

    union node* vname = n->ndup.vname;

    arg_TYPE tgt;
    if (vname == NULL) {
        tgt = newCharList ();

        int dupfd = n->ndup.dupfd;

        if (dupfd == -1) {
            appendCharList_arg_char (tgt, newArgCharC ('-'));
        } else {
            // List.map (fun c -> C c) (explode (string_of_int dupfd))
            // Use string instead of list

            char dupfd_str [640];
            sprintf (dupfd_str, "%d", dupfd);

            int i = 0;
            while ((i < 100) && (dupfd_str [i] != '\0')) {
                appendCharList_arg_char (tgt, newArgCharC (dupfd_str [i]));
                i ++;
            }
        }
    } else {
        tgt = to_arg (&(vname->narg)); // Not used!
    }

    redirection->type = REDIRECTION_TYPE_DUP;

    redirection->dup.dup_type = ty;
    redirection->dup.fd       = n->ndup.fd;
    redirection->dup.tgt      = tgt;
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


/*
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
    struct redirectionList* headRL = NULL;
    struct redirectionList* lastRL = NULL;

    while (n != NULL) {
        struct redirection_TYPE* redirection = malloc (sizeof (struct redirection_TYPE));
        assert (redirection != NULL);

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

                redirection->type = REDIRECTION_TYPE_FILE;
                mk_file (redirection, n, REDIR_TYPE_CLOBBER);
            }
                break;

            case NFROM: {
                // (* NFROM *)
                // | 18 -> mk_file From
                redirection->type = REDIRECTION_TYPE_FILE;
                mk_file (redirection, n, REDIR_TYPE_FROM);
            }
                break;

            case NFROMTO: {
                // (* NFROMTO *)
                // | 19 -> mk_file FromTo
                redirection->type = REDIRECTION_TYPE_FILE;
                mk_file (redirection, n, REDIR_TYPE_FROMTO);
            }
                break;

            case NAPPEND: {
                // (* NAPPEND *)
                // | 20 -> mk_file Append
                redirection->type = REDIRECTION_TYPE_FILE;
                mk_file (redirection, n, REDIR_TYPE_APPEND);
            }
                break;

            case NTOFD: {
                // (* NTOFD *)
                // | 21 -> mk_dup ToFD

                redirection->type = REDIRECTION_TYPE_DUP;
                mk_dup (redirection, n, DUP_TYPE_TOFD);
            }
                break;

            case NFROMFD: {
                // (* NFROMFD *)
                // | 22 -> mk_dup FromFD

                redirection->type = REDIRECTION_TYPE_DUP;
                mk_dup (redirection, n, DUP_TYPE_FROMFD);
            }
                break;

            case NHERE: {
                // (* NHERE quoted heredoc---no expansion)*)
                // | 23 -> mk_here Here

                redirection->type = REDIRECTION_TYPE_HEREDOC;
                mk_here (redirection, n, HEREDOC_TYPE_HERE);
            }
                break;

            case NXHERE: {
                // (* NXHERE unquoted heredoc (param/command/arith expansion) *)
                // | 24 -> mk_here XHere
                redirection->type = REDIRECTION_TYPE_HEREDOC;
                mk_here (redirection, n, HEREDOC_TYPE_XHERE);
            }
                break;

            default: {
                // | nt -> failwith ("unexpected node_type in redirlist: " ^ string_of_int nt)
                assert (! "Invalid redirs type");
            }
                break;
        }

        struct redirectionList* newRL = malloc (sizeof (struct redirectionList));
        assert (newRL != NULL);
        newRL->redir = redirection;
        newRL->next = NULL;

        if (headRL == NULL) {
            headRL = newRL;
            lastRL = newRL;
        } else {
            lastRL->next = newRL;
            lastRL = newRL;
        }

        n = n->nfile.next;
    }

    return headRL;
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

    destroyStack (stack);

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
                  ...
               in

               arg_char v s bqlist stack
        */

        charListTail (s);
        char t = charListHead_char (s);
        char tM = t & 0xF;
        charListTail (s);

        CharList var_name = split_at (s, '=');

        struct arg_char_TYPE* v = NULL;

        if (   (tM == 0x1)
            && (charListLength (s) >= 1)
            && (charListHead_char (s) == '=')) {
            /*
                  (* VSNORMAL and VSLENGTH get special treatment
                  neither ever gets VSNUL
                  VSNORMAL is terminated just with the =, without a CTLENDVAR *)
                  (* VSNORMAL *)
                  | 0x1,'='::s ->
                     V (Normal,false,implode var_name,[]),s,bqlist,stack
            */

            charListTail (s);

            v = newArgCharV (VAR_TYPE_NORMAL, FALSE, implode (var_name), newCharList ());
        } else if (   (tM == 0xA)
                   && (charListLength (s) >= 2)
                   && (charListHead_char (s) == '=')
                   && (charListSecond_char (s) == CTLENDVAR)) {
            /*
                  (* VSLENGTH *)
                  | 0xa,'='::'\131'::s ->
                     V (Length,false,implode var_name,[]),s,bqlist,stack
            */

            charListTail (s);
            charListTail (s); // Drop two

            v = newArgCharV (VAR_TYPE_LENGTH, FALSE, implode (var_name), newCharList ());
        } else if ((tM == 0x1 || tM == 0xA) && (charListLength (s) >= 0)) {
            /*
                  | 0x1,c::_ | 0xa,c::_ ->
                     failwith ("Missing CTLENDVAR for VSNORMAL/VSLENGTH, found " ^ Char.escaped c)
            */

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
            char vstype = tM;

            pushStack (stack, STACK_CTLVar);
            arg_TYPE a = parse_arg (s, bqlist, stack);

            v = newArgCharV (var_type (vstype), ((t & 0x10) == 0x10), implode (var_name), a);
        } else if (charListLength (s) >= 1) {
            // | _,c::_ -> failwith ("Expected '=' terminating variable name, found " ^ Char.escaped c)
            assert (! "Expected '=' terminating variable name");
        } else {
            // | _,[] -> failwith "Expected '=' terminating variable name, found EOF"
            assert ("Expected '=' terminating variable name, found EOF");
        }

        // 'implode' has already destroyed the list
        // destroyCharList (var_name);

        return (arg_char_FUNC (v, s, bqlist, stack));
    } else if (   FALSE /* match Pash's version of libdash */
               && (! isCharListEmpty (s))
               && (charListHead_char (s) == CTLVAR)) {
        // | '\130'::_, _ -> raise (ParseException "bad substitution (missing variable name in ${}?")
        assert (! "bad substitution (missing variable name in ${}?");
    } else if (   (! isCharListEmpty (s))
               && (charListHead_char (s) == CTLENDVAR)) {
        if (isStackEmpty (stack)) {
            // | '\131'::_,[] -> failwith "Saw CTLENDVAR outside of CTLVAR"
            assert (! "Saw CTLENDVAR outside of CTLVAR");
        } else if (topStack (stack) == STACK_CTLVar) {
            // | '\131'::s,`CTLVar::stack' -> [],s,bqlist,stack'
            charListTail (s);
            popStack (stack);

            return newCharList ();
        } else if (topStack (stack) == STACK_CTLAri) {
            // | '\131'::_,`CTLAri::_ -> failwith "Saw CTLENDVAR before CTLENDARI"
            assert (! "Saw CTLENDVAR before CTLENDARI");
        } else if (topStack (stack) == STACK_CTLQuo) {
            // | '\131'::_,`CTLQuo::_ -> failwith "Saw CTLENDVAR before CTLQUOTEMARK"
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
            *bqlist = (*bqlist)->next;

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
        if (isStackEmpty (stack)) {
            // | '\135'::_,[] -> failwith "Saw CTLENDARI outside of CTLARI"
            assert (! "Saw CTLENDARI outside of CTLARI");
        } else if (topStack (stack) == STACK_CTLAri) {
            // | '\135'::s,`CTLAri::stack' -> [],s,bqlist,stack'
            charListTail (s);
            popStack (stack);

            return newCharList ();
        } else if (topStack (stack) == STACK_CTLVar) {
            // | '\135'::_,`CTLVar::_' -> failwith "Saw CTLENDARI before CTLENDVAR"
            assert (! "Saw CTLENDARI before CTLENDVAR");
        } else if (topStack (stack) == STACK_CTLQuo) {
            // | '\135'::_,`CTLQuo::_' -> failwith "Saw CTLENDARI before CTLQUOTEMARK"
            assert (! "Saw CTLENDARI before CTLQUOTEMARK");
        } else {
            assert (! "Unexpected stack contents");
        }
    } else if (   (! isCharListEmpty (s))
               && (charListHead_char (s) == CTLQUOTEMARK)) {
        charListTail (s); // See below

        if ((! isStackEmpty (stack)) && (topStack (stack) == STACK_CTLQuo)) {
            // | '\136'::s,`CTLQuo::stack' -> [],s,bqlist,stack'

            popStack (stack);

            return newCharList ();
        } else {
            /*
                (* CTLQUOTEMARK *)
                | '\136'::s,_ ->
                   let a,s,bqlist,stack' = parse_arg s bqlist (`CTLQuo::stack) in
                   assert (stack' = stack);
                   arg_char (Q a) s bqlist stack'
            */

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

    return NULL; // Reachable if NDEBUG
}


static char* implodeOrNull (CharList acc) {
    if (isCharListEmpty (acc)) {
        return NULL;
    } else {
        return implode (acc);
    }
}


// Note: we encode "None" as NULL, and "Some String" as char*
//       The OCaml-style serialization will be handled by ast2json.c.
/*
    parse_tilde acc =
      let ret = if acc = [] then None else Some (implode acc) in
      function
      ...
*/
char* parse_tilde (CharList acc, CharList s) {
    // C does not have lazy evaluation
    // hence we can't afford to define 'ret' here

    if (isCharListEmpty (s)) {
        // |     [] -> (ret , [])
        return implodeOrNull (acc);
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
        return implodeOrNull (acc);
    } else {
        // (* ordinary char *)
        // (* TODO 2019-01-03 only characters from the portable character set *)
        // |     c::s' -> parse_tilde (acc @ [c]) s'

        char c = charListHead_char (s);
        charListTail (s);
        appendCharList_char (acc, c);

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
        prependCharList_arg_char (a, (void*) c); // Intentional preprend
    }

    return (a);
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
        struct arg_char_TYPE* c = charListHead_arg_char (ca);
        charListTail (ca);
        CharList a = ca;

        if (c->type != TYPE_ARG_CHAR_C) {
            assert (! "Unexpected special character in assignment");
        }

        if (c->C.c == '=') {
            struct assign_TYPE* assign = malloc (sizeof (struct assign_TYPE));
            assert (assign != NULL);

            assign->string = implode (v);
            assign->arg    = a;

            return (assign);
        } else {
            appendCharList_char (v, c->C.c);

            return to_assign (v, a);
        }
    }

    return NULL; // Reachable if NDEBUG
}


// to_assigns n = List.map (to_assign []) (to_args n)
struct assign_list* to_assigns (union node* n) {
    struct args_TYPE* args = to_args (n);
    struct assign_list* assigns = NULL;

    struct assign_list* assignsLast = NULL;

    while (args != NULL) {
        struct assign_TYPE* assign = to_assign (newCharList (), args->arg);

        if (assigns == NULL) {
            assigns = malloc (sizeof (struct assign_list));
            assert (assigns != NULL);

            assignsLast = assigns;
        } else {
            assignsLast->next = malloc (sizeof (struct assign_list));
            assert (assignsLast->next != NULL);

            assignsLast = assignsLast->next;
        }

        assignsLast->assign = assign;
        assignsLast->next   = NULL;

        struct args_TYPE* next = args->next;
        free (args);
        args = next;
    }

    free (args);

    return assigns;
}


/*
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
