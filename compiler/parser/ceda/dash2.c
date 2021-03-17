#include <malloc.h>
#include <assert.h>

#include "shell.h"

#include "alias.h"
#include "init.h"
#include "input.h"
#include "main.h"
#include "memalloc.h"
#include "mystring.h"
#include "nodes.h"
#include "parser.h"
#include "redir.h"
#include "var.h"

#include "dash2.h"


// As a poor man's namespace, we prepend everything with "Dash_"
// Sometimes, Dash_x is the same as x (e.g., popfile), but sometimes
// they are subtly different (e.g., setvar), hence we always make a
// separate Dash_ shim.


struct stackmark* Dash_init_stack (void) {
    /* memalloc.h

       struct stackmark {
          struct stack_block *stackp;
          char *stacknxt;
          size_t stacknleft;
       };
    */

    struct stackmark* smark = malloc (sizeof (struct stackmark));
    assert (smark != NULL);

    setstackmark (smark);

    return (smark);
}


void Dash_pop_stack (struct stackmark* smark) {
    popstackmark (smark);
}



char* Dash_alloc_stack_string (char* str) {
    return sstrdup (str);
}


void Dash_free_stack_string (char* str) {
    stunalloc (str);
}


// Closely map to libdash/ocaml/dash.ml
void Dash_dash_init (void) {
    init ();
}


void Dash_initialize_dash_errno (void) {
    initialize_dash_errno ();
}


void Dash_initialize (void) {
    Dash_initialize_dash_errno ();

    Dash_dash_init ();
}


void Dash_popfile (void) {
    popfile ();
}


void Dash_setinputstring (char* str) {
    setinputstring (str);
}


void Dash_setinputtostdin (void) {
    setinputfd (0, 0); // fd, push
}


int Dash_setinputfile (char* str, int push) {
    return setinputfile (str, push ? 1 : 0);
}


void* Dash_setvar (const char* name, const char* val) {
    // struct var *setvar(const char *name, const char *val, int flags)
    return setvar (name, val, 0);
}


void Dash_setalias (const char *name, const char *val) {
    setalias (name, val);
}


void Dash_unalias (const char* name) {
    unalias (name);
}


int Dash_freshfd_ge10 (int fd) {
    return freshfd_ge10 (fd);
}


union node* Dash_parsecmd_safe (int interact) {
    // let parsecmd_safe : int -> node union ptr =
    //   foreign "parsecmd_safe" (int @-> returning (ptr node))

    return parsecmd_safe (interact);
}


/*
let parse_next ?interactive:(i=false) () =
  let n = parsecmd_safe (if i then 1 else 0) in
  if eqptr n neof
  then Done
  else if eqptr n nerr
  then Error
  else if nullptr n
  then Null (* comment or blank line or error ... *)
  else Parsed n
            
*/
union node* Dash_parse_next (int interactive) {
    union node* n = Dash_parsecmd_safe (interactive);

    // We use the parser.h types directly instead of using the OCaml types
    // e.g.,
    //     let nerr : node union ptr = foreign_value "lasttoken" node
    //     parser.h:#define NERR ((union node *)&lasttoken)

    return n;
}


/*
let explode s =
  let rec exp i l =
    if i < 0 then l else exp (i - 1) (s.[i] :: l) in
  exp (String.length s - 1) []
*/
CharList explode (char* str) {
    assert (str != NULL);

    CharList list = newCharList ();

    for (int i = strlen (str) - 1; i >= 0; i--) {
        prependCharList_char (list, str [i]);
    }

    return (list);
}


/*
let implode l =
  let s = Bytes.create (List.length l) in
  let rec imp i l =
    match l with
    | []  -> ()
    | (c::l) -> (Bytes.set s i c; imp (i+1) l)
  in
  imp 0 l;
  Bytes.unsafe_to_string s
*/
char* implode (CharList myList) {
    assert (myList != NULL);

    char* str = malloc (charListLength (myList) + 1);
    assert (str != NULL);

    int i = 0;
    while (! isCharListEmpty (myList)) {
        str [i] = charListHead_char (myList);
        charListTail (myList);

        i ++;
    }

    str [i] = '\0';

    destroyCharList (myList);

    return (str);
}


