(* dash internals 

   call initialize before doing anything!
*)

val initialize : unit -> unit

(* stackmark discipline:

   (init_stack parse_next [process AST] pop_stack[deallocates dash AST])*

   see libdash/test/test.ml for an example usage in parse_all
*)
type stackmark
val init_stack : unit -> stackmark Ctypes.structure
val pop_stack : stackmark Ctypes.structure -> unit

val alloc_stack_string : string -> (char Ctypes.ptr)
val free_stack_string : (char Ctypes.ptr) -> unit
  
val popfile : unit -> unit
val setinputstring : (char Ctypes.ptr) -> unit
val setinputtostdin : unit -> unit
val setinputfile : ?push:bool -> string -> unit

val setvar : string -> string -> unit
val setalias : string -> string -> unit
val unalias : string -> unit

(* returns -1 when fd was closed; -2 on other errors *)
val freshfd_ge10 : int -> int 
  
(* Ctypes mappings of the node types *)
type node
val node : node Ctypes.union Ctypes.typ
val node_type : (int, node Ctypes.union) Ctypes.field

type nodelist
val nodelist_next :
  (nodelist Ctypes.structure Ctypes_static.ptr, nodelist Ctypes.structure)
  Ctypes.field
val nodelist_n :
  (node Ctypes.union Ctypes_static.ptr, nodelist Ctypes.structure)
  Ctypes.field

type ncmd
val ncmd : ncmd Ctypes.structure Ctypes.typ
val ncmd_type : (int, ncmd Ctypes.structure) Ctypes.field
val ncmd_linno : (int, ncmd Ctypes.structure) Ctypes.field
val ncmd_assign :
  (node Ctypes.union Ctypes_static.ptr, ncmd Ctypes.structure) Ctypes.field
val ncmd_args :
  (node Ctypes.union Ctypes_static.ptr, ncmd Ctypes.structure) Ctypes.field
val ncmd_redirect :
  (node Ctypes.union Ctypes_static.ptr, ncmd Ctypes.structure) Ctypes.field
val node_ncmd : (ncmd Ctypes.structure, node Ctypes.union) Ctypes.field

type npipe
val npipe : npipe Ctypes.structure Ctypes.typ
val npipe_type : (int, npipe Ctypes.structure) Ctypes.field
val npipe_backgnd : (int, npipe Ctypes.structure) Ctypes.field
val npipe_cmdlist :
  (nodelist Ctypes.structure Ctypes_static.ptr, npipe Ctypes.structure)
  Ctypes.field
val node_npipe : (npipe Ctypes.structure, node Ctypes.union) Ctypes.field

type nredir
val nredir : nredir Ctypes.structure Ctypes.typ
val nredir_type : (int, nredir Ctypes.structure) Ctypes.field
val nredir_linno : (int, nredir Ctypes.structure) Ctypes.field
val nredir_n :
  (node Ctypes.union Ctypes_static.ptr, nredir Ctypes.structure) Ctypes.field
val nredir_redirect :
  (node Ctypes.union Ctypes_static.ptr, nredir Ctypes.structure) Ctypes.field
val node_nredir : (nredir Ctypes.structure, node Ctypes.union) Ctypes.field

type nbinary
val nbinary : nbinary Ctypes.structure Ctypes.typ
val nbinary_type : (int, nbinary Ctypes.structure) Ctypes.field
val nbinary_ch1 :
  (node Ctypes.union Ctypes_static.ptr, nbinary Ctypes.structure)
  Ctypes.field
val nbinary_ch2 :
  (node Ctypes.union Ctypes_static.ptr, nbinary Ctypes.structure)
  Ctypes.field
val node_nbinary : (nbinary Ctypes.structure, node Ctypes.union) Ctypes.field

type nif
val nif : nif Ctypes.structure Ctypes.typ
val nif_type : (int, nif Ctypes.structure) Ctypes.field
val nif_test :
  (node Ctypes.union Ctypes_static.ptr, nif Ctypes.structure) Ctypes.field
val nif_ifpart :
  (node Ctypes.union Ctypes_static.ptr, nif Ctypes.structure) Ctypes.field
val nif_elsepart :
  (node Ctypes.union Ctypes_static.ptr, nif Ctypes.structure) Ctypes.field
val node_nif : (nif Ctypes.structure, node Ctypes.union) Ctypes.field

type nfor
val nfor : nfor Ctypes.structure Ctypes.typ
val nfor_type : (int, nfor Ctypes.structure) Ctypes.field
val nfor_linno : (int, nfor Ctypes.structure) Ctypes.field
val nfor_args :
  (node Ctypes.union Ctypes_static.ptr, nfor Ctypes.structure) Ctypes.field
val nfor_body :
  (node Ctypes.union Ctypes_static.ptr, nfor Ctypes.structure) Ctypes.field
val nfor_var : (string, nfor Ctypes.structure) Ctypes.field
val node_nfor : (nfor Ctypes.structure, node Ctypes.union) Ctypes.field

type ncase
val ncase : ncase Ctypes.structure Ctypes.typ
val ncase_type : (int, ncase Ctypes.structure) Ctypes.field
val ncase_linno : (int, ncase Ctypes.structure) Ctypes.field
val ncase_expr :
  (node Ctypes.union Ctypes_static.ptr, ncase Ctypes.structure) Ctypes.field
val ncase_cases :
  (node Ctypes.union Ctypes_static.ptr, ncase Ctypes.structure) Ctypes.field
val node_ncase : (ncase Ctypes.structure, node Ctypes.union) Ctypes.field

type nclist
val nclist : nclist Ctypes.structure Ctypes.typ
val nclist_type : (int, nclist Ctypes.structure) Ctypes.field
val nclist_next :
  (node Ctypes.union Ctypes_static.ptr, nclist Ctypes.structure) Ctypes.field
val nclist_pattern :
  (node Ctypes.union Ctypes_static.ptr, nclist Ctypes.structure) Ctypes.field
val nclist_body :
  (node Ctypes.union Ctypes_static.ptr, nclist Ctypes.structure) Ctypes.field
val node_nclist : (nclist Ctypes.structure, node Ctypes.union) Ctypes.field

type ndefun
val ndefun : ndefun Ctypes.structure Ctypes.typ
val ndefun_type : (int, ndefun Ctypes.structure) Ctypes.field
val ndefun_linno : (int, ndefun Ctypes.structure) Ctypes.field
val ndefun_text : (string, ndefun Ctypes.structure) Ctypes.field
val ndefun_body :
  (node Ctypes.union Ctypes_static.ptr, ndefun Ctypes.structure) Ctypes.field
val node_ndefun : (ndefun Ctypes.structure, node Ctypes.union) Ctypes.field

type narg
val narg : narg Ctypes.structure Ctypes.typ
val narg_type : (int, narg Ctypes.structure) Ctypes.field
val narg_next :
  (node Ctypes.union Ctypes_static.ptr, narg Ctypes.structure) Ctypes.field
val narg_text : (string, narg Ctypes.structure) Ctypes.field
val narg_backquote :
  (nodelist Ctypes.structure Ctypes_static.ptr, narg Ctypes.structure)
  Ctypes.field
val node_narg : (narg Ctypes.structure, node Ctypes.union) Ctypes.field

type nfile
val nfile : nfile Ctypes.structure Ctypes.typ
val nfile_type : (int, nfile Ctypes.structure) Ctypes.field
val nfile_next :
  (node Ctypes.union Ctypes_static.ptr, nfile Ctypes.structure) Ctypes.field
val nfile_fd : (int, nfile Ctypes.structure) Ctypes.field
val nfile_fname :
  (node Ctypes.union Ctypes_static.ptr, nfile Ctypes.structure) Ctypes.field
val nfile_expfname : (string, nfile Ctypes.structure) Ctypes.field
val node_nfile : (nfile Ctypes.structure, node Ctypes.union) Ctypes.field

type ndup
val ndup : ndup Ctypes.structure Ctypes.typ
val ndup_type : (int, ndup Ctypes.structure) Ctypes.field
val ndup_next :
  (node Ctypes.union Ctypes_static.ptr, ndup Ctypes.structure) Ctypes.field
val ndup_fd : (int, ndup Ctypes.structure) Ctypes.field
val ndup_dupfd : (int, ndup Ctypes.structure) Ctypes.field
val ndup_vname :
  (node Ctypes.union Ctypes_static.ptr, ndup Ctypes.structure) Ctypes.field
val node_ndup : (ndup Ctypes.structure, node Ctypes.union) Ctypes.field

type nhere
val nhere : nhere Ctypes.structure Ctypes.typ
val nhere_type : (int, nhere Ctypes.structure) Ctypes.field
val nhere_next :
  (node Ctypes.union Ctypes_static.ptr, nhere Ctypes.structure) Ctypes.field
val nhere_fd : (int, nhere Ctypes.structure) Ctypes.field
val nhere_doc :
  (node Ctypes.union Ctypes_static.ptr, nhere Ctypes.structure) Ctypes.field
val node_nhere : (nhere Ctypes.structure, node Ctypes.union) Ctypes.field

type nnot
val nnot : nnot Ctypes.structure Ctypes.typ
val nnot_type : (int, nnot Ctypes.structure) Ctypes.field
val nnot_com :
  (node Ctypes.union Ctypes_static.ptr, nnot Ctypes.structure) Ctypes.field
val node_nnot : (nnot Ctypes.structure, node Ctypes.union) Ctypes.field

val ( @-> ) :
  ('b, 'c) Ctypes.structured Ctypes.ptr ->
  ('a, ('b, 'c) Ctypes.structured) Ctypes.field -> 'a
val arglist : narg Ctypes.structure -> narg Ctypes.structure list
val nodelist :
  nodelist Ctypes.structure Ctypes.ptr -> node Ctypes.union Ctypes.ptr list
val redirlist :
  node Ctypes.union Ctypes.ptr ->
  [> `Dup of int * string * ndup Ctypes.structure
   | `File of int * string * nfile Ctypes.structure
   | `Here of int * string * bool * nhere Ctypes.structure ]
  list
val caselist :
  node Ctypes.union Ctypes.ptr ->
  (node Ctypes.union Ctypes_static.ptr * node Ctypes.union Ctypes_static.ptr)
  list

(* useful functions for working with the Ctypes AST *)
val addrof : 'a Ctypes.ptr -> nativeint
val eqptr : 'a Ctypes.ptr -> 'b Ctypes.ptr -> bool
val nullptr : 'a Ctypes.ptr -> bool

(* useful functions for pretty printing *)
val explode : string -> char list
val implode : char list -> string
val intercalate : string -> string list -> string
val lines : string -> string list
val split_at : ('a -> bool) -> 'a list -> 'a list * 'a list

(* shell-specific functions for pretty printing *)
val braces : string -> string
val parens : string -> string
val fresh_marker : string list -> string -> string

(* parser *)
type parse_result =
    Done
  | Error
  | Null
  | Parsed of node Ctypes.union Ctypes.ptr
val parse_next : ?interactive:bool -> unit -> parse_result

(* native pretty printer *)
val show : node Ctypes.union Ctypes.ptr -> string
