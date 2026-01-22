type linno = int

type t =
    Command of (linno * assign list * args * redirection list)
  | Pipe of (bool * t list)
  | Redir of (linno * t * redirection list)
  | Background of (linno * t * redirection list)
  | Subshell of (linno * t * redirection list)
  | And of (t * t)
  | Or of (t * t)
  | Not of t
  | Semi of (t * t)
  | If of (t * t * t)
  | While of (t * t)
  | For of (linno * arg * t * string)
  | Case of (linno * arg * case list)
  | Defun of (linno * string * t)
and assign = (string * arg)
and redirection =
    File of (redir_type * int * arg)
  | Dup of (dup_type * int * arg)
  | Heredoc of (heredoc_type * int * arg)
and redir_type = To | Clobber | From | FromTo | Append
and dup_type = ToFD | FromFD
and heredoc_type = Here | XHere
and args = arg list
and arg = arg_char list
and arg_char =
    C of char
  | E of char
  | T of string option
  | A of arg
  | V of (var_type * bool * string * arg)
  | Q of arg
  | B of t
and var_type =
    Normal
  | Minus
  | Plus
  | Question
  | Assign
  | TrimR
  | TrimRMax
  | TrimL
  | TrimLMax
  | Length
and case = { cpattern : arg list; cbody : t; }

val of_node : Dash.node Ctypes.union Ctypes.ptr -> t

(* command that does nothing *)
val skip : t

(* render to string *)
val to_string : t -> string
