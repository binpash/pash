(* Auto-generated from "ast_atd.atd" *)
[@@@ocaml.warning "-27-32-35-39"]

type char = Ast_atd_t.char

type dup_type = Ast_atd_t.dup_type

type heredoc_type = Ast_atd_t.heredoc_type

type linno = Ast_atd_t.linno

type redir_type = Ast_atd_t.redir_type

type var_type = Ast_atd_t.var_type

type arg = Ast_atd_t.arg

and arg_char = Ast_atd_t.arg_char

and args = Ast_atd_t.args

and assign = Ast_atd_t.assign

and redirection = Ast_atd_t.redirection

and t = Ast_atd_t.t

val write_char :
  Bi_outbuf.t -> char -> unit
  (** Output a JSON value of type {!char}. *)

val string_of_char :
  ?len:int -> char -> string
  (** Serialize a value of type {!char}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_char :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> char
  (** Input JSON data of type {!char}. *)

val char_of_string :
  string -> char
  (** Deserialize JSON data of type {!char}. *)

val write_dup_type :
  Bi_outbuf.t -> dup_type -> unit
  (** Output a JSON value of type {!dup_type}. *)

val string_of_dup_type :
  ?len:int -> dup_type -> string
  (** Serialize a value of type {!dup_type}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_dup_type :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> dup_type
  (** Input JSON data of type {!dup_type}. *)

val dup_type_of_string :
  string -> dup_type
  (** Deserialize JSON data of type {!dup_type}. *)

val write_heredoc_type :
  Bi_outbuf.t -> heredoc_type -> unit
  (** Output a JSON value of type {!heredoc_type}. *)

val string_of_heredoc_type :
  ?len:int -> heredoc_type -> string
  (** Serialize a value of type {!heredoc_type}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_heredoc_type :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> heredoc_type
  (** Input JSON data of type {!heredoc_type}. *)

val heredoc_type_of_string :
  string -> heredoc_type
  (** Deserialize JSON data of type {!heredoc_type}. *)

val write_linno :
  Bi_outbuf.t -> linno -> unit
  (** Output a JSON value of type {!linno}. *)

val string_of_linno :
  ?len:int -> linno -> string
  (** Serialize a value of type {!linno}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_linno :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> linno
  (** Input JSON data of type {!linno}. *)

val linno_of_string :
  string -> linno
  (** Deserialize JSON data of type {!linno}. *)

val write_redir_type :
  Bi_outbuf.t -> redir_type -> unit
  (** Output a JSON value of type {!redir_type}. *)

val string_of_redir_type :
  ?len:int -> redir_type -> string
  (** Serialize a value of type {!redir_type}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_redir_type :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> redir_type
  (** Input JSON data of type {!redir_type}. *)

val redir_type_of_string :
  string -> redir_type
  (** Deserialize JSON data of type {!redir_type}. *)

val write_var_type :
  Bi_outbuf.t -> var_type -> unit
  (** Output a JSON value of type {!var_type}. *)

val string_of_var_type :
  ?len:int -> var_type -> string
  (** Serialize a value of type {!var_type}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_var_type :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> var_type
  (** Input JSON data of type {!var_type}. *)

val var_type_of_string :
  string -> var_type
  (** Deserialize JSON data of type {!var_type}. *)

val write_arg :
  Bi_outbuf.t -> arg -> unit
  (** Output a JSON value of type {!arg}. *)

val string_of_arg :
  ?len:int -> arg -> string
  (** Serialize a value of type {!arg}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_arg :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> arg
  (** Input JSON data of type {!arg}. *)

val arg_of_string :
  string -> arg
  (** Deserialize JSON data of type {!arg}. *)

val write_arg_char :
  Bi_outbuf.t -> arg_char -> unit
  (** Output a JSON value of type {!arg_char}. *)

val string_of_arg_char :
  ?len:int -> arg_char -> string
  (** Serialize a value of type {!arg_char}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_arg_char :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> arg_char
  (** Input JSON data of type {!arg_char}. *)

val arg_char_of_string :
  string -> arg_char
  (** Deserialize JSON data of type {!arg_char}. *)

val write_args :
  Bi_outbuf.t -> args -> unit
  (** Output a JSON value of type {!args}. *)

val string_of_args :
  ?len:int -> args -> string
  (** Serialize a value of type {!args}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_args :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> args
  (** Input JSON data of type {!args}. *)

val args_of_string :
  string -> args
  (** Deserialize JSON data of type {!args}. *)

val write_assign :
  Bi_outbuf.t -> assign -> unit
  (** Output a JSON value of type {!assign}. *)

val string_of_assign :
  ?len:int -> assign -> string
  (** Serialize a value of type {!assign}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_assign :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> assign
  (** Input JSON data of type {!assign}. *)

val assign_of_string :
  string -> assign
  (** Deserialize JSON data of type {!assign}. *)

val write_redirection :
  Bi_outbuf.t -> redirection -> unit
  (** Output a JSON value of type {!redirection}. *)

val string_of_redirection :
  ?len:int -> redirection -> string
  (** Serialize a value of type {!redirection}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_redirection :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> redirection
  (** Input JSON data of type {!redirection}. *)

val redirection_of_string :
  string -> redirection
  (** Deserialize JSON data of type {!redirection}. *)

val write_t :
  Bi_outbuf.t -> t -> unit
  (** Output a JSON value of type {!t}. *)

val string_of_t :
  ?len:int -> t -> string
  (** Serialize a value of type {!t}
      into a JSON string.
      @param len specifies the initial length
                 of the buffer used internally.
                 Default: 1024. *)

val read_t :
  Yojson.Safe.lexer_state -> Lexing.lexbuf -> t
  (** Input JSON data of type {!t}. *)

val t_of_string :
  string -> t
  (** Deserialize JSON data of type {!t}. *)

