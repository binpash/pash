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

and case = Ast_atd_t.case = { cpattern: arg list; cbody: t }

and redirection = Ast_atd_t.redirection

and t = Ast_atd_t.t

let write__7 = (
  Atdgen_runtime.Oj_run.write_option (
    Yojson.Safe.write_string
  )
)
let string_of__7 ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write__7 ob x;
  Bi_outbuf.contents ob
let read__7 = (
  fun p lb ->
    Yojson.Safe.read_space p lb;
    match Yojson.Safe.start_any_variant p lb with
      | `Edgy_bracket -> (
          match Yojson.Safe.read_ident p lb with
            | "None" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              (None : _ option)
            | "Some" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  Atdgen_runtime.Oj_run.read_string
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              (Some x : _ option)
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Double_quote -> (
          match Yojson.Safe.finish_string p lb with
            | "None" ->
              (None : _ option)
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Square_bracket -> (
          match Atdgen_runtime.Oj_run.read_string p lb with
            | "Some" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  Atdgen_runtime.Oj_run.read_string
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              (Some x : _ option)
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
)
let _7_of_string s =
  read__7 (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
let write_char = (
  Yojson.Safe.write_int
)
let string_of_char ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_char ob x;
  Bi_outbuf.contents ob
let read_char = (
  Atdgen_runtime.Oj_run.read_int
)
let char_of_string s =
  read_char (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
let write_dup_type = (
  fun ob x ->
    match x with
      | `ToFD -> Bi_outbuf.add_string ob "<\"ToFD\">"
      | `FromFD -> Bi_outbuf.add_string ob "<\"FromFD\">"
)
let string_of_dup_type ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_dup_type ob x;
  Bi_outbuf.contents ob
let read_dup_type = (
  fun p lb ->
    Yojson.Safe.read_space p lb;
    match Yojson.Safe.start_any_variant p lb with
      | `Edgy_bracket -> (
          match Yojson.Safe.read_ident p lb with
            | "ToFD" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `ToFD
            | "FromFD" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `FromFD
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Double_quote -> (
          match Yojson.Safe.finish_string p lb with
            | "ToFD" ->
              `ToFD
            | "FromFD" ->
              `FromFD
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Square_bracket -> (
          match Atdgen_runtime.Oj_run.read_string p lb with
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
)
let dup_type_of_string s =
  read_dup_type (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
let write_heredoc_type = (
  fun ob x ->
    match x with
      | `Here -> Bi_outbuf.add_string ob "<\"Here\">"
      | `XHere -> Bi_outbuf.add_string ob "<\"XHere\">"
)
let string_of_heredoc_type ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_heredoc_type ob x;
  Bi_outbuf.contents ob
let read_heredoc_type = (
  fun p lb ->
    Yojson.Safe.read_space p lb;
    match Yojson.Safe.start_any_variant p lb with
      | `Edgy_bracket -> (
          match Yojson.Safe.read_ident p lb with
            | "Here" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Here
            | "XHere" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `XHere
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Double_quote -> (
          match Yojson.Safe.finish_string p lb with
            | "Here" ->
              `Here
            | "XHere" ->
              `XHere
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Square_bracket -> (
          match Atdgen_runtime.Oj_run.read_string p lb with
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
)
let heredoc_type_of_string s =
  read_heredoc_type (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
let write_linno = (
  Yojson.Safe.write_int
)
let string_of_linno ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_linno ob x;
  Bi_outbuf.contents ob
let read_linno = (
  Atdgen_runtime.Oj_run.read_int
)
let linno_of_string s =
  read_linno (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
let write_redir_type = (
  fun ob x ->
    match x with
      | `To -> Bi_outbuf.add_string ob "<\"To\">"
      | `Clobber -> Bi_outbuf.add_string ob "<\"Clobber\">"
      | `From -> Bi_outbuf.add_string ob "<\"From\">"
      | `FromTo -> Bi_outbuf.add_string ob "<\"FromTo\">"
      | `Append -> Bi_outbuf.add_string ob "<\"Append\">"
)
let string_of_redir_type ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_redir_type ob x;
  Bi_outbuf.contents ob
let read_redir_type = (
  fun p lb ->
    Yojson.Safe.read_space p lb;
    match Yojson.Safe.start_any_variant p lb with
      | `Edgy_bracket -> (
          match Yojson.Safe.read_ident p lb with
            | "To" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `To
            | "Clobber" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Clobber
            | "From" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `From
            | "FromTo" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `FromTo
            | "Append" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Append
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Double_quote -> (
          match Yojson.Safe.finish_string p lb with
            | "To" ->
              `To
            | "Clobber" ->
              `Clobber
            | "From" ->
              `From
            | "FromTo" ->
              `FromTo
            | "Append" ->
              `Append
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Square_bracket -> (
          match Atdgen_runtime.Oj_run.read_string p lb with
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
)
let redir_type_of_string s =
  read_redir_type (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
let write_var_type = (
  fun ob x ->
    match x with
      | `Normal -> Bi_outbuf.add_string ob "<\"Normal\">"
      | `Minus -> Bi_outbuf.add_string ob "<\"Minus\">"
      | `Plus -> Bi_outbuf.add_string ob "<\"Plus\">"
      | `Question -> Bi_outbuf.add_string ob "<\"Question\">"
      | `Assign -> Bi_outbuf.add_string ob "<\"Assign\">"
      | `TrimR -> Bi_outbuf.add_string ob "<\"TrimR\">"
      | `TrimRMax -> Bi_outbuf.add_string ob "<\"TrimRMax\">"
      | `TrimL -> Bi_outbuf.add_string ob "<\"TrimL\">"
      | `TrimLMax -> Bi_outbuf.add_string ob "<\"TrimLMax\">"
      | `Length -> Bi_outbuf.add_string ob "<\"Length\">"
)
let string_of_var_type ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_var_type ob x;
  Bi_outbuf.contents ob
let read_var_type = (
  fun p lb ->
    Yojson.Safe.read_space p lb;
    match Yojson.Safe.start_any_variant p lb with
      | `Edgy_bracket -> (
          match Yojson.Safe.read_ident p lb with
            | "Normal" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Normal
            | "Minus" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Minus
            | "Plus" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Plus
            | "Question" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Question
            | "Assign" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Assign
            | "TrimR" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `TrimR
            | "TrimRMax" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `TrimRMax
            | "TrimL" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `TrimL
            | "TrimLMax" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `TrimLMax
            | "Length" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Length
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Double_quote -> (
          match Yojson.Safe.finish_string p lb with
            | "Normal" ->
              `Normal
            | "Minus" ->
              `Minus
            | "Plus" ->
              `Plus
            | "Question" ->
              `Question
            | "Assign" ->
              `Assign
            | "TrimR" ->
              `TrimR
            | "TrimRMax" ->
              `TrimRMax
            | "TrimL" ->
              `TrimL
            | "TrimLMax" ->
              `TrimLMax
            | "Length" ->
              `Length
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Square_bracket -> (
          match Atdgen_runtime.Oj_run.read_string p lb with
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
)
let var_type_of_string s =
  read_var_type (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
let rec write__1 ob x = (
  Atdgen_runtime.Oj_run.write_list (
    write_assign
  )
) ob x
and string_of__1 ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write__1 ob x;
  Bi_outbuf.contents ob
and write__2 ob x = (
  Atdgen_runtime.Oj_run.write_list (
    write_redirection
  )
) ob x
and string_of__2 ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write__2 ob x;
  Bi_outbuf.contents ob
and write__3 ob x = (
  Atdgen_runtime.Oj_run.write_list (
    write_t
  )
) ob x
and string_of__3 ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write__3 ob x;
  Bi_outbuf.contents ob
and write__4 ob x = (
  Atdgen_runtime.Oj_run.write_list (
    write_case
  )
) ob x
and string_of__4 ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write__4 ob x;
  Bi_outbuf.contents ob
and write__5 ob x = (
  Atdgen_runtime.Oj_run.write_list (
    write_arg
  )
) ob x
and string_of__5 ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write__5 ob x;
  Bi_outbuf.contents ob
and write__6 ob x = (
  Atdgen_runtime.Oj_run.write_list (
    write_arg_char
  )
) ob x
and string_of__6 ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write__6 ob x;
  Bi_outbuf.contents ob
and write_arg ob x = (
  write__6
) ob x
and string_of_arg ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_arg ob x;
  Bi_outbuf.contents ob
and write_arg_char = (
  fun ob x ->
    match x with
      | `C x ->
        Bi_outbuf.add_string ob "<\"C\":";
        (
          write_char
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `E x ->
        Bi_outbuf.add_string ob "<\"E\":";
        (
          write_char
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `T x ->
        Bi_outbuf.add_string ob "<\"T\":";
        (
          write__7
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `A x ->
        Bi_outbuf.add_string ob "<\"A\":";
        (
          write_arg
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `V x ->
        Bi_outbuf.add_string ob "<\"V\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _, _ = x in
            (
              write_var_type
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _, _ = x in
            (
              Yojson.Safe.write_bool
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x, _ = x in
            (
              Yojson.Safe.write_string
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, _, x = x in
            (
              write_arg
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Q x ->
        Bi_outbuf.add_string ob "<\"Q\":";
        (
          write_arg
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `B x ->
        Bi_outbuf.add_string ob "<\"B\":";
        (
          write_t
        ) ob x;
        Bi_outbuf.add_char ob '>'
)
and string_of_arg_char ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_arg_char ob x;
  Bi_outbuf.contents ob
and write_args ob x = (
  write__5
) ob x
and string_of_args ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_args ob x;
  Bi_outbuf.contents ob
and write_assign = (
  fun ob x ->
    Bi_outbuf.add_char ob '(';
    (let x, _ = x in
    (
      Yojson.Safe.write_string
    ) ob x
    );
    Bi_outbuf.add_char ob ',';
    (let _, x = x in
    (
      write_arg
    ) ob x
    );
    Bi_outbuf.add_char ob ')';
)
and string_of_assign ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_assign ob x;
  Bi_outbuf.contents ob
and write_case : _ -> case -> _ = (
  fun ob x ->
    Bi_outbuf.add_char ob '{';
    let is_first = ref true in
    if !is_first then
      is_first := false
    else
      Bi_outbuf.add_char ob ',';
    Bi_outbuf.add_string ob "\"cpattern\":";
    (
      write__5
    )
      ob x.cpattern;
    if !is_first then
      is_first := false
    else
      Bi_outbuf.add_char ob ',';
    Bi_outbuf.add_string ob "\"cbody\":";
    (
      write_t
    )
      ob x.cbody;
    Bi_outbuf.add_char ob '}';
)
and string_of_case ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_case ob x;
  Bi_outbuf.contents ob
and write_redirection = (
  fun ob x ->
    match x with
      | `File x ->
        Bi_outbuf.add_string ob "<\"File\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _ = x in
            (
              write_redir_type
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _ = x in
            (
              Yojson.Safe.write_int
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x = x in
            (
              write_arg
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Dup x ->
        Bi_outbuf.add_string ob "<\"Dup\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _ = x in
            (
              write_dup_type
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _ = x in
            (
              Yojson.Safe.write_int
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x = x in
            (
              write_arg
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Heredoc x ->
        Bi_outbuf.add_string ob "<\"Heredoc\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _ = x in
            (
              write_heredoc_type
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _ = x in
            (
              Yojson.Safe.write_int
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x = x in
            (
              write_arg
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
)
and string_of_redirection ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_redirection ob x;
  Bi_outbuf.contents ob
and write_t = (
  fun ob x ->
    match x with
      | `Command x ->
        Bi_outbuf.add_string ob "<\"Command\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _, _ = x in
            (
              write_linno
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _, _ = x in
            (
              write__1
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x, _ = x in
            (
              write_args
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, _, x = x in
            (
              write__2
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Pipe x ->
        Bi_outbuf.add_string ob "<\"Pipe\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _ = x in
            (
              Yojson.Safe.write_bool
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x = x in
            (
              write__3
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Redir x ->
        Bi_outbuf.add_string ob "<\"Redir\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _ = x in
            (
              write_linno
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _ = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x = x in
            (
              write__2
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Background x ->
        Bi_outbuf.add_string ob "<\"Background\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _ = x in
            (
              write_linno
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _ = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x = x in
            (
              write__2
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Subshell x ->
        Bi_outbuf.add_string ob "<\"Subshell\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _ = x in
            (
              write_linno
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _ = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x = x in
            (
              write__2
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `And x ->
        Bi_outbuf.add_string ob "<\"And\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _ = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Or x ->
        Bi_outbuf.add_string ob "<\"Or\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _ = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Not x ->
        Bi_outbuf.add_string ob "<\"Not\":";
        (
          write_t
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Semi x ->
        Bi_outbuf.add_string ob "<\"Semi\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _ = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `If x ->
        Bi_outbuf.add_string ob "<\"If\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _ = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _ = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `While x ->
        Bi_outbuf.add_string ob "<\"While\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _ = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `For x ->
        Bi_outbuf.add_string ob "<\"For\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _, _ = x in
            (
              write_linno
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _, _ = x in
            (
              write_arg
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x, _ = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, _, x = x in
            (
              Yojson.Safe.write_string
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Case x ->
        Bi_outbuf.add_string ob "<\"Case\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _ = x in
            (
              write_linno
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _ = x in
            (
              write_arg
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x = x in
            (
              write__4
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
      | `Defun x ->
        Bi_outbuf.add_string ob "<\"Defun\":";
        (
          fun ob x ->
            Bi_outbuf.add_char ob '(';
            (let x, _, _ = x in
            (
              write_linno
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, x, _ = x in
            (
              Yojson.Safe.write_string
            ) ob x
            );
            Bi_outbuf.add_char ob ',';
            (let _, _, x = x in
            (
              write_t
            ) ob x
            );
            Bi_outbuf.add_char ob ')';
        ) ob x;
        Bi_outbuf.add_char ob '>'
)
and string_of_t ?(len = 1024) x =
  let ob = Bi_outbuf.create len in
  write_t ob x;
  Bi_outbuf.contents ob
let rec read__1 p lb = (
  Atdgen_runtime.Oj_run.read_list (
    read_assign
  )
) p lb
and _1_of_string s =
  read__1 (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read__2 p lb = (
  Atdgen_runtime.Oj_run.read_list (
    read_redirection
  )
) p lb
and _2_of_string s =
  read__2 (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read__3 p lb = (
  Atdgen_runtime.Oj_run.read_list (
    read_t
  )
) p lb
and _3_of_string s =
  read__3 (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read__4 p lb = (
  Atdgen_runtime.Oj_run.read_list (
    read_case
  )
) p lb
and _4_of_string s =
  read__4 (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read__5 p lb = (
  Atdgen_runtime.Oj_run.read_list (
    read_arg
  )
) p lb
and _5_of_string s =
  read__5 (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read__6 p lb = (
  Atdgen_runtime.Oj_run.read_list (
    read_arg_char
  )
) p lb
and _6_of_string s =
  read__6 (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read_arg p lb = (
  read__6
) p lb
and arg_of_string s =
  read_arg (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read_arg_char = (
  fun p lb ->
    Yojson.Safe.read_space p lb;
    match Yojson.Safe.start_any_variant p lb with
      | `Edgy_bracket -> (
          match Yojson.Safe.read_ident p lb with
            | "C" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  read_char
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `C x
            | "E" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  read_char
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `E x
            | "T" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  read__7
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `T x
            | "A" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  read_arg
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `A x
            | "V" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_var_type
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_bool
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_string
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x3 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2, x3)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2; 3 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `V x
            | "Q" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  read_arg
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Q x
            | "B" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  read_t
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `B x
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Double_quote -> (
          match Yojson.Safe.finish_string p lb with
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Square_bracket -> (
          match Atdgen_runtime.Oj_run.read_string p lb with
            | "C" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  read_char
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `C x
            | "E" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  read_char
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `E x
            | "T" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  read__7
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `T x
            | "A" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  read_arg
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `A x
            | "V" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_var_type
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_bool
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_string
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x3 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2, x3)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2; 3 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `V x
            | "Q" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  read_arg
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Q x
            | "B" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  read_t
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `B x
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
)
and arg_char_of_string s =
  read_arg_char (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read_args p lb = (
  read__5
) p lb
and args_of_string s =
  read_args (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read_assign = (
  fun p lb ->
    Yojson.Safe.read_space p lb;
    let std_tuple = Yojson.Safe.start_any_tuple p lb in
    let len = ref 0 in
    let end_of_tuple = ref false in
    (try
      let x0 =
        let x =
          (
            Atdgen_runtime.Oj_run.read_string
          ) p lb
        in
        incr len;
        Yojson.Safe.read_space p lb;
        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
        x
      in
      let x1 =
        let x =
          (
            read_arg
          ) p lb
        in
        incr len;
        (try
          Yojson.Safe.read_space p lb;
          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
        with Yojson.End_of_tuple -> end_of_tuple := true);
        x
      in
      if not !end_of_tuple then (
        try
          while true do
            Yojson.Safe.skip_json p lb;
            Yojson.Safe.read_space p lb;
            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
          done
        with Yojson.End_of_tuple -> ()
      );
      (x0, x1)
    with Yojson.End_of_tuple ->
      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1 ]);
)
and assign_of_string s =
  read_assign (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read_case = (
  fun p lb ->
    Yojson.Safe.read_space p lb;
    Yojson.Safe.read_lcurl p lb;
    let field_cpattern = ref (Obj.magic (Sys.opaque_identity 0.0)) in
    let field_cbody = ref (Obj.magic (Sys.opaque_identity 0.0)) in
    let bits0 = ref 0 in
    try
      Yojson.Safe.read_space p lb;
      Yojson.Safe.read_object_end lb;
      Yojson.Safe.read_space p lb;
      let f =
        fun s pos len ->
          if pos < 0 || len < 0 || pos + len > String.length s then
            invalid_arg "out-of-bounds substring position or length";
          match len with
            | 5 -> (
                if String.unsafe_get s pos = 'c' && String.unsafe_get s (pos+1) = 'b' && String.unsafe_get s (pos+2) = 'o' && String.unsafe_get s (pos+3) = 'd' && String.unsafe_get s (pos+4) = 'y' then (
                  1
                )
                else (
                  -1
                )
              )
            | 8 -> (
                if String.unsafe_get s pos = 'c' && String.unsafe_get s (pos+1) = 'p' && String.unsafe_get s (pos+2) = 'a' && String.unsafe_get s (pos+3) = 't' && String.unsafe_get s (pos+4) = 't' && String.unsafe_get s (pos+5) = 'e' && String.unsafe_get s (pos+6) = 'r' && String.unsafe_get s (pos+7) = 'n' then (
                  0
                )
                else (
                  -1
                )
              )
            | _ -> (
                -1
              )
      in
      let i = Yojson.Safe.map_ident p f lb in
      Atdgen_runtime.Oj_run.read_until_field_value p lb;
      (
        match i with
          | 0 ->
            field_cpattern := (
              (
                read__5
              ) p lb
            );
            bits0 := !bits0 lor 0x1;
          | 1 ->
            field_cbody := (
              (
                read_t
              ) p lb
            );
            bits0 := !bits0 lor 0x2;
          | _ -> (
              Yojson.Safe.skip_json p lb
            )
      );
      while true do
        Yojson.Safe.read_space p lb;
        Yojson.Safe.read_object_sep p lb;
        Yojson.Safe.read_space p lb;
        let f =
          fun s pos len ->
            if pos < 0 || len < 0 || pos + len > String.length s then
              invalid_arg "out-of-bounds substring position or length";
            match len with
              | 5 -> (
                  if String.unsafe_get s pos = 'c' && String.unsafe_get s (pos+1) = 'b' && String.unsafe_get s (pos+2) = 'o' && String.unsafe_get s (pos+3) = 'd' && String.unsafe_get s (pos+4) = 'y' then (
                    1
                  )
                  else (
                    -1
                  )
                )
              | 8 -> (
                  if String.unsafe_get s pos = 'c' && String.unsafe_get s (pos+1) = 'p' && String.unsafe_get s (pos+2) = 'a' && String.unsafe_get s (pos+3) = 't' && String.unsafe_get s (pos+4) = 't' && String.unsafe_get s (pos+5) = 'e' && String.unsafe_get s (pos+6) = 'r' && String.unsafe_get s (pos+7) = 'n' then (
                    0
                  )
                  else (
                    -1
                  )
                )
              | _ -> (
                  -1
                )
        in
        let i = Yojson.Safe.map_ident p f lb in
        Atdgen_runtime.Oj_run.read_until_field_value p lb;
        (
          match i with
            | 0 ->
              field_cpattern := (
                (
                  read__5
                ) p lb
              );
              bits0 := !bits0 lor 0x1;
            | 1 ->
              field_cbody := (
                (
                  read_t
                ) p lb
              );
              bits0 := !bits0 lor 0x2;
            | _ -> (
                Yojson.Safe.skip_json p lb
              )
        );
      done;
      assert false;
    with Yojson.End_of_object -> (
        if !bits0 <> 0x3 then Atdgen_runtime.Oj_run.missing_fields p [| !bits0 |] [| "cpattern"; "cbody" |];
        (
          {
            cpattern = !field_cpattern;
            cbody = !field_cbody;
          }
         : case)
      )
)
and case_of_string s =
  read_case (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read_redirection = (
  fun p lb ->
    Yojson.Safe.read_space p lb;
    match Yojson.Safe.start_any_variant p lb with
      | `Edgy_bracket -> (
          match Yojson.Safe.read_ident p lb with
            | "File" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_redir_type
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_int
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `File x
            | "Dup" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_dup_type
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_int
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Dup x
            | "Heredoc" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_heredoc_type
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_int
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Heredoc x
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Double_quote -> (
          match Yojson.Safe.finish_string p lb with
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Square_bracket -> (
          match Atdgen_runtime.Oj_run.read_string p lb with
            | "File" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_redir_type
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_int
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `File x
            | "Dup" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_dup_type
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_int
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Dup x
            | "Heredoc" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_heredoc_type
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_int
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Heredoc x
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
)
and redirection_of_string s =
  read_redirection (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
and read_t = (
  fun p lb ->
    Yojson.Safe.read_space p lb;
    match Yojson.Safe.start_any_variant p lb with
      | `Edgy_bracket -> (
          match Yojson.Safe.read_ident p lb with
            | "Command" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read__1
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_args
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x3 =
                        let x =
                          (
                            read__2
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2, x3)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2; 3 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Command x
            | "Pipe" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_bool
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read__3
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Pipe x
            | "Redir" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read__2
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Redir x
            | "Background" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read__2
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Background x
            | "Subshell" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read__2
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Subshell x
            | "And" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `And x
            | "Or" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Or x
            | "Not" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  read_t
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Not x
            | "Semi" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Semi x
            | "If" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `If x
            | "While" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `While x
            | "For" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x3 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_string
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2, x3)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2; 3 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `For x
            | "Case" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read__4
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Case x
            | "Defun" ->
              Atdgen_runtime.Oj_run.read_until_field_value p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_string
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_gt p lb;
              `Defun x
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Double_quote -> (
          match Yojson.Safe.finish_string p lb with
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
      | `Square_bracket -> (
          match Atdgen_runtime.Oj_run.read_string p lb with
            | "Command" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read__1
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_args
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x3 =
                        let x =
                          (
                            read__2
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2, x3)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2; 3 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Command x
            | "Pipe" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_bool
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read__3
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Pipe x
            | "Redir" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read__2
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Redir x
            | "Background" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read__2
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Background x
            | "Subshell" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read__2
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Subshell x
            | "And" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `And x
            | "Or" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Or x
            | "Not" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  read_t
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Not x
            | "Semi" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Semi x
            | "If" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `If x
            | "While" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `While x
            | "For" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x3 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_string
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2, x3)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2; 3 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `For x
            | "Case" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            read_arg
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read__4
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Case x
            | "Defun" ->
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_comma p lb;
              Yojson.Safe.read_space p lb;
              let x = (
                  fun p lb ->
                    Yojson.Safe.read_space p lb;
                    let std_tuple = Yojson.Safe.start_any_tuple p lb in
                    let len = ref 0 in
                    let end_of_tuple = ref false in
                    (try
                      let x0 =
                        let x =
                          (
                            read_linno
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x1 =
                        let x =
                          (
                            Atdgen_runtime.Oj_run.read_string
                          ) p lb
                        in
                        incr len;
                        Yojson.Safe.read_space p lb;
                        Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        x
                      in
                      let x2 =
                        let x =
                          (
                            read_t
                          ) p lb
                        in
                        incr len;
                        (try
                          Yojson.Safe.read_space p lb;
                          Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                        with Yojson.End_of_tuple -> end_of_tuple := true);
                        x
                      in
                      if not !end_of_tuple then (
                        try
                          while true do
                            Yojson.Safe.skip_json p lb;
                            Yojson.Safe.read_space p lb;
                            Yojson.Safe.read_tuple_sep2 p std_tuple lb;
                          done
                        with Yojson.End_of_tuple -> ()
                      );
                      (x0, x1, x2)
                    with Yojson.End_of_tuple ->
                      Atdgen_runtime.Oj_run.missing_tuple_fields p !len [ 0; 1; 2 ]);
                ) p lb
              in
              Yojson.Safe.read_space p lb;
              Yojson.Safe.read_rbr p lb;
              `Defun x
            | x ->
              Atdgen_runtime.Oj_run.invalid_variant_tag p x
        )
)
and t_of_string s =
  read_t (Yojson.Safe.init_lexer ()) (Lexing.from_string s)
