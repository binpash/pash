#!/usr/bin/python3


# let separated f l = intercalate " " (List.map f l)
# 
# let show_unless expected actual =
#   if expected = actual
#   then ""
#   else string_of_int actual
# 
# let background s = "{ " ^ s ^ " & }"
#                             


def separated (f, l):
    return " ".join (map (f, l));



def to_string (ast):
    print (ast);

    if (len (ast) == 0):
        pass;
    else:
        assert (len (ast) == 2);

        type = ast [0];
        params = ast [1]

        if (type == "Command"):
            print ("Command");

            (_, assigns, cmds, redirs) = params;
            str = separated (string_of_assign, assigns);
            if ((len (assigns) == 0) or (len (cmds) == 0)):
                pass;
            else:
                str = str + " ";
            str = str + separated (string_of_arg, cmds) + string_of_redirs (redirs);

            return (str);
        elif (type == "Pipe"):
            print ("Type")
        elif (type == "If"):
            (c, t, e) = params;
            print ("If")

            return string_of_if (c, t, e);

    return "";


def string_of_if (c, t, e):
    str = "if " + to_string (c) + \
          "; then " + to_string (t);

    # TODO: uncommon cases
    str = str + "; else " + to_string (e) + "; fi";

    return (str);



# let rec to_string = function
#   | Command (_,assigns,cmds,redirs) ->
#      separated string_of_assign assigns ^
#      (if List.length assigns = 0 || List.length cmds = 0 then "" else " ") ^
#      separated string_of_arg cmds ^ string_of_redirs redirs
#   | Pipe (bg,ps) ->
#      let p = intercalate " | " (List.map to_string ps) in
#      if bg then background p else p
#   | Redir (_,a,redirs) ->
#      to_string a ^ string_of_redirs redirs
#   | Background (_,a,redirs) ->
#      (* we translate 
#            cmds... &
#         to
#            { cmds & }
#         this avoids issues with parsing; in particular,
#           cmd1 & ; cmd2 & ; cmd3
#         doesn't parse; it must be:
#           cmd1 & cmd2 & cmd3
#         it's a little too annoying to track "was the last thing
#         backgrounded?" so the braces resolve the issue. testing
#         indicates that they're semantically equivalent.
#       *)
#      background (to_string a ^ string_of_redirs redirs)
#   | Subshell (_,a,redirs) ->
#      parens (to_string a ^ string_of_redirs redirs)
#   | And (a1,a2) -> to_string  a1 ^ " && " ^ to_string a2
#   | Or (a1,a2) -> to_string a1 ^ " || " ^ to_string a2
#   | Not a -> "! " ^ braces (to_string a)
#   | Semi (a1,a2) -> to_string a1 ^ " ; " ^ to_string a2
#   | If (c,t,e) -> string_of_if c t e
#   | While (Not t,b) ->
#      "until " ^ to_string t ^ "; do " ^ to_string b ^ "; done "
#   | While (t,b) ->
#      "while " ^ to_string t ^ "; do " ^ to_string b ^ "; done "
#   | For (_,a,body,var) ->
#      "for " ^ var ^ " in " ^ string_of_arg a ^ "; do " ^
#      to_string body ^ "; done"
#   | Case (_,a,cs) ->
#      "case " ^ string_of_arg a ^ " in " ^
#      separated string_of_case cs ^ " esac"
#   | Defun (_,name,body) -> name ^ "() {\n" ^ to_string body ^ "\n}"





# and string_of_if c t e =
#   "if " ^ to_string c ^
#   "; then " ^ to_string t ^
#   (match e with
#    | Command (-1,[],[],[]) -> "; fi" (* one-armed if *)
#    | If (c,t,e) -> "; el" ^ string_of_if c t e
#    | _ -> "; else " ^ to_string e ^ "; fi")
#                                                  
# and string_of_arg_char = function
#   | E '\'' -> "\\'"
#   | E '\"' -> "\\\""
#   | E '(' -> "\\("
#   | E ')' -> "\\)"
#   | E '{' -> "\\{"
#   | E '}' -> "\\}"
#   | E '$' -> "\\$"
#   | E '!' -> "\\!"
#   | E '&' -> "\\&"
#   | E '|' -> "\\|" 
#   | E ';' -> "\\;" 
#   | C c -> String.make 1 c
#   | E c -> Char.escaped c
#   | T None -> "~"
#   | T (Some u) -> "~" ^ u
#   | A a -> "$((" ^ string_of_arg a ^ "))"
#   | V (Length,_,name,_) -> "${#" ^ name ^ "}"
#   | V (vt,nul,name,a) ->
#      "${" ^ name ^ (if nul then ":" else "") ^ string_of_var_type vt ^ string_of_arg a ^ "}"
#   | Q a -> "\"" ^ string_of_arg a ^ "\""
#   | B t -> "$(" ^ to_string t ^ ")"


# and string_of_arg = function
#   | [] -> ""
#   | c :: a -> string_of_arg_char c ^ string_of_arg a
def string_of_arg (arg):
    return "TODO";


# and string_of_assign (v,a) = v ^ "=" ^ string_of_arg a
def string_of_assign (both):
    (v, a) = both;
    return v + "=" + string_of_arg (a);


#                                                    
# and string_of_case c =
#   let pats = List.map string_of_arg c.cpattern in
#   intercalate "|" pats ^ ") " ^ to_string c.cbody ^ ";;"
# 
# and string_of_redir = function
#   | File (To,fd,a)      -> show_unless 1 fd ^ ">" ^ string_of_arg a
#   | File (Clobber,fd,a) -> show_unless 1 fd ^ ">|" ^ string_of_arg a
#   | File (From,fd,a)    -> show_unless 0 fd ^ "<" ^ string_of_arg a
#   | File (FromTo,fd,a)  -> show_unless 0 fd ^ "<>" ^ string_of_arg a
#   | File (Append,fd,a)  -> show_unless 1 fd ^ ">>" ^ string_of_arg a
#   | Dup (ToFD,fd,tgt)   -> show_unless 1 fd ^ ">&" ^ string_of_arg tgt
#   | Dup (FromFD,fd,tgt) -> show_unless 0 fd ^ "<&" ^ string_of_arg tgt
#   | Heredoc (t,fd,a) ->
#      let heredoc = string_of_arg a in
#      let marker = fresh_marker (lines heredoc) "EOF" in
#      show_unless 0 fd ^ "<<" ^
#      (if t = XHere then marker else "'" ^ marker ^ "'") ^ "\n" ^ heredoc ^ marker ^ "\n"
#                                                                                
# and string_of_redirs rs =
#   let ss = List.map string_of_redir rs in
#   (if List.length ss > 0 then " " else "") ^ intercalate " " ss
def string_of_redirs (rs):
    return "TODO"

