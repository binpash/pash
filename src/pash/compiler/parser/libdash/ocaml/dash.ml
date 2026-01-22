open Printf
open Ctypes
open Ctypes_types
open Foreign

(* First, some dash trivia. *)
   
type stackmark
              
let stackmark : stackmark structure typ = structure "stackmark"
let stackp = field stackmark "stackp" (ptr void)
let nxt = field stackmark "nxt" string
let size = field stackmark "stacknleft" PosixTypes.size_t
let () = seal stackmark

let init_stack () =
  let stack = make stackmark in
  foreign "setstackmark" (ptr stackmark @-> returning void) (addr stack);
  stack

let pop_stack stack =
  foreign "popstackmark" (ptr stackmark @-> returning void) (addr stack)

let alloc_stack_string =
  foreign "sstrdup" (string @-> returning (ptr char))

let free_stack_string s =
  foreign "stunalloc" (ptr char @-> returning void) s
  
let dash_init : unit -> unit = foreign "init" (void @-> returning void)
let initialize_dash_errno : unit -> unit = 
  foreign "initialize_dash_errno" (void @-> returning void)

let root_stackmark = ref None
let initialize () =
  initialize_dash_errno ();
  dash_init ();
  root_stackmark := Some (init_stack ())

let popfile : unit -> unit =
  foreign "popfile" (void @-> returning void)
                                  
let setinputstring : char ptr -> unit =
  foreign "setinputstring" (ptr char @-> returning void)

let setinputtostdin () : unit =
  foreign "setinputfd" (int @-> int @-> returning void) 0 0 (* don't both pushing the file *)

let setinputfile ?push:(push=false) (s : string) : unit =
  let _ = foreign "setinputfile" (string @-> int @-> returning int) s (if push then 1 else 0) in
  ()

let setvar (x : string) (v : string) : unit =
  let _ = foreign "setvar" (string @-> string @-> int @-> returning (ptr void)) x v 0 in
  ()
          
let setalias (name : string) (mapping : string) : unit =
  foreign "setalias" (string @-> string @-> returning void) name mapping

let unalias (name : string) : unit =
  foreign "unalias" (string @-> returning void) name

(* Next, a utility function that isn't in Unix or ExtUnix. *)

let freshfd_ge10 (fd : int) : int =
  foreign "freshfd_ge10" (int @-> returning int) fd
  
(* Actual AST stuff begins here. *)
(* first, we define the node type... *)
          
type node       
let node : node union typ = union "node"
let node_type = field node "type" int
(* but we don't seal it yet! *)

type nodelist
let nodelist : nodelist structure typ = structure "nodelist"       
let nodelist_next = field nodelist "next" (ptr nodelist)
let nodelist_n = field nodelist "n" (ptr node)
let () = seal nodelist
                       
type ncmd

let ncmd : ncmd structure typ = structure "ncmd"
let ncmd_type = field ncmd "type" int
let ncmd_linno = field ncmd "linno" int
let ncmd_assign = field ncmd "assign" (ptr node)
let ncmd_args = field ncmd "args" (ptr node)
let ncmd_redirect = field ncmd "redirect" (ptr node)
let () = seal ncmd

let node_ncmd = field node "ncmd" ncmd

type npipe

let npipe : npipe structure typ = structure "npipe"
let npipe_type = field npipe "type" int
let npipe_backgnd = field npipe "backgnd" int
let npipe_cmdlist = field npipe "cmdlist" (ptr nodelist)
let () = seal npipe

let node_npipe = field node "npipe" npipe
                           
type nredir

let nredir : nredir structure typ = structure "nredir"
let nredir_type = field nredir "type" int
let nredir_linno = field nredir "linno" int
let nredir_n = field nredir "n" (ptr node)
let nredir_redirect = field nredir "redirect" (ptr node)
let () = seal nredir

let node_nredir = field node "nredir" nredir

type nbinary

let nbinary : nbinary structure typ = structure "nbinary"
let nbinary_type = field nbinary "type" int
let nbinary_ch1 = field nbinary "ch1" (ptr node)
let nbinary_ch2 = field nbinary "ch2" (ptr node)
let () = seal nbinary

let node_nbinary = field node "nbinary" nbinary

type nif

let nif : nif structure typ = structure "nif"
let nif_type = field nif "type" int
let nif_test = field nif "test" (ptr node)
let nif_ifpart = field nif "ifpart" (ptr node)
let nif_elsepart = field nif "elsepart" (ptr node)
let () = seal nif

let node_nif = field node "nif" nif

type nfor

let nfor : nfor structure typ = structure "nfor"
let nfor_type = field nfor "type" int
let nfor_linno = field nfor "linno" int
let nfor_args = field nfor "args" (ptr node)
let nfor_body = field nfor "body" (ptr node)
let nfor_var = field nfor "var" string
let () = seal nfor

let node_nfor = field node "nfor" nfor

type ncase

let ncase : ncase structure typ = structure "ncase"
let ncase_type = field ncase "type" int
let ncase_linno = field ncase "linno" int
let ncase_expr = field ncase "expr" (ptr node)
let ncase_cases = field ncase "cases" (ptr node)
let () = seal ncase

let node_ncase = field node "ncase" ncase

type nclist

let nclist : nclist structure typ = structure "nclist"
let nclist_type = field nclist "type" int
let nclist_next = field nclist "next" (ptr node)
let nclist_pattern = field nclist "pattern" (ptr node)
let nclist_body = field nclist "body" (ptr node)
let () = seal nclist

let node_nclist = field node "nclist" nclist

type ndefun

let ndefun : ndefun structure typ = structure "ndefun"
let ndefun_type = field ndefun "type" int
let ndefun_linno = field ndefun "linno" int
let ndefun_text = field ndefun "text" string
let ndefun_body = field ndefun "body" (ptr node)
let () = seal ndefun

let node_ndefun = field node "ndefun" ndefun

type narg

let narg : narg structure typ = structure "narg"
let narg_type = field narg "type" int
let narg_next = field narg "next" (ptr node)
let narg_text = field narg "text" string
let narg_backquote = field narg "backquote" (ptr nodelist)
let () = seal narg

let node_narg = field node "narg" narg

type nfile

let nfile : nfile structure typ = structure "nfile"
let nfile_type = field nfile "type" int
let nfile_next = field nfile "next" (ptr node)
let nfile_fd = field nfile "fd" int
let nfile_fname = field nfile "fname" (ptr node)
let nfile_expfname = field nfile "expfname" string
let () = seal nfile

let node_nfile = field node "nfile" nfile

type ndup

let ndup : ndup structure typ = structure "ndup"
let ndup_type = field ndup "type" int
let ndup_next = field ndup "next" (ptr node)
let ndup_fd = field ndup "fd" int
let ndup_dupfd = field ndup "dupfd" int
let ndup_vname = field ndup "vname" (ptr node)
let () = seal ndup

let node_ndup = field node "ndup" ndup

type nhere

let nhere : nhere structure typ = structure "nhere"
let nhere_type = field nhere "type" int
let nhere_next = field nhere "next" (ptr node)
let nhere_fd = field nhere "fd" int
let nhere_doc = field nhere "doc" (ptr node)
let () = seal nhere

let node_nhere = field node "nhere" nhere

type nnot

let nnot : nnot structure typ = structure "nnot"
let nnot_type = field nnot "type" int
let nnot_com = field nnot "com" (ptr node)
let () = seal nnot

let node_nnot = field node "nnot" nnot
let () = seal node
              
let parsecmd_safe : int -> node union ptr =
  foreign "parsecmd_safe" (int @-> returning (ptr node))
          
let parse s =
  setinputstring s; (* TODO set stack mark? *)
  parsecmd_safe 0

let neof : node union ptr = foreign_value "tokpushback" node
let nerr : node union ptr = foreign_value "lasttoken" node

let addrof p = raw_address_of_ptr (to_voidp p)

let eqptr p1 p2 = addrof p1 = addrof p2
                                  
let nullptr (p : 'a ptr) = addrof p = Nativeint.zero

type parse_result = Done | Error | Null | Parsed of (node union ptr)

let parse_next ?interactive:(i=false) () =
  let n = parsecmd_safe (if i then 1 else 0) in
  if eqptr n neof
  then Done
  else if eqptr n nerr
  then
    begin
      begin
        match !root_stackmark with
        | None -> failwith "!!! missing root stackmark"
        | Some smark -> pop_stack smark
      end;
      Error
    end
  else if nullptr n
  then Null (* comment or blank line or error ... *)
  else Parsed n
            
let (@->) (s : ('b, 'c) structured ptr) (f : ('a, ('b, 'c) structured) field) =
  getf (!@ s) f

let rec arglist (n : narg structure) : (narg structure) list =
  let next = getf n narg_next in
  if nullptr next
  then [n] 
  else
    (assert (next @-> node_type = 15);
     n::arglist (next @-> node_narg))

let rec nodelist (n : nodelist structure ptr) : (node union ptr) list =
  if nullptr n
  then []
  else (n @-> nodelist_n)::nodelist (n @-> nodelist_next)
                  
let rec redirlist (n : node union ptr) =
  if nullptr n
  then []
  else
    let h = match n @-> node_type with
      (* NTO *)
      | 16 -> `File (1,">",n @-> node_nfile)
      (* NCLOBBER *)
      | 17 -> `File (1,">|",n @-> node_nfile)
      (* NFROM *)
      | 18 -> `File (0,"<",n @-> node_nfile)
      (* NFROMTO *)
      | 19 -> `File (0,"<>",n @-> node_nfile)
      (* NAPPEND *)
      | 20 -> `File (1,">>",n @-> node_nfile)
      (* NTOFD *)      
      | 21 -> `Dup (1,">&",n @-> node_ndup)
      (* NFROMFD *)              
      | 22 -> `Dup (0,"<&",n @-> node_ndup)
      (* NHERE quoted heredoc---no expansion)*)
      | 23 -> `Here (0,"<<",false,n @-> node_nhere)
      (* NXHERE unquoted heredoc (param/command/arith expansion) *)
      | 24 -> `Here (0,"<<",true,n @-> node_nhere)
      | nt -> failwith ("unexpected node_type in redirlist: " ^ string_of_int nt)
    in
    h :: redirlist (getf (n @-> node_nfile) nfile_next)

let rec caselist (n : node union ptr) =
  if nullptr n
  then []
  else    
    let n = n @-> node_nclist in
    assert (getf n nclist_type = 13); (* NCLIST *)
    (getf n nclist_pattern, getf n nclist_body)::caselist (getf n nclist_next)
                   
let explode s =
  let rec exp i l =
    if i < 0 then l else exp (i - 1) (s.[i] :: l) in
  exp (String.length s - 1) []

let implode l =
  let s = Bytes.create (List.length l) in
  let rec imp i l =
    match l with
    | []  -> ()
    | (c::l) -> (Bytes.set s i c; imp (i+1) l)
  in
  imp 0 l;
  Bytes.unsafe_to_string s
                   
let rec intercalate p ss =
  match ss with
  | [] -> ""
  | [s] -> s
  | s::ss -> s ^ p ^ intercalate p ss          

let lines = Str.split (Str.regexp "[\n\r]+")

let rec fresh_marker ls s =
  if List.mem s ls
  then fresh_marker ls (s ^ (String.sub s (String.length s - 1) 1))
  else s
                      
let rec split_at p xs =
  match xs with
  | [] -> ([],[])
  | x::xs ->
     if p x
     then ([],x::xs)
     else let (xs,ys) = split_at p xs in
          (x::xs, ys)

let string_of_vs = function
  | 0x1 -> (* VSNORMAL ${var} *) []
  | 0x2 -> (* VSMINUS ${var-text} *) ['-']
  | 0x3 -> (* VSPLUS ${var+text} *) ['+']
  | 0x4 -> (* VSQUESTION ${var?message} *) ['?']
  | 0x5 -> (* VSASSIGN ${var=text} *) ['=']
  | 0x6 -> (* VSTRIMRIGHT ${var%pattern} *) ['%']
  | 0x7 -> (* VSTRIMRIGHTMAX ${var%%pattern} *) ['%';'%']
  | 0x8 -> (* VSTRIMLEFT ${var#pattern} *) ['#']
  | 0x9 -> (* VSTRIMLEFTMAX ${var##pattern} *) ['#';'#']
  | vs -> failwith ("Unknown VSTYPE: " ^ string_of_int vs)
                   
let braces s = "{ " ^ s ^ " ; }"
let parens s = "( " ^ s ^ " )"
                  
let rec show (n : node union ptr) : string =
  match (n @-> node_type) with
  (* NCMD *)
  | 0  ->
     let n = n @-> node_ncmd in
     let raw_cmd = intercalate " " (List.map sharg (arglist (getf n ncmd_args @-> node_narg))) in
     let vars = if nullptr (getf n ncmd_assign) then "" else intercalate " " (List.map sharg (arglist (getf n ncmd_assign @-> node_narg))) ^ " " in
     vars ^ raw_cmd ^ shredir (getf n ncmd_redirect)
  (* NPIPE *)
  | 1  ->
     let n = n @-> node_npipe in
     let cmds = nodelist (getf n npipe_cmdlist) in
     intercalate " | " (List.map show cmds) ^ if (getf n npipe_backgnd) = 0 then "" else " &"
  (* NREDIR *)
  | 2  -> shnredir braces n 
  (* NBACKGND *)
  | 3  -> shnredir braces n ^ " &"
  (* NSUBSHELL *)
  | 4  -> shnredir parens n
  (* NAND *)
  | 5  -> shbinary "&&" (n @-> node_nbinary)
  (* NOR *)
  | 6  -> shbinary "||" (n @-> node_nbinary)
  (* NSEMI *)
  | 7  -> shbinary ";" (n @-> node_nbinary)
  (* NIF *)
  | 8  -> shif (n @-> node_nif)
  (* NWHILE *)
  | 9  ->
     let n = n @-> node_nbinary in
     "while " ^ show (getf n nbinary_ch1) ^ "; do " ^ show (getf n nbinary_ch2) ^ "; done"
  (* NUNTIL *)
  | 10 ->
     let n = n @-> node_nbinary in
     "until " ^ show (getf n nbinary_ch1) ^ "; do " ^ show (getf n nbinary_ch2) ^ "; done"
  (* NFOR *)
  | 11 ->
     let n = n @-> node_nfor in
     "for " ^ (getf n nfor_var) ^ " in " ^ sharg (getf n nfor_args @-> node_narg) ^ "; do " ^ show (getf n nfor_body) ^ "; done"
  (* NCASE *)
  | 12 ->
     let n = n @-> node_ncase in
     "case " ^ sharg (getf n ncase_expr @-> node_narg) ^ " in " ^ shclist (getf n ncase_cases) ^ " esac"
  (* NDEFUN *)
  | 14 ->
     let n = n @-> node_ndefun in
     (getf n ndefun_text) ^ "() " ^ braces (show (getf n ndefun_body))
  (* NARG *)
  | 15 -> failwith "Didn't expect narg at the top-level"
  (* NNOT *)
  | 25 -> "! { " ^ show (getf (n @-> node_nnot) nnot_com) ^ " }"
  | nt -> failwith ("unexpected node_type " ^ string_of_int nt)

and shbinary (op : string) (n : nbinary structure) : string =
  show (getf n nbinary_ch1) ^ " " ^ op ^ " " ^ show (getf n nbinary_ch2)

and shnredir parenthesize n =
  let nr = n @-> node_nredir in
  parenthesize (show (getf nr nredir_n)) ^ shredir (getf nr nredir_redirect)

and shif n =
  "if " ^ show (getf n nif_test) ^
  "; then " ^ show (getf n nif_ifpart) ^
  (let else_part = getf n nif_elsepart in
   if nullptr else_part
   then "; fi"
   else if (else_part @-> node_type = 8)
   then "; el" ^ shif (else_part @-> node_nif)
   else "; else " ^ show else_part ^ "; fi")

and shclist clist = intercalate " " (List.map shcase (caselist clist)) (* handles NCLIST = 13 *)
    
and shcase (pat,body) =
  assert (pat @-> node_type = 15);
  sharg (pat @-> node_narg) ^ ") " ^ show body ^ ";;"
    
and shredir (n : node union ptr) : string =
  let redirs = redirlist n in
  if redirs = []
  then ""
  else " " ^ intercalate " " (List.map show_redir redirs)
and show_redir n : string =
  match n with
  | `File (src,sym,f) -> show_redir_src (getf f nfile_fd) src ^ sym ^ sharg ((getf f nfile_fname) @-> node_narg)
  | `Dup (src,sym,d) -> 
      let vname = getf d ndup_vname in
      let tgt =
        if nullptr vname
        then string_of_int (getf d ndup_dupfd)
        else sharg (vname @-> node_narg)
      in
     show_redir_src (getf d ndup_fd) src ^ sym ^ tgt
  | `Here (src,sym,exp,h) ->
     let heredoc = sharg ((getf h nhere_doc) @-> node_narg) in
     let marker = fresh_marker (lines heredoc) "EOF" in
     show_redir_src (getf h nhere_fd) src ^ sym ^ (if exp then marker else "'" ^ marker ^ "'") ^ "\n" ^ heredoc ^ marker
and show_redir_src actual expected =
  if actual = expected
  then ""
  else string_of_int actual
                                                    
and sharg (n : narg structure) : string =
  let str,s',bqlist,stack = show_arg (explode (getf n narg_text)) (getf n narg_backquote) [] in
  (* we should have used up the string and have no backquotes left in our list *)
  assert (s' = []);
  assert (nullptr bqlist);
  assert (stack = []);
  str    
and show_arg (s : char list) (bqlist : nodelist structure ptr) stack =
  (* we have to look at the string and interpret control characters... *)
  match s,stack with
  | [],[] -> "",[],bqlist,[]
  | [],`CTLVar::stack' -> failwith "End of string before CTLENDVAR"
  | [],`CTLAri::stack' -> failwith "End of string before CTLENDARI"
  | [],`CTLQuo::stack' -> failwith "End of string before CTLQUOTEMARK"
  (* CTLESC *)
  | '\129'::c::s',_ -> 
     let str,s'',bqlist',stack' = show_arg s' bqlist stack in
     let c' = match c with
      | '\'' -> "\\'"
      | '\"' -> "\\\""
      | _ -> String.make 1 c
     in
     c' ^ str,s'',bqlist',stack'
  (* CTLVAR *)
  | '\130'::t::s',_ -> 
     let v,s'',bqlist',stack' = show_var (int_of_char t) s' bqlist stack in
     assert (stack = stack');
     let str,s''',bqlist'',stack'' = show_arg s'' bqlist' stack' in
     "${" ^ v ^ "}" ^ str, s''', bqlist'', stack''
  (* CTLENDVAR *)
  | '\131'::s',`CTLVar::stack' -> "",[],bqlist,stack' (* s' gets handled by CTLVAR *)
  | '\131'::s',`CTLAri::stack' -> failwith "Saw CTLENDVAR before CTLENDARI"
  | '\131'::s',`CTLQuo::stack' -> failwith "Saw CTLENDVAR before CTLQUOTEMARK"
  | '\131'::s',[] -> failwith "Saw CTLENDVAR outside of CTLVAR"
  (* CTLBACKQ *)
  | '\132'::s',_ ->
     if nullptr bqlist
     then failwith "Saw CTLBACKQ but bqlist was null"
     else
       let n = bqlist @-> nodelist_n in
       (* MMG: !!! dash has a bug in its sharg function... it doesn't advance the list! *)
       let bqlist' = bqlist @-> nodelist_next in
       let str,s'',bqlist'',stack' = show_arg s' bqlist' stack in
       "$(" ^ show n ^ ")" ^ str,s'',bqlist'',stack'
  (* CTLARI *)
  | '\134'::s',_ ->
     let ari,s'',bqlist',stack' = show_arg s' bqlist (`CTLAri::stack) in
     assert (stack = stack');
     let str,s''',bqlist'',stack'' = show_arg s'' bqlist' stack' in
     "$((" ^ ari ^ "))" ^ str, s''', bqlist'', stack''
  (* CTLENDARI *)
  | '\135'::s',`CTLAri::stack' -> "",s',bqlist,stack'
  | '\135'::s',`CTLVar::stack' -> failwith "Saw CTLENDARI before CTLENDVAR"
  | '\135'::s',`CTLQuo::stack' -> failwith "Saw CTLENDARI before CTLQUOTEMARK"
  | '\135'::s',[] -> failwith "Saw CTLENDARI outside of CTLARI"
  (* CTLQUOTEMARK *)
  | '\136'::s',[`CTLQuo] -> "",s',bqlist,[]
  | '\136'::s',_ ->
     let quoted,s'',bqlist',stack' = show_arg  s' bqlist [`CTLQuo] in
     assert (stack' = []);
     let str,s''',bqlist'',stack'' = show_arg s'' bqlist' stack in
     "\"" ^ quoted ^ "\"" ^ str, s''', bqlist'', stack''
  (* ordinary character *)
  | c::s',_ -> 
     let str,s',bqlist',stack' = show_arg s' bqlist stack in
     let c' = match c with
      | '\'' -> "\\'"
      | '\"' -> "\\\""
      | _ -> String.make 1 c
     in
     c' ^ str,s',bqlist',stack'
and show_var (t : int) (s : char list) (bqlist : nodelist structure ptr) stack =
  let var_name,s' = split_at (fun c -> c = '=') s in
  (* mask out VSNUL, check VSTYPE *)
  match t land 0x0f, s' with
  (* VSNORMAL and VSLENGTH get special treatment

     neither ever gets VSNUL
     VSNORMAL is terminated just with the =, without a CTLENDVAR *)
  (* VSNORMAL *)
  | 0x1,'='::s'' -> implode var_name, s'', bqlist, stack
  (* VSLENGTH *)
  | 0xa,'='::'\131'::s'' -> implode (['#'] @ var_name), s'', bqlist, stack
  | 0x1,c::_ | 0xa,c::_ -> failwith ("Missing CTLENDVAR for VSNORMAL/VSLENGTH, found " ^ Char.escaped c)
  (* every other VSTYPE takes mods before CTLENDVAR *)
  | vstype,'='::s' ->
     (* check VSNUL *)
     let vsnul = if t land 0x10 = 1 then [] else [':'] in
     let mods,s'',bqlist',stack' = show_arg s' bqlist (`CTLVar::stack) in
     implode (var_name @ vsnul @ string_of_vs vstype) ^ mods, s'', bqlist', stack'
  | _,c::s' -> failwith ("Expected '=' terminating variable name, found " ^ Char.escaped c)
  | _,[] -> failwith "Expected '=' terminating variable name, found EOF"

